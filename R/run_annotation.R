#' @title Run Reference-Based Cell Type Annotation
#' @description Automated cell type annotation using Reynolds 2021 reference
#' @param query_file Path to query Seurat object (rds)
#' @param reference_file Path to reference h5ad file
#' @param output_dir Output directory for results
#' @param n_hvg Number of highly variable genes (default: 3000)
#' @param n_pcs Number of PCA dimensions (default: 30)
#' @param k_anchor Number of anchors (default: 20)
#' @param k_score Scoring neighbors (default: 30)
#' @return Annotated Seurat object with predicted cell types
#' @export
#' @import Seurat anndata dplyr ggplot2
run_reference_annotation <- function(
    query_file,
    reference_file,
    output_dir,
    n_hvg = 3000,
    n_pcs = 30,
    k_anchor = 20,
    k_score = 30
) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("===============================================\n")
  cat("Reference-Based Cell Type Annotation Pipeline\n")
  cat("Using Reynolds et al. 2021 Reference\n")
  cat("===============================================\n\n")
  
  set.seed(42)
  
  # ============================================================================
  # STEP 1: Load Reference
  # ============================================================================
  cat("Step 1: Loading reference data...\n")
  cat("File:", reference_file, "\n")
  
  h5ad_data <- anndata::read_h5ad(reference_file)
  cat("Reference dimensions:", dim(h5ad_data)[1], "genes x", dim(h5ad_data)[2], "cells\n")
  
  obs_df <- as.data.frame(h5ad_data$obs)
  cat("Cell types in reference:", length(unique(obs_df$full_clustering)), "\n")
  
  # Convert to matrix and create Seurat
  counts <- t(h5ad_data$X)
  if (is(counts, "scipy.sparse.csr.csr_matrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  rm(h5ad_data); gc()
  
  ref_seurat <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = obs_df,
    project = "Reynolds2021_Reference"
  )
  rm(counts, obs_df); gc()
  
  ref_seurat$celltype <- ref_seurat$full_clustering
  Seurat::Idents(ref_seurat) <- "celltype"
  
  cat("Reference created:", Seurat::ncol(ref_seurat), "cells,", Seurat::nrow(ref_seurat), "genes\n\n")
  
  # ============================================================================
  # STEP 2: Prepare Reference
  # ============================================================================
  cat("Step 2: Preparing reference...\n")
  
  Seurat::DefaultAssay(ref_seurat) <- "RNA"
  ref_seurat <- Seurat::NormalizeData(ref_seurat, verbose = FALSE)
  ref_seurat <- Seurat::FindVariableFeatures(ref_seurat, nfeatures = min(2000, nrow(ref_seurat)), verbose = FALSE)
  ref_seurat <- Seurat::ScaleData(ref_seurat, verbose = FALSE)
  ref_seurat <- Seurat::RunPCA(ref_seurat, npcs = n_pcs, verbose = FALSE)
  ref_seurat <- Seurat::RunUMAP(ref_seurat, dims = 1:n_pcs, verbose = FALSE)
  
  # Plot reference UMAP
  p_ref <- Seurat::DimPlot(ref_seurat, group.by = "celltype", label = TRUE, repel = TRUE, label.size = 2) + 
    ggplot2::theme(legend.text = ggplot2::element_text(size = 6)) +
    Seurat::NoLegend() + ggplot2::ggtitle("Reynolds 2021 Reference Cell Types")
  
  ggplot2::ggsave(
    file.path(output_dir, "reference_celltypes.pdf"), 
    plot = p_ref, 
    width = 14, height = 10
  )
  cat("Reference UMAP saved\n\n")
  
  # ============================================================================
  # STEP 3: Load Query
  # ============================================================================
  cat("Step 3: Loading query data...\n")
  cat("File:", query_file, "\n")
  
  query_seurat <- readRDS(query_file)
  
  # Handle Seurat v5 multi-layer format
  if (Seurat::IsSeuratV5(query_seurat) || "SCT" %in% names(query_seurat@assays)) {
    if (length(query_seurat@assays$SCT@layers) > 1) {
      query_seurat <- Seurat::JoinLayers(query_seurat)
    }
  } else if (length(query_seurat@assays$RNA@layers) > 1) {
    query_seurat <- Seurat::JoinLayers(query_seurat)
  }
  
  cat("Query loaded:", Seurat::ncol(query_seurat), "cells\n\n")
  
  # ============================================================================
  # STEP 4: Prepare Query
  # ============================================================================
  cat("Step 4: Preparing query...\n")
  
  Seurat::DefaultAssay(query_seurat) <- "RNA"
  query_seurat <- Seurat::NormalizeData(query_seurat, verbose = FALSE)
  query_seurat <- Seurat::FindVariableFeatures(query_seurat, nfeatures = n_hvg, verbose = FALSE)
  query_seurat <- Seurat::ScaleData(query_seurat, verbose = FALSE)
  query_seurat <- Seurat::RunPCA(query_seurat, npcs = n_pcs, verbose = FALSE)
  
  cat("Query PCA completed\n\n")
  
  # ============================================================================
  # STEP 5: Find Integration Anchors
  # ============================================================================
  cat("Step 5: Finding integration anchors...\n")
  
  anchors <- Seurat::FindTransferAnchors(
    reference = ref_seurat,
    query = query_seurat,
    normalization.method = "LogNormalize",
    reference.assay = "RNA",
    query.assay = "RNA",
    dims = 1:n_pcs,
    k.filter = NA,
    k.anchor = k_anchor,
    k.score = k_score,
    verbose = TRUE
  )
  
  cat("Found", length(anchors@anchors), "anchors\n\n")
  
  # ============================================================================
  # STEP 6: Transfer Labels
  # ============================================================================
  cat("Step 6: Transferring cell type labels...\n")
  
  predictions <- Seurat::TransferData(
    anchorset = anchors,
    refdata = ref_seurat$celltype,
    dims = 1:n_pcs,
    weight.reduction = query_seurat[["pca"]],
    sd.weight = 1,
    verbose = TRUE
  )
  
  query_seurat <- Seurat::AddMetaData(query_seurat, metadata = predictions)
  
  # ============================================================================
  # STEP 7: Calculate Prediction Scores
  # ============================================================================
  cat("\nStep 7: Calculating prediction statistics...\n")
  
  score_cols <- grep("prediction.score", colnames(query_seurat@meta.data), value = TRUE)
  score_cols <- score_cols[score_cols != "prediction.score.max"]
  max_scores <- apply(query_seurat@meta.data[, score_cols, drop = FALSE], 1, max, na.rm = TRUE)
  query_seurat$prediction.score.max <- max_scores
  
  cat("Mean prediction score:", round(mean(max_scores, na.rm = TRUE), 3), "\n")
  cat("Median prediction score:", round(median(max_scores, na.rm = TRUE), 3), "\n")
  cat("Low confidence (<0.5):", sum(max_scores < 0.5, na.rm = TRUE), "\n\n")
  
  # ============================================================================
  # STEP 8: Create UMAP
  # ============================================================================
  cat("Step 8: Creating UMAP for query...\n")
  query_seurat <- Seurat::RunUMAP(query_seurat, dims = 1:n_pcs, verbose = FALSE)
  cat("UMAP completed\n\n")
  
  # ============================================================================
  # STEP 9: Visualizations
  # ============================================================================
  cat("Step 9: Creating visualizations...\n")
  
  # Main annotation plot
  p1 <- Seurat::DimPlot(query_seurat, reduction = "umap", group.by = "predicted.id", 
                        label = TRUE, repel = TRUE, label.size = 2.5) + 
    Seurat::NoLegend() + ggplot2::ggtitle("Predicted Cell Types (Reynolds 2021)")
  
  p2 <- Seurat::DimPlot(query_seurat, reduction = "umap", group.by = "condition") + 
    ggplot2::ggtitle("Condition")
  
  p3 <- Seurat::FeaturePlot(query_seurat, features = "prediction.score.max", reduction = "umap", 
                            cols = c("lightgrey", "blue")) + ggplot2::ggtitle("Prediction Confidence")
  
  ggplot2::ggsave(
    file.path(output_dir, "reference_annotation_results.pdf"), 
    plot = p1 | p2 | p3, 
    width = 18, height = 6
  )
  
  # By condition
  p4 <- Seurat::DimPlot(query_seurat, reduction = "umap", split.by = "condition",
                        group.by = "predicted.id", label = TRUE, repel = TRUE, label.size = 2) + 
    Seurat::NoLegend() + ggplot2::ggtitle("Cell Types by Condition")
  
  ggplot2::ggsave(
    file.path(output_dir, "reference_annotation_by_condition.pdf"), 
    plot = p4, width = 16, height = 6
  )
  
  # Violin plot
  p5 <- Seurat::VlnPlot(query_seurat, features = "prediction.score.max", group.by = "predicted.id", 
                        pt.size = 0) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "Prediction Score by Cell Type", y = "Max Prediction Score") +
    ggplot2::ylim(0, 1.05)
  
  ggplot2::ggsave(
    file.path(output_dir, "reference_prediction_violin.pdf"), 
    plot = p5, width = 14, height = 6
  )
  
  cat("Visualizations saved\n\n")
  
  # ============================================================================
  # STEP 10: Composition Analysis
  # ============================================================================
  cat("Step 10: Composition analysis...\n")
  
  if ("condition" %in% colnames(query_seurat@meta.data)) {
    prop_table <- table(query_seurat$predicted.id, query_seurat$condition)
    
    prop_df <- as.data.frame(prop_table)
    colnames(prop_df) <- c("CellType", "Condition", "Count")
    
    p6 <- ggplot2::ggplot(prop_df, ggplot2::aes(x = CellType, y = Count, fill = Condition)) + 
      ggplot2::geom_bar(stat = "identity", position = "fill") + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "Cell Type Proportion by Condition", y = "Proportion", x = "")
    
    p7 <- ggplot2::ggplot(prop_df, ggplot2::aes(x = CellType, y = Count, fill = Condition)) + 
      ggplot2::geom_bar(stat = "identity", position = "stack") + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "Cell Type Count by Condition", y = "Cell Count", x = "")
    
    ggplot2::ggsave(
      file.path(output_dir, "reference_celltype_composition.pdf"), 
      plot = p6 / p7, width = 14, height = 8
    )
    
    write.csv(prop_table, file.path(output_dir, "reference_celltype_composition.csv"))
    cat("Composition analysis saved\n\n")
  }
  
  # ============================================================================
  # STEP 11: Marker Validation
  # ============================================================================
  cat("Step 11: Marker validation...\n")
  
  marker_genes <- c(
    "KRT5", "KRT14", "KRT1", "KRT10",           # Keratinocytes
    "COL1A1", "COL1A2", "PDGFRA",               # Fibroblasts
    "VWF", "PECAM1", "FLT1", "FLT4", "PROX1",   # Endothelial
    "CD8A", "CD8B", "CD4", "FOXP3",              # T cells
    "CD68", "CD163", "CD14",                      # Macrophages
    "CD1C", "FCER1A", "CD207",                   # DCs
    "PMEL", "MITF", "SOX10",                      # Melanocytes/Schwann
    "TPSAB1", "KIT", "MS4A1", "JCHAIN"           # Mast/B cells
  )
  marker_genes <- marker_genes[marker_genes %in% rownames(query_seurat)]
  
  p8 <- Seurat::DotPlot(query_seurat, features = marker_genes, group.by = "predicted.id") + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
                   axis.text.y = ggplot2::element_text(size = 9))
  
  ggplot2::ggsave(
    file.path(output_dir, "reference_dotplot.pdf"), 
    plot = p8, width = 16, height = 10
  )
  
  cat("Marker validation saved\n\n")
  
  # ============================================================================
  # STEP 12: Save Results
  # ============================================================================
  cat("Step 12: Saving results...\n")
  
  saveRDS(query_seurat, file.path(output_dir, "seurat_annotated_reference.rds"))
  saveRDS(ref_seurat, file.path(output_dir, "reference_prepared.rds"))
  
  cat("Objects saved:\n")
  cat(" - seurat_annotated_reference.rds\n")
  cat(" - reference_prepared.rds\n\n")
  
  # ============================================================================
  # Summary
  # ============================================================================
  cat("===============================================\n")
  cat("Annotation Complete!\n")
  cat("===============================================\n")
  cat("\nSummary:\n")
  cat("- Reference:", Seurat::ncol(ref_seurat), "cells,", length(unique(ref_seurat$celltype)), "cell types\n")
  cat("- Query:", Seurat::ncol(query_seurat), "cells\n")
  cat("- Predicted:", length(unique(query_seurat$predicted.id)), "cell types\n")
  cat("- Mean confidence:", round(mean(max_scores, na.rm = TRUE), 3), "\n")
  cat("\nOutput directory:", output_dir, "\n")
  
  return(query_seurat)
}
