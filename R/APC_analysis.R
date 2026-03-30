#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(ggplot2)

cat("===============================================\n")
cat("Antigen-Presenting Cell (APC) Subset Analysis\n")
cat("===============================================\n\n")

set.seed(42)

work_dir <- "/Users/Henning/OpenMind Folders/Project_2/GSM files"
data_dir <- file.path(work_dir, "QC_Results_v2")
out_dir <- file.path(data_dir, "APC_analysis")

if (!dir.exists(out_dir)) dir.create(out_dir)

cat("Step 1: Load annotated object\n")
cat("===============================================\n")

obj <- readRDS(file.path(data_dir, "seurat_annotated_reference.rds"))
obj <- JoinLayers(obj)
cat("Full object:", ncol(obj), "cells\n")

apc_types <- c("LC", "DC1", "DC2", "MigDC", "Macro_1", "Macro_2", "Inf_mac", "Mono_mac", "moDC")

apc_obj <- subset(obj, cells = WhichCells(obj, expression = predicted.id %in% apc_types))

cat("\nAPC subset:", ncol(apc_obj), "cells\n")
print(table(apc_obj$predicted.id, apc_obj$condition))

cat("\nStep 2: Re-normalize and find variable features\n")
cat("===============================================\n")

DefaultAssay(apc_obj) <- "RNA"
apc_obj <- NormalizeData(apc_obj, verbose = FALSE)
apc_obj <- FindVariableFeatures(apc_obj, nfeatures = 2000, verbose = FALSE)

hvg <- VariableFeatures(apc_obj)
cat("HVGs identified:", length(hvg), "\n")

top_hvg <- head(hvg, 10)
cat("\nTop 10 HVGs:\n")
cat(paste(top_hvg, collapse = ", "), "\n")

p_hvg <- VariableFeaturePlot(apc_obj)
p_hvg <- LabelPoints(plot = p_hvg, points = top_hvg, repel = TRUE, xnudge = 0, ynudge = 0)
pdf(file.path(out_dir, "APC_hvg.pdf"), width = 10, height = 6)
print(p_hvg)
dev.off()

cat("\nStep 3: Scale and run PCA\n")
cat("===============================================\n")

apc_obj <- ScaleData(apc_obj, features = hvg, verbose = FALSE)
apc_obj <- RunPCA(apc_obj, features = hvg, npcs = 30, verbose = FALSE)

pdf(file.path(out_dir, "APC_elbow.pdf"), width = 8, height = 5)
print(ElbowPlot(apc_obj, ndims = 30))
dev.off()

cat("\nStep 4: Find neighbors and cluster\n")
cat("===============================================\n")

apc_obj <- FindNeighbors(apc_obj, reduction = "pca", dims = 1:30, verbose = FALSE)

for (res in c(0.2, 0.4, 0.6, 0.8)) {
  apc_obj <- FindClusters(apc_obj, resolution = res, algorithm = 4, verbose = FALSE)
  n_clust <- length(unique(apc_obj@meta.data[[paste0("RNA_snn_res.", res)]]))
  cat("Resolution", res, "->", n_clust, "clusters\n")
}

res_use <- 0.6
apc_obj$apc_clusters <- apc_obj@meta.data[[paste0("RNA_snn_res.", res_use)]]
n_clusters <- length(unique(apc_obj$apc_clusters))

apc_obj$apc_clusters <- factor(paste0("C", apc_obj$apc_clusters))

cat("\nUsing resolution", res_use, "with", n_clusters, "clusters\n")

cat("\nStep 5: Run UMAP\n")
cat("===============================================\n")

apc_obj <- RunUMAP(apc_obj, reduction = "pca", dims = 1:30, verbose = FALSE)

cat("\nStep 6: Visualizations\n")
cat("===============================================\n")

p1 <- DimPlot(apc_obj, reduction = "umap", group.by = "apc_clusters", 
              label = TRUE, repel = TRUE) + NoLegend() + 
  ggtitle(paste0("APC Clusters (n=", n_clusters, ")"))

p2 <- DimPlot(apc_obj, reduction = "umap", group.by = "predicted.id", 
              label = TRUE, repel = TRUE) + NoLegend() + 
  ggtitle("Original Cell Type Labels")

p3 <- DimPlot(apc_obj, reduction = "umap", group.by = "condition") + 
  ggtitle("Condition")

pdf(file.path(out_dir, "APC_umap_clusters.pdf"), width = 16, height = 5)
print(p1 | p2 | p3)
dev.off()

p4 <- DimPlot(apc_obj, reduction = "umap", split.by = "condition",
              group.by = "apc_clusters", label = TRUE, label.size = 3) + 
  NoLegend() + ggtitle("APC Clusters by Condition")

pdf(file.path(out_dir, "APC_umap_by_condition.pdf"), width = 14, height = 5)
print(p4)
dev.off()

cat("\nStep 7: Find cluster markers\n")
cat("===============================================\n")

Idents(apc_obj) <- "apc_clusters"

markers <- FindAllMarkers(apc_obj, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.5, verbose = FALSE)

markers <- markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC)

cat("\nTop markers per cluster:\n")
for (i in sort(unique(markers$cluster))) {
  cat("\n--- Cluster", i, "---\n")
  cat(paste(head(markers$gene[markers$cluster == i], 8), collapse = ", "), "\n")
}

write.csv(markers, file.path(out_dir, "APC_cluster_markers.csv"))

cat("\nStep 8: Cell type composition of clusters\n")
cat("===============================================\n")

ct_clust <- table(apc_obj$predicted.id, apc_obj$apc_clusters)
write.csv(ct_clust, file.path(out_dir, "APC_celltype_by_cluster.csv"))

cat("\nCell type composition per cluster:\n")
print(ct_clust)

p5 <- ggplot(as.data.frame(ct_clust), aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cluster", y = "Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cell Type Composition by Cluster")

pdf(file.path(out_dir, "APC_celltype_composition.pdf"), width = 12, height = 6)
print(p5)
dev.off()

cat("\nStep 9: APC marker validation\n")
cat("===============================================\n")

apc_markers <- c(
  "CD207", "EPCAM", "LY75",
  "CLEC9A", "XCR1", "IRF8",
  "CD1C", "FCER1A", "CLEC10A",
  "CCR7", "CD83", "LAMP3",
  "CD68", "CD163", "CD14",
  "FCN1", "VCAN", "S100A8",
  "TREM2", "APOE", "C1QC"
)
apc_markers <- apc_markers[apc_markers %in% rownames(apc_obj)]

p6 <- FeaturePlot(apc_obj, features = apc_markers, reduction = "umap", ncol = 4)

pdf(file.path(out_dir, "APC_marker_expression.pdf"), width = 14, height = 12)
print(p6)
dev.off()

p7 <- DotPlot(apc_obj, features = apc_markers, group.by = "apc_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9)) +
  labs(title = "APC Marker Expression by Cluster")

pdf(file.path(out_dir, "APC_dotplot.pdf"), width = 12, height = 8)
print(p7)
dev.off()

cat("\nStep 10: GA vs Healthy comparison\n")
cat("===============================================\n")

ga_healthy <- table(apc_obj$predicted.id, apc_obj$condition)
ga_healthy_pct <- prop.table(ga_healthy, margin = 2) * 100

cat("\nCell counts:\n")
print(ga_healthy)

cat("\nPercentages:\n")
print(round(ga_healthy_pct, 1))

p8 <- ggplot(as.data.frame(ga_healthy_pct), aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cell Type", y = "Percentage (%)", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("APC Composition: GA vs Healthy")

pdf(file.path(out_dir, "APC_GA_vs_Healthy.pdf"), width = 12, height = 6)
print(p8)
dev.off()

write.csv(ga_healthy, file.path(out_dir, "APC_GA_vs_Healthy_counts.csv"))
write.csv(ga_healthy_pct, file.path(out_dir, "APC_GA_vs_Healthy_percentage.csv"))

cat("\nStep 11: Save results\n")
cat("===============================================\n")

saveRDS(apc_obj, file.path(out_dir, "APC_seurat.rds"))

cat("\n===============================================\n")
cat("APC Analysis Complete!\n")
cat("===============================================\n")

cat("\nSummary:\n")
cat("- Total APCs:", ncol(apc_obj), "cells\n")
cat("- Clusters:", n_clusters, "at resolution", res_use, "\n")
cat("- Cell types:", length(unique(apc_obj$predicted.id)), "\n")
cat("\nOutput directory:", out_dir, "\n")
