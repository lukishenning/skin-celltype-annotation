# Example usage of the skin-celltype-annotation pipeline

# Load required libraries
library(Seurat)
library(anndata)

# Source the annotation function
source("R/run_annotation.R")

# Define paths
work_dir <- "/path/to/skin-celltype-annotation"
data_dir <- "/path/to/your/data"
ref_dir <- file.path(work_dir, "data/reference")
output_dir <- file.path(data_dir, "annotation_results")

# Create output directory
dir.create(output_dir, showWarnings = FALSE)

# Download reference data (if not already present)
ref_file <- file.path(ref_dir, "healthy.h5ad")
if (!file.exists(ref_file)) {
  dir.create(ref_dir, recursive = TRUE)
  download.file(
    "https://zenodo.org/record/4569496/files/healthy.h5ad",
    destfile = ref_file,
    method = "curl"
  )
}

# Run annotation
annotated <- run_reference_annotation(
  query_file = file.path(data_dir, "seurat_integrated.rds"),
  reference_file = ref_file,
  output_dir = output_dir,
  n_hvg = 3000,
  n_pcs = 30
)

# Analyze results
cat("\nPredicted cell types:\n")
print(table(annotated$predicted.id))

cat("\nComposition by condition:\n")
print(table(annotated$predicted.id, annotated$condition))

cat("\nPrediction confidence:\n")
print(summary(annotated$prediction.score.max))
