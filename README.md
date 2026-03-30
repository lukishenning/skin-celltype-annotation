# Skin Cell Type Annotation Pipeline

Automated cell type annotation pipeline for GA (Gingival/Autoimmune) vs Healthy skin scRNA-seq data using the Reynolds et al. 2021 (Science) human skin cell atlas reference.

## Overview

This pipeline performs reference-based cell type annotation on human skin scRNA-seq datasets using the comprehensive Reynolds 2021 skin cell atlas (528,253 cells, 34 cell types).

## Features

- Reference-based annotation using label transfer
- Memory-efficient workflow for large datasets
- Comprehensive visualizations (UMAP, composition, markers)
- Support for both GA (diseased) and Healthy skin comparisons
- Integration with Seurat v5
- Subset analysis for specific cell populations (e.g., APCs)

## Requirements

### R Packages
```r
install.packages(c("Seurat", "hdf5r", "anndata", "dplyr", "ggplot2"))
```

### Reference Data
Download the Reynolds 2021 healthy skin reference from:
- **Zenodo**: https://zenodo.org/record/4569496 (healthy.h5ad, ~1.2GB)
- **ArrayExpress**: E-MTAB-8142

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/lukishenning/skin-celltype-annotation.git
cd skin-celltype-annotation
```

### 2. Download reference data
```bash
mkdir -p data/reference
# Download healthy.h5ad from Zenodo
wget https://zenodo.org/record/4569496/files/healthy.h5ad -O data/reference/healthy.h5ad
```

### 3. Prepare your query data
Ensure your query Seurat object has:
- Normalized data in the `RNA` assay
- A `condition` column in metadata (e.g., "GA_skin", "Healthy_skin")
- Cells as columns, genes as rows

```r
library(Seurat)
query <- readRDS("your_integrated_object.rds")
query <- JoinLayers(query)  # Required for Seurat v5
```

### 4. Run annotation
```r
source("R/run_annotation.R")

# Set paths
work_dir <- "/path/to/skin-celltype-annotation"
data_dir <- "/path/to/your/data"
ref_dir <- file.path(work_dir, "data/reference")

# Run the complete pipeline
run_reference_annotation(
  query_file = file.path(data_dir, "your_query.rds"),
  reference_file = file.path(ref_dir, "healthy.h5ad"),
  output_dir = data_dir
)
```

## Workflows

### Workflow 1: Reference-Based Cell Type Annotation

The main pipeline consists of these steps:

1. **Load Reference**: Read h5ad file and convert to Seurat
2. **Prepare Reference**: Normalize, find variable features, scale, PCA, UMAP
3. **Load Query**: Load your integrated scRNA-seq data
4. **Find Anchors**: Identify integration anchors between reference and query
5. **Transfer Labels**: Transfer cell type labels from reference to query
6. **Visualize**: Generate UMAPs, composition plots, marker validation
7. **Save Results**: Export annotated object and summary statistics

### Workflow 2: APC (Antigen-Presenting Cell) Subset Analysis

After the main annotation, you can perform focused analysis on specific cell populations:

```r
source("R/APC_analysis.R")
```

This script:
1. Subsets APCs (DC1, DC2, LC, MigDC, Macro_1, Macro_2, Inf_mac, Mono_mac, moDC)
2. Re-normalizes and finds variable features
3. Re-clusters at multiple resolutions (0.2, 0.4, 0.6, 0.8)
4. Generates UMAP visualizations
5. Identifies cluster markers
6. Compares GA vs Healthy composition
7. Validates APC-specific markers

## Output Files

### Main Annotation
| File | Description |
|------|-------------|
| `seurat_annotated_reference.rds` | Query object with predicted cell types |
| `reference_celltypes.pdf` | UMAP of reference cell types |
| `reference_annotation_results.pdf` | Query UMAP with predictions |
| `reference_annotation_by_condition.pdf` | Split UMAP by condition |
| `reference_celltype_composition.pdf` | Cell type proportions |
| `reference_celltype_composition.csv` | Raw composition data |
| `reference_dotplot.pdf` | Marker validation dotplot |

### APC Analysis
| File | Description |
|------|-------------|
| `APC_seurat.rds` | APC subset Seurat object |
| `APC_umap_clusters.pdf` | UMAP with clusters, labels, and condition |
| `APC_umap_by_condition.pdf` | UMAP split by condition |
| `APC_cluster_markers.csv` | Top markers per cluster |
| `APC_celltype_by_cluster.csv` | Cell type composition per cluster |
| `APC_celltype_composition.pdf` | Bar plot of composition |
| `APC_marker_expression.pdf` | Feature plots of APC markers |
| `APC_dotplot.pdf` | Marker expression dotplot |
| `APC_GA_vs_Healthy.pdf` | GA vs Healthy comparison |
| `APC_GA_vs_Healthy_counts.csv` | Raw cell counts |
| `APC_GA_vs_Healthy_percentage.csv` | Percentages |

## Cell Types Annotated

The pipeline identifies these cell types from Reynolds 2021:

| Category | Cell Types |
|----------|------------|
| Keratinocytes | Differentiated_KC, Undifferentiated_KC |
| Fibroblasts | F1, F2, F3 |
| T cells | Th, Tc, Treg |
| Macrophages | Macro_1, Macro_2, Inf_mac, Mono_mac |
| Dendritic Cells | DC1, DC2, LC, migLC, MigDC, moDC |
| Endothelial | VE1, VE2, VE3, LE1, LE2 |
| Other | Melanocyte, Mast_cell, Schwann_1, Schwann_2, NK, ILC1_3, ILC1_NK, ILC2, Plasma, Pericyte_1, Pericyte_2 |

## Example Results

### Overall Composition by Condition
GA skin shows:
- Expansion of inflammatory fibroblasts (F2) and macrophages (Macro_1)
- Increased T cell infiltration (Tc: 8x, Treg: 3x)
- Plasma cell expansion (40x higher)
- Decreased melanocytes and ILC populations

### APC Composition by Condition
- **Macro_1** dramatically expanded in GA (53.5% vs 8.4%)
- **DC2** cells reduced in GA (7% vs 35.8%)
- **moDC** enriched in healthy skin

### APC Cluster Markers
| Cluster | Key Markers | Likely Identity |
|---------|-------------|----------------|
| C5 | CD1C, FCER1A, CLEC10A | cDC2 |
| C11 | XCR1, CLEC9A | cDC1 |
| C13 | LAMP3, IDO2 | mregDC/MigDC |
| C18 | CD207, EPCAM | Langerhans cells |
| C9 | CCL18, CCL13, CD209 | M2-like macrophages |

### Prediction Confidence
- Mean score: 0.80
- Median score: 0.865
- Low confidence (<0.5): ~10% of cells

## Troubleshooting

### Memory Issues
If you encounter memory errors:
- Reduce the number of variable features (`nfeatures = 1000`)
- Process samples in batches
- Use a machine with more RAM (recommended: 64GB+)

### SCTransform Errors
If SCTransform fails, the pipeline automatically falls back to standard LogNormalize normalization.

### Low Prediction Confidence
Cells with low confidence scores may represent:
- Rare cell types not well-represented in reference
- Novel cell states specific to your condition
- Poor quality cells requiring QC filtering

## Citation

If you use this pipeline, please cite:

```bibtex
Reynolds G, Vegh P, et al. (2021) Developmental cell programs are 
co-opted in inflammatory skin disease. Science 371(6527)
```

## License

MIT License

## Author

Lukas Henning
