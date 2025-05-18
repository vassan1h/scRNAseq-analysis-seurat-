# # scRNA-seq Analysis of GSE197177 using seurat 
**Author:** Vassanth Mathan  
**Date:** April 2025  

## ⚙️ Setup & Dependencies

This pipeline is built in **R (≥ 4.3)** and depends heavily on the Seurat ecosystem and several additional packages. This workflow analyzes single-cell RNA-seq data from [GSE197177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197177), focusing on epithelial and immune populations across eight samples. It implements per-sample and integrated workflows using Seurat, Harmony, Slingshot, and SingleR to derive cluster structure, pseudotemporal trajectories, marker gene signatures, and proportional analysis.

### Required R Packages:
Install all required packages using:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyverse", "ggrepel", "pheatmap", "ggpubr", 
  "reshape2", "RColorBrewer"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Seurat", "SingleCellExperiment", "slingshot", "SingleR", 
  "celldex", "scDblFinder", "EnhancedVolcano"
))

# Harmony (from GitHub)
devtools::install_github("immunogenomics/harmony")
```

## Workflow Summary  

### 1. **Input & QC**
- 8 samples loaded using `Read10X()`, converted to Seurat objects.
- QC metrics computed (`nFeature_RNA`, `nCount_RNA`, `percent.mt`).
- Violin plots and scatter plots for initial inspection.
- Cells filtered (>200 genes, <20% mitochondrial reads).

### 2. **Normalization & Feature Selection**
- Each sample normalized using `LogNormalize`.
- Top 2000 variable features selected via VST method.

### 3. **Dimensionality Reduction & Clustering**
- PCA and clustering (resolution = 0.5) performed **per sample**.
- UMAPs visualized for each object.
- Clustering results summarized.

### 4. **Doublet Detection**
- `scDblFinder` run on SCE-converted Seurat objects.
- Summary table generated per sample.

### 5. **Integration (Harmony)**
- Normalized samples merged into a single object.
- Scaled, PCA run on merged object.
- Batch correction performed using `RunHarmony()`.
- Post-Harmony clustering and UMAP performed on first 25 Harmony dimensions.

### 6. **Pseudotime Analysis**
- Seurat object converted to `SingleCellExperiment`.
- Slingshot run using UMAP reduction and cluster labels.
- Pseudotime plotted over UMAP.

### 7. **Marker Identification & Annotation**
- `FindAllMarkers()` run on sample 1.
- Top markers visualized via DotPlot and FeaturePlot.
- Automated cell type annotation performed using `SingleR` with Blueprint+ENCODE reference.

### 8. **Cell Type Proportion**
- Cluster frequencies per sample calculated and plotted as bar charts.

### 9. **Reproducibility**
- Full session information logged via:
  ```r
  sink("sessionInfo.txt")
  sessionInfo()
  sink()
  ```
---
***Reference:***
Zhang, S., Fang, W., Zhou, S. et al. Single cell transcriptomic analyses implicate an immunosuppressive tumor microenvironment in pancreatic cancer liver metastasis. Nat Commun 14, 5123 (2023). https://doi.org/10.1038/s41467-023-40727-7
---
