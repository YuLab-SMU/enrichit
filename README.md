# enrichit: C++ Implementations of Functional Enrichment Analysis

The `enrichit` package provides fast, efficient, and lightweight implementations of common functional enrichment analysis methods, including **Over-Representation Analysis (ORA)** and **Gene Set Enrichment Analysis (GSEA)**. The core algorithms are implemented in C++ using `Rcpp` to ensure high performance, making it suitable for analyzing large datasets or running simulations.

## Installation

You can install the development version of `enrichit` from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("YuLab-SMU/enrichit")
```

## Features

- **High Performance**: Core calculations are written in C++.
- **ORA**: Standard hypergeometric test for over-representation analysis.
- **GSEA**: 
  - **Multilevel**: Efficient p-value estimation for high-significance results (similar to `fgsea`).
  - **Permutation**: Standard permutation-based p-value calculation.
  - **Adaptive**: Adaptive permutation approach.
- **GSON Support**: Native support for `GSON` objects for gene set management.
- **Standardized Output**: Returns `enrichResult` and `gseaResult` objects compatible with the `clusterProfiler` ecosystem.

## Usage

### Over-Representation Analysis (ORA)

```r
library(enrichit)

# Example gene sets
gene_sets <- list(
  pathway1 = paste0("Gene", 1:50),
  pathway2 = paste0("Gene", 51:100)
)

# Define a universe and a list of significant genes
universe <- paste0("Gene", 1:1000)
sig_genes <- paste0("Gene", 1:20) # Significant genes

# Run ORA
ora_res <- ora(gene = sig_genes, 
               gene_sets = gene_sets, 
               universe = universe)

print(ora_res)
```

### Gene Set Enrichment Analysis (GSEA)

```r
library(enrichit)

# Generate a ranked gene list
set.seed(123)
geneList <- sort(rnorm(1000), decreasing = TRUE)
names(geneList) <- paste0("Gene", 1:1000)

# Define gene sets
gene_sets <- list(
  pathway1 = paste0("Gene", 1:50),  # Enriched at top
  pathway2 = paste0("Gene", 951:1000) # Enriched at bottom
)

# Run GSEA
gsea_res <- gsea(geneList = geneList, 
                 gene_sets = gene_sets, 
                 method = "multilevel")

print(gsea_res)
```
