# enrichit: C++ Implementations of Functional Enrichment Analysis

`enrichit` is part of the `clusterProfiler` family, serving as the underlying algorithm implementation layer. It focuses on fast core computation, standardized result objects, and reusable data-preparation layers for downstream visualization packages such as `enrichplot`.

The package now covers not only classical enrichment workflows such as **ORA** and **GSEA**, but also **weighted enrichment**, **network propagation-based enrichment**, **multi-omics early/late fusion**, **multi-layer topology fusion**, and **explanation-ready data extraction**.

## Installation

You can install the development version of `enrichit` from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("YuLab-SMU/enrichit")
```

## Scope

`enrichit` is designed around four layers:

- **Core enrichment engines**: ORA, GSEA, weighted ORA/GSEA, and GSON-aware variants.
- **Network-aware enrichment**: single-layer `nsea()` and multi-layer `mnsea()` workflows based on Random Walk with Restart.
- **Multi-omics integration**: early fusion at the feature level and late fusion at the pathway level.
- **Explanation-ready outputs**: contribution tables and topology-aware extraction helpers prepared for visualization in `enrichplot`.

## Feature Map

- **High performance core**: key algorithms are implemented in `C++` via `Rcpp`, with sparse network propagation powered by `RcppEigen`.
- **ORA**: standard hypergeometric ORA with optional weighted ORA through Wallenius' noncentral hypergeometric distribution.
- **GSEA**: multilevel, permutation, and adaptive strategies for ranked enrichment analysis.
- **GSON support**: native `ora_gson()` and `gsea_gson()` interfaces for structured gene set collections.
- **NSEA**: `nsea()` and `nsea_gson()` for network-ranked enrichment on a single graph, including `mode = "signed"` for bidirectional propagation.
- **Multi-layer topology fusion**: `mnsea()` and `mnsea_gson()` for multiplex or heterogeneous network propagation across multiple layers.
- **Multi-omics early fusion**: `aggregate_omics()`, `harmonize_ids()`, and `select_features_for_ora()` for feature-level integration before enrichment.
- **Multi-omics late fusion**: `aggregate_enrichment()` for pathway-level aggregation of multiple enrichment results.
- **Contribution tracing**: `get_omics_contribution()`, `classify_omics_pattern()`, and `get_mnsea_contribution()` for explanation-oriented summaries.
- **Topology-aware extraction**: `extract_mnsea_subnetwork()` for pathway-specific node/edge tables that can be passed to downstream visualization packages.
- **Bayesian compression**: `bayes_enrich()` and `bayes_summary()` for posterior-based term prioritization.

## Main APIs

### Classical enrichment

- `ora()`, `ora_gson()`
- `gsea()`, `gsea_gson()`
- `gseaScores()`

### Weighted enrichment

- `ora(..., weight = )`
- `ora_gson(..., weight = )`
- `gsea(..., weight = )`
- `gsea_gson(..., weight = )`

### Network-aware enrichment

- `prepare_network()`
- `nsea()`, `nsea_gson()`
- `prepare_multilayer_network()`
- `propagate_multilayer()`
- `collapse_multilayer_scores()`
- `mnsea()`, `mnsea_gson()`

### Multi-omics integration

- `aggregate_omics()`
- `harmonize_ids()`
- `select_features_for_ora()`
- `aggregate_enrichment()`

### Explanation helpers

- `get_omics_contribution()`
- `classify_omics_pattern()`
- `get_mnsea_contribution()`
- `extract_mnsea_subnetwork()`

## Result Objects

The package provides the standard enrichment result object model used by the `clusterProfiler` family:

- `enrichResult` for ORA-like workflows
- `gseaResult` for ranked enrichment workflows
- `nseaResult` for single-network propagation plus enrichment
- `mnseaResult` for multi-layer propagation, collapsed scores, and cached explanation tables

These objects are intended to support a clean separation of concerns across the `clusterProfiler` family:

- `enrichit` handles core computation, algorithm implementation, and explanation-ready data preparation
- `clusterProfiler` provides high-level biological interpretation workflows and general enrichment analysis interfaces
- `enrichplot` handles visualization
- `gson` provides a structured gene set resource layer for managing and exchanging gene set collections across the family
- knowledge-base-oriented downstream packages such as `DOSE`, `ReactomePA`, `meshes`, and `MicrobiomeProfiler` provide domain-specific annotation and interpretation layers

## Design Notes

- **Computation first**: this package prioritizes fast and robust numerical routines over plot helpers.
- **Decoupled architecture**: integration layers, propagation layers, and enrichment layers are exposed separately where useful.
- **Stable downstream interface**: explanation helpers return standard tables so that plotting logic can evolve independently in downstream packages.
