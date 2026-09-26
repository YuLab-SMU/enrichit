# enrichit: C++ Implementations of Functional Enrichment Analysis

`enrichit` provides C++ implementations of functional enrichment analysis methods and S4 result classes used by the `clusterProfiler` family. It supports ORA, GSEA, weighted enrichment, network-based enrichment, multilayer network workflows, and multi-omics aggregation and contribution analysis.

## Installation

You can install the development version of `enrichit` from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("YuLab-SMU/enrichit")
```

## Main components

`enrichit` organizes its functions around four components:

- **Core enrichment engines**: ORA, GSEA, weighted ORA/GSEA, and GSON-aware variants.
- **Network-aware enrichment**: single-layer `nsea()` and multi-layer `mnsea()` workflows based on Random Walk with Restart.
- **Multi-omics integration**: early fusion at the feature level and late fusion at the pathway level.
- **Contribution and topology outputs**: contribution tables and topology-aware extraction helpers for visualization in `enrichplot`.

## Implementation and features

- **Implementation**: Algorithms are implemented in `C++` via `Rcpp`, with sparse network propagation using `RcppEigen`.
- **ORA**: standard hypergeometric ORA with optional weighted ORA through Wallenius' noncentral hypergeometric distribution.
- **GSEA**: multilevel, permutation, and adaptive strategies for ranked enrichment analysis.
- **GSON support**: native `ora_gson()` and `gsea_gson()` interfaces for structured gene set collections.
- **NSEA**: `nsea()` and `nsea_gson()` for network-ranked enrichment on a single graph, including `mode = "signed"` for bidirectional propagation.
- **Multi-layer topology fusion**: `mnsea()` and `mnsea_gson()` for multiplex or heterogeneous network propagation across multiple layers.
- **Multi-omics early fusion**: `aggregate_omics()`, `harmonize_ids()`, and `select_features_for_ora()` for feature-level integration before enrichment.
- **Multi-omics late fusion**: `aggregate_enrichment()` for pathway-level aggregation of multiple enrichment results.
- **Contribution tracing**: `get_omics_contribution()`, `classify_omics_pattern()`, and `get_mnsea_contribution()` for contribution summaries.
- **Topology-aware extraction**: `extract_mnsea_subnetwork()` for pathway-specific node/edge tables that can be passed to downstream visualization packages.
- **Bayesian term selection**: `bayes_enrich()` and `bayes_summary()` for posterior-based term prioritization.

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

The package returns the following S4 result classes:

- `enrichResult` for ORA-like workflows
- `gseaResult` for ranked enrichment workflows
- `nseaResult` for single-network propagation plus enrichment
- `mnseaResult` for multi-layer propagation, collapsed scores, and cached explanation tables

These classes are used across the `clusterProfiler` family:

- `enrichit` handles core computation, algorithm implementation, and contribution data preparation
- `clusterProfiler` provides high-level biological interpretation workflows and general enrichment analysis interfaces
- `enrichplot` handles visualization
- `gson` provides a structured gene set resource layer for managing and exchanging gene set collections across the family
- knowledge-base-oriented downstream packages such as `DOSE`, `ReactomePA`, `meshes`, and `MicrobiomeProfiler` provide domain-specific annotation and interpretation layers
