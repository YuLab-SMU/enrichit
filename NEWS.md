# enrichit 0.2.1

+ fix `gsea()` to intersect gene sets with `names(geneList)` before applying `minGSSize`/`maxGSSize`, so the size filter constrains the actual overlap rather than the raw gene set size; also guard `gsea_gson()` against `NA` pvalue rows leaking into the result table (2026-08-04, Tue, clusterProfiler#824)

# enrichit 0.2.0

+ extend multi-omics integration with pathway-level and topology-level workflows (2026-06-24, Wed)
    - add `aggregate_enrichment()` for multi-omics Late Fusion at the pathway level by aggregating multiple `enrichResult`/`gseaResult` objects
    - support correlation-aware p-value aggregation with `method = "brown"` in `aggregate_omics()`
    - support `method = "weighted_mean"` in `aggregate_omics()` for signed statistics with layer-specific weights
    - add `mnsea()` and `mnsea_gson()` for multi-layer network-based enrichment
    - add `prepare_multilayer_network()`, `propagate_multilayer()`, and `collapse_multilayer_scores()` for the multi-layer propagation pipeline
    - add `mnseaResult` to store multi-layer diffusion results, collapsed scores, layer weights, and cached explanation tables
    - precompute pathway-level and feature-level explanation caches inside `mnseaResult`
    - add `get_mnsea_contribution()` and `extract_mnsea_subnetwork()` for explanation-ready data extraction
    - provide the generic engine layer for downstream high-level wrappers such as `clusterProfiler::mnseGO()`, `mnseKEGG()`, `mnseMKEGG()`, and `mnseWP()`
+ add `aggregate_omics()`, `harmonize_ids()` and `select_features_for_ora()` to support Multi-omics Early Integration
    - add `conflict_policy` parameter ("keep_all", "strict", "penalty") to handle directional conflicts in signed statistics (2026-06-23, Tue)
    - support early fusion before `ora()`, `gsea()`, and `nsea()` through a decoupled aggregation layer
    - these workflow helpers, together with `aggregate_enrichment()`, are designed to be reused by downstream packages for high-level multi-omics analysis
+ add `get_omics_contribution()` and `classify_omics_pattern()` for Multi-omics contribution tracing
+ implement Weighted Enrichment Analysis (2026-06-23, Tue)
    - add `weight` parameter to `ora()`, `ora_gson()`, `gsea()`, and `gsea_gson()`
    - support Weighted ORA using Wallenius' noncentral hypergeometric distribution via the `BiasedUrn` package
    - support Weighted GSEA by fusing external weights with ranked statistics
+ implement Network-based Set Enrichment Analysis (NSEA) (2026-06-23, Tue)
    - add `nsea()` and `nsea_gson()` for network-ranked GSEA based on Random Walk with Restart (RWR)
    - add `mode = "signed"` support in `nsea()` and `nsea_gson()` for bidirectional network propagation using signed statistics
    - add `prepare_network()` for parsing and normalizing edge lists or sparse matrices
    - implement extremely fast RWR using `RcppEigen` sparse matrix multiplication
    - introduce zero-dependency integration strategy for network propagation followed by multilevel GSEA
    - provide the generic engine layer for downstream high-level wrappers such as `clusterProfiler::nseGO()`, `nseKEGG()`, `nseMKEGG()`, and `nseWP()`
+ align multilevel GSEA rank scaling with `fgsea::prepareStats()` to reduce result drift relative to the long-used fgsea backend (2026-06-22, Mon)
    - replace the fixed `* 1e6` scaling in `prepare_gsea_inputs()` with fgsea-style total-weight normalization and integer rounding
    - add a regression test that compares `gsea(method = "multilevel")` against `fgsea::fgseaMultilevel()` on the same ranked input
+ exclude zero-overlap gene sets from `ora_gson()` before multiple-testing correction (2026-06-22, Mon)
    - keep `Count = 0` rows out of `p.adjust`/`qvalue` so ORA results match historical `DOSE`/`clusterProfiler` behavior
    - resolves inflated adjustment in downstream `clusterProfiler::enricher()`, `enrichKEGG()`, and `compareCluster()` workflows (e.g. compound KEGG analyses, #821 & #819 of 'clusterProfiler')

# enrichit 0.1.5

+ add Bayesian term selection for ORA results (2026-06-16, Tue)
    - add `bayes_enrich()` for estimating posterior probabilities of active explanatory terms from `enrichResult` objects
    - add `bayes_summary()` to return Bayesian enrichment results ordered by posterior probability
    - support candidate term spaces from top ranked result rows, significant `as.data.frame()` rows, all `x@result` rows, or explicit term IDs
    - keep the output as an `enrichResult` object with posterior, posterior odds, Bayesian rank, active flag, and covered-gene columns

# enrichit 0.1.4

+ fix GSON helper functions (2026-04-06, Mon)
    - correct `TERM2NAME()` for `GSON` objects to map `gsid` to term names properly
    - correct `TERMID2EXTID()` for `GSON` objects to return gene vectors by requested term order
    - return `character(0)` for missing terms to keep downstream behavior stable
    - resolves malformed `Description.*` columns in downstream `clusterProfiler::groupGO()`

# enrichit 0.1.3

+ fix `gsea_gson()` (2026-03-10, Tue)
    - force sorting of `geneList` to ensure result object is consistent with input
    - resolves issue with `enrichplot::gseaplot` showing incorrect metric/color alignment

# enrichit 0.1.2

+ improve robustness of `calculate_qvalue()` (2026-02-02, Mon)
    - handle missing `qvalue` package gracefully
    - retry with different parameters if `qvalue()` fails
    - handle invalid p-values (NA, infinite, out of range)
+ validate p-values in `ora_gson()` (2026-02-02, Mon)
    - ensure p-values are within [0, 1] range
    - warn and report invalid p-values
+ optimize ORA p-value calculation in C++ (2026-02-02, Mon)
    - use `phyper` instead of summing `dhyper` for better performance and precision

# enrichit 0.1.1

+ update `setReadable()` to support converting gene ID to other types (not limited to SYMBOL) (2026-01-21, Wed)
+ add `organism` slot in `compareClusterResult` (2026-01-20, Tue)

# enrichit 0.1.0

+ fix bugs in `gsea_gson()` and `ora_gson()` (2026-01-11, Sun)
    - handle missing columns (e.g., qvalues) gracefully by filling with NA
    - handle `NA` or duplicate gene set IDs in result rownames to prevent errors
+ improve robustness of `calculate_qvalue()` (2026-01-11, Sun)
    - return NA instead of NULL when qvalue calculation fails
+ update `ora_gson()` output columns (2026-01-11, Sun)

# enrichit 0.0.9

+ add leading edge analysis for GSEA (2026-01-10, Sat)

# enrichit 0.0.8

+ fixed bugs of multilevel GSEA in p value calculation (2025-12-10, Wed)
    - by learning the source code of 'fgsea'

# enrichit 0.0.7

+ add `gseaScores` function (2025-12-07, Sun)
    - to calculate GSEA scores for a single gene set

# enrichit 0.0.6

+ add vignette (2025-12-07, Sun)

# enrichit 0.0.5

+ implement multi-level GSEA algorithm (2025-12-06, Sat)
    
# enrichit 0.0.4

+ implement a simplified adaptive early-stopping GSEA algorithm (2025-12-05, Fri)
    - For each gene set:
    - 1. Run initial batch (e.g., 1000 permutations)
    - 2. If p-value > threshold (e.g., 0.05), stop
    - 3. If significant, increase permutations geometrically (2x, 4x, 8x...)
    - 4. Continue until p-value stabilizes or max permutations reached

# enrichit 0.0.3

+ implement `ora_gson` and `gsea_gson` (2025-12-05, Fri)
    - as replacement for `enricher_internal` and `GSEA_internal`
+ mv helper functions and class definitions from `DOSE` to `enrichit` (2025-12-05, Fri)
    - to extend this package as the base package for the `clusterProfiler` family

# enrichplot 0.0.2

+ `gsea` function (2025-12-04, Thu)
    - Gene Set Enrichment Analysis (GSEA) using C++ via Rcpp.

# enrichit 0.0.1

+ `ora` function (2025-12-03, Wed)
    - Fast Over-Representation Analysis (ORA) using C++ via Rcpp.
