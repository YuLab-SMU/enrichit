# enrichit 0.2.5

- ORA results now carry an **`oddsRatio`** column (Fisher's exact 2x2 odds ratio, placed next to `FoldEnrichment`). `FoldEnrichment` is a ratio of proportions, whereas the odds ratio is the effect size that the hypergeometric/Fisher test is built on; note it is the plain cross-product odds ratio and differs slightly from `fisher.test()$estimate`, which reports the conditional MLE (2026-09-22, Tue)
- report non-finite gene statistics instead of crashing with `missing value where TRUE/FALSE needed`: `gsea_gson()` checked `is.unsorted()` before validating the input, and a bare `if (NA)` aborted there, so an `NA`/`NaN`/`Inf` in the gene statistics never reached the finiteness check that explains the problem; `gseaScores()` now validates too and returns `NA` rather than aborting when the running score is not finite (2026-09-22, Tue)
- add generic constructors `as_enrichResult()` / `as_gseaResult()`: build 'enrichit' result objects from result tables of external enrichment tools (enrichr, g:Profiler, WebGestalt, fgsea, ...), with canonical column aliases, ID de-duplication, p-value validation, derived statistics (`RichFactor` / `FoldEnrichment` / `zScore`) when the query and background are known, and missing GSEA detail columns recomputed via `gsea_leading_edge_details()` (2026-09-21, Sun)

# enrichit 0.2.4

- stop emitting spurious `no package '<...>' was found` warnings when no input gene can be mapped (2026-09-17, Thu)
  - `check_gene_id()` routed its informational notices through `yulab.utils::yulab_msg()`, which builds the package *citation banner* and expects a package name; it therefore called `packageDescription()` on the notice text itself, so every run with unmappable input emitted three warnings such as `no package '--> No gene can be mapped....' was found`
  - the notices now go through `message()`, and `yulab.utils::yulab_msg()` is no longer imported
  - the sample of expected gene IDs shown in the notice no longer contains `NA` when a gene set has fewer than 100 genes

# enrichit 0.2.3

- calibrate NSEA significance testing by switching the default to a whole-pipeline permutation null instead of GSEA's label-permutation test (2026-08-23, Sun)
  - the legacy `nsea()`/`nsea_gson()` pipeline ran GSEA's label-permutation test on the network-diffused scores; because network diffusion induces strong autocorrelation between neighbouring genes' scores, that test violates the exchangeability assumption of GSEA's permutation null and produces anti-conservative p-values (empirically \~2.5x the nominal false-positive rate)
  - `significance = "whole_pipeline"` (new default) builds the null distribution of the enrichment score by re-running the *entire* pipeline under the null: permute gene labels -> re-diffuse over the network with RWR -> recompute the enrichment score; p-values and NES are then derived from this null, automatically accounting for the smoothing induced by diffusion
  - `significance = "internal"` retains the legacy behaviour (GSEA's internal label-permutation test on the diffused scores), kept for comparison/debugging
- add `significance`, `nPerm`, `seed` and `exponent` arguments to `nsea()` and `nsea_gson()`; expose `pvalueCutoff` and `pAdjustMethod` in `nsea_gson()`
  - `nPerm` (default 1000) sets the number of whole-pipeline permutations; the smallest estimable p-value is `1/(nPerm + 1)`
  - `seed` (default NULL) makes the whole-pipeline null reproducible; set a numeric seed for deterministic results
  - `exponent` (default 1) controls the weight of each step of the running enrichment score
  - `nsea_gson()` now applies `pvalueCutoff` to both the raw and the adjusted p-value (matching `gsea_gson()`), and accepts `pAdjustMethod`
- extend `prepare_network()` to accept a `mechgraph` object as well as an edge-list `data.frame`/`matrix` or a sparse matrix
  - the `mechgraph` `edges` table must have `from`/`to` columns; a numeric `score` or `weight` column, when present, supplies the edge weight, otherwise unit weights are used
  - the class is duck-typed by name, so `enrichit` does not hard-depend on the `mechgraph` package
- add a fast C\++ kernel `gsea_es_all_cpp()` that computes the classic weighted enrichment score for many gene sets at once from the hit positions only, making it cheap enough to evaluate on every whole-pipeline permutation
- harden the NSEA internals
  - robust node-name extraction guards against Matrix accessor quirks on serialized sparse matrices
  - share single implementations of the RWR diffusion (`.nsea_diffuse()`) and TF-IDF specificity weighting (`.nsea_specific_weights()`); in "signed" mode positive and negative seeds are propagated in a single linear solve
- add regression tests for whole-pipeline reproducibility and for p-value calibration on null data

# enrichit 0.2.2

- align `gsea_gson()` p-value filtering with the historical clusterProfiler/DOSE behavior: `pvalueCutoff` now requires both the raw p-value and the adjusted p-value (`p.adjust`) to pass the cutoff (previously only the raw p-value was filtered), restoring significant-pathway counts comparable to clusterProfiler <= 4.18.x (2026-08-14, Fri)
- clarify and harden the `seed` interface of GSEA for reproducibility (2026-08-14, Thu)
  - `gsea()` now treats `seed = TRUE` as a fixed default seed (consistent with the C\++ default) instead of silently coercing it to the integer `1`
  - `gsea_gson()` now exposes an explicit `seed` argument (previously only reachable through `...`) and forwards it to `gsea()`
  - document that `seed = FALSE` (default) draws a fresh seed from R's RNG on each run, so results may vary between runs, while a numeric seed (or `set.seed()` before the call) makes the result reproducible; the C\++ engine seeds its own RNG with this value
  - add regression tests asserting identical results across runs with a fixed seed

# enrichit 0.2.1

- fix `gsea()` to intersect gene sets with `names(geneList)` before applying `minGSSize`/`maxGSSize`, so the size filter constrains the actual overlap rather than the raw gene set size; also guard `gsea_gson()` against `NA` pvalue rows leaking into the result table (2026-08-04, Tue, clusterProfiler#824)

# enrichit 0.2.0

- extend multi-omics integration with pathway-level and topology-level workflows (2026-06-24, Wed)
  - add `aggregate_enrichment()` for multi-omics Late Fusion at the pathway level by aggregating multiple `enrichResult`/`gseaResult` objects
  - support correlation-aware p-value aggregation with `method = "brown"` in `aggregate_omics()`
  - support `method = "weighted_mean"` in `aggregate_omics()` for signed statistics with layer-specific weights
  - add `mnsea()` and `mnsea_gson()` for multi-layer network-based enrichment
  - add `prepare_multilayer_network()`, `propagate_multilayer()`, and `collapse_multilayer_scores()` for the multi-layer propagation pipeline
  - add `mnseaResult` to store multi-layer diffusion results, collapsed scores, layer weights, and cached explanation tables
  - precompute pathway-level and feature-level explanation caches inside `mnseaResult`
  - add `get_mnsea_contribution()` and `extract_mnsea_subnetwork()` for explanation-ready data extraction
  - provide the generic engine layer for downstream high-level wrappers such as `clusterProfiler::mnseGO()`, `mnseKEGG()`, `mnseMKEGG()`, and `mnseWP()`
- add `aggregate_omics()`, `harmonize_ids()` and `select_features_for_ora()` to support Multi-omics Early Integration
  - add `conflict_policy` parameter ("keep_all", "strict", "penalty") to handle directional conflicts in signed statistics (2026-06-23, Tue)
  - support early fusion before `ora()`, `gsea()`, and `nsea()` through a decoupled aggregation layer
  - these workflow helpers, together with `aggregate_enrichment()`, are designed to be reused by downstream packages for high-level multi-omics analysis
- add `get_omics_contribution()` and `classify_omics_pattern()` for Multi-omics contribution tracing
- implement Weighted Enrichment Analysis (2026-06-23, Tue)
  - add `weight` parameter to `ora()`, `ora_gson()`, `gsea()`, and `gsea_gson()`
  - support Weighted ORA using Wallenius' noncentral hypergeometric distribution via the `BiasedUrn` package
  - support Weighted GSEA by fusing external weights with ranked statistics
- implement Network-based Set Enrichment Analysis (NSEA) (2026-06-23, Tue)
  - add `nsea()` and `nsea_gson()` for network-ranked GSEA based on Random Walk with Restart (RWR)
  - add `mode = "signed"` support in `nsea()` and `nsea_gson()` for bidirectional network propagation using signed statistics
  - add `prepare_network()` for parsing and normalizing edge lists or sparse matrices
  - implement extremely fast RWR using `RcppEigen` sparse matrix multiplication
  - introduce zero-dependency integration strategy for network propagation followed by multilevel GSEA
  - provide the generic engine layer for downstream high-level wrappers such as `clusterProfiler::nseGO()`, `nseKEGG()`, `nseMKEGG()`, and `nseWP()`
- align multilevel GSEA rank scaling with `fgsea::prepareStats()` to reduce result drift relative to the long-used fgsea backend (2026-06-22, Mon)
  - replace the fixed `* 1e6` scaling in `prepare_gsea_inputs()` with fgsea-style total-weight normalization and integer rounding
  - add a regression test that compares `gsea(method = "multilevel")` against `fgsea::fgseaMultilevel()` on the same ranked input
- exclude zero-overlap gene sets from `ora_gson()` before multiple-testing correction (2026-06-22, Mon)
  - keep `Count = 0` rows out of `p.adjust`/`qvalue` so ORA results match historical `DOSE`/`clusterProfiler` behavior
  - resolves inflated adjustment in downstream `clusterProfiler::enricher()`, `enrichKEGG()`, and `compareCluster()` workflows (e.g. compound KEGG analyses, #821 & #819 of 'clusterProfiler')

# enrichit 0.1.5

- add Bayesian term selection for ORA results (2026-06-16, Tue)
  - add `bayes_enrich()` for estimating posterior probabilities of active explanatory terms from `enrichResult` objects
  - add `bayes_summary()` to return Bayesian enrichment results ordered by posterior probability
  - support candidate term spaces from top ranked result rows, significant `as.data.frame()` rows, all `x@result` rows, or explicit term IDs
  - keep the output as an `enrichResult` object with posterior, posterior odds, Bayesian rank, active flag, and covered-gene columns

# enrichit 0.1.4

- fix GSON helper functions (2026-04-06, Mon)
  - correct `TERM2NAME()` for `GSON` objects to map `gsid` to term names properly
  - correct `TERMID2EXTID()` for `GSON` objects to return gene vectors by requested term order
  - return `character(0)` for missing terms to keep downstream behavior stable
  - resolves malformed `Description.*` columns in downstream `clusterProfiler::groupGO()`

# enrichit 0.1.3

- fix `gsea_gson()` (2026-03-10, Tue)
  - force sorting of `geneList` to ensure result object is consistent with input
  - resolves issue with `enrichplot::gseaplot` showing incorrect metric/color alignment

# enrichit 0.1.2

- improve robustness of `calculate_qvalue()` (2026-02-02, Mon)
  - handle missing `qvalue` package gracefully
  - retry with different parameters if `qvalue()` fails
  - handle invalid p-values (NA, infinite, out of range)
- validate p-values in `ora_gson()` (2026-02-02, Mon)
  - ensure p-values are within [0, 1] range
  - warn and report invalid p-values
- optimize ORA p-value calculation in C\++ (2026-02-02, Mon)
  - use `phyper` instead of summing `dhyper` for better performance and precision

# enrichit 0.1.1

- update `setReadable()` to support converting gene ID to other types (not limited to SYMBOL) (2026-01-21, Wed)
- add `organism` slot in `compareClusterResult` (2026-01-20, Tue)

# enrichit 0.1.0

- fix bugs in `gsea_gson()` and `ora_gson()` (2026-01-11, Sun)
  - handle missing columns (e.g., qvalues) gracefully by filling with NA
  - handle `NA` or duplicate gene set IDs in result rownames to prevent errors
- improve robustness of `calculate_qvalue()` (2026-01-11, Sun)
  - return NA instead of NULL when qvalue calculation fails
- update `ora_gson()` output columns (2026-01-11, Sun)

# enrichit 0.0.9

- add leading edge analysis for GSEA (2026-01-10, Sat)

# enrichit 0.0.8

- fixed bugs of multilevel GSEA in p value calculation (2025-12-10, Wed)
  - by learning the source code of 'fgsea'

# enrichit 0.0.7

- add `gseaScores` function (2025-12-07, Sun)
  - to calculate GSEA scores for a single gene set

# enrichit 0.0.6

- add vignette (2025-12-07, Sun)

# enrichit 0.0.5

- implement multi-level GSEA algorithm (2025-12-06, Sat)

# enrichit 0.0.4

- implement a simplified adaptive early-stopping GSEA algorithm (2025-12-05, Fri)
  - For each gene set:
  - 1. Run initial batch (e.g., 1000 permutations)
  - 2. If p-value > threshold (e.g., 0.05), stop
  - 3. If significant, increase permutations geometrically (2x, 4x, 8x...)
  - 4. Continue until p-value stabilizes or max permutations reached

# enrichit 0.0.3

- implement `ora_gson` and `gsea_gson` (2025-12-05, Fri)
  - as replacement for `enricher_internal` and `GSEA_internal`
- mv helper functions and class definitions from `DOSE` to `enrichit` (2025-12-05, Fri)
  - to extend this package as the base package for the `clusterProfiler` family

# enrichplot 0.0.2

- `gsea` function (2025-12-04, Thu)
  - Gene Set Enrichment Analysis (GSEA) using C\++ via Rcpp.

# enrichit 0.0.1

- `ora` function (2025-12-03, Wed)
  - Fast Over-Representation Analysis (ORA) using C\++ via Rcpp.

