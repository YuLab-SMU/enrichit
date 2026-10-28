# enrichit 0.2.7

- GSEA no longer turns a degenerate input into a silently useless table: `eps` outside `[0, 1]` is
  rejected (it returned `p = eps > 1`), an all-zero ranking is an error instead of a frame of
  `ES = NaN`, a gene set covering the whole ranked list is excluded as fgsea does, and gene sets
  dropped for having no estimable p-value are counted instead of vanishing (2026-10-27)
- the multilevel "too few same-side permutations" guard is a rate, not a count, so `nPermSimple`
  below 10 no longer turns every pathway `NA`; default results are unchanged (2026-10-27)
- `qvalueCutoff` is applied **row by row**: one unestimable Storey q-value used to disable the
  cutoff for the whole table, so `qvalueCutoff = 0.2` returned terms with q = 0.9
  (clusterProfiler#710, #584, #614, #332). Related: a failed estimate is now a `message()` rather
  than a `warning()`, and `NA` `pvalue`/`p.adjust` no longer yield all-`NA` rows (2026-10-27)
- NSEA's permutation null no longer counts skipped permutations (`nPerm` reports what ran).
  Not reachable from ordinary input, so this is hardening (2026-10-27)
- ORA results keep `DESize`, `UniverseSize` and `SetSize` as columns, so `GeneRatio == Count/DESize`
  and friends are checkable from the table (clusterProfiler#703, #523, #618, #723) (2026-10-27)
- `bayes_enrich()`: default `false_negative` is 0.9, not 0.1, which drove every posterior to 0 on
  real KEGG data. Also rejects `burnin + thin > n_iter` and warns when no candidate reaches
  `posterior_cutoff` (2026-10-26)
- `dsea()` and `dsea_scores()`: design-aware GSEA and per-sample ssGSEA scores, with the contrast
  never guessed; `as_ranked()` and `as_genes()` adapt `DESeq2` objects for it (2026-10-18)
- `filter_gene_sets()` is vectorized and `ora_gson()` skips the background intersection when the
  background was not restricted: 14,209 GO BP sets 14.0 s -> 0.9 s, single `enrichGO()` 19 s ->
  4.2 s, 13-cluster `compareCluster()` 381 s -> 197 s. `ora()` returns `geneID` in input order
  instead of hash-set order (2026-10-18)
- remove the golden-reference script and data for the retired `camera`/`roast`/`fry` engines, and
  the design note for the never-implemented `dsea_nb` set test (2026-10-18)

# enrichit 0.2.6

- fix permutation-based GSEA p-values in the fixed and adaptive engines
  (`method = "sample"` / `"permute"`) for the default two-sided
  (`scoreType = "std"`) calculation: the numerator counted only same-sign
  permutations while the denominator divided by all of them, so p-values were too
  small under a symmetric null. The denominator now counts only permutations on
  the same side of zero as the observed ES, following the main fgsea convention
  (2026-10-06, closes #3)
- update the `RcppEigen` sparse-matrix interface in `rwr_eigen_cpp()` from
  `MappedSparseMatrix` to `Map<SparseMatrix>` for Eigen 5 compatibility
  (2026-10-06, PR #2)
- add a `reportNA` argument to `gsea_gson()` so pathways with `NA` multilevel
  p-values can be retained (2026-10-04, closes DOSE#88)
- fix GSEA leading-edge edge cases with one `rank = 0L` sentinel for zero-overlap,
  `N_R == 0` and non-finite running-score paths; align `as_gseaResult()` with the
  same contract and document the `gsea_leading_edge_details()` return value
  (2026-10-04–2026-10-05, closes DOSE#46)
- keep `qvalue` as `NA` when estimation fails instead of falling back to
  `p.adjust`, and document the `pvalue` / `p.adjust` / `qvalue` semantics in
  `as_enrichResult()` and `as_gseaResult()` (2026-10-05, addresses DOSE#24)

# enrichit 0.2.5

- ORA results carry an **`oddsRatio`** column (Fisher's exact 2x2 odds ratio, next
  to `FoldEnrichment`). `FoldEnrichment` is a ratio of proportions; the odds ratio
  is the effect size the hypergeometric test is built on. It is the plain
  cross-product odds ratio and differs slightly from `fisher.test()$estimate`,
  which reports the conditional MLE (2026-09-22)
- report non-finite gene statistics instead of aborting with
  `missing value where TRUE/FALSE needed`. `gsea_gson()` checked `is.unsorted()`
  before validating its input, and a bare `if (NA)` stopped there, so an
  `NA` / `NaN` / `Inf` never reached the finiteness check that explains the
  problem. `gseaScores()` now validates too and returns `NA` rather than aborting
  when the running score is not finite (2026-09-22)
- add generic constructors `as_enrichResult()` / `as_gseaResult()` to build
  'enrichit' result objects from external result tables (enrichr, g:Profiler,
  WebGestalt, fgsea, ...), with canonical column aliases, ID de-duplication,
  p-value validation, derived statistics (`RichFactor` / `FoldEnrichment` /
  `zScore`) when the query and background are known, and missing GSEA detail
  columns recomputed via `gsea_leading_edge_details()` (2026-09-21)

# enrichit 0.2.4

- stop emitting spurious `no package '<...>' was found` warnings when no input
  gene can be mapped (2026-09-17). `check_gene_id()` routed its notices through
  `yulab.utils::yulab_msg()`, which builds the package *citation banner* and
  expects a package name, so it called `packageDescription()` on the notice text
  itself and emitted three such warnings per run. The notices now go through
  `message()`, `yulab.utils::yulab_msg()` is no longer imported, and the sample of
  expected gene IDs no longer contains `NA` for sets with fewer than 100 genes

# enrichit 0.2.3

- calibrate NSEA significance testing by switching the default to a
  whole-pipeline permutation null instead of GSEA's label-permutation test
  (2026-08-23). The legacy pipeline ran the label-permutation test on the
  network-diffused scores; diffusion induces strong autocorrelation between
  neighbouring genes, which violates the exchangeability assumption of that null
  and produces anti-conservative p-values (~2.5x the nominal false-positive rate).
  `significance = "whole_pipeline"` rebuilds the null of the enrichment score by
  permuting gene labels, re-diffusing with RWR and recomputing the score;
  `significance = "internal"` keeps the earlier formulation
- add `significance`, `nPerm`, `seed` and `exponent` to `nsea()` and
  `nsea_gson()`, and expose `pvalueCutoff` / `pAdjustMethod` in `nsea_gson()`.
  `nPerm` (default 1000) sets the number of whole-pipeline permutations with a
  smallest estimable p-value of `1/(nPerm + 1)`; `seed` makes the null
  reproducible; `nsea_gson()` now applies `pvalueCutoff` to both the raw and the
  adjusted p-value, matching `gsea_gson()`
- extend `prepare_network()` to accept a `mechgraph` object as well as an
  edge-list `data.frame` / `matrix` or a sparse matrix. The `edges` table must
  have `from` / `to` columns; a numeric `score` or `weight` column supplies the
  edge weight, otherwise unit weights are used. The class is duck-typed by name,
  so `enrichit` does not hard-depend on the `mechgraph` package
- add a fast C\++ kernel `gsea_es_all_cpp()` that computes the classic weighted
  enrichment score for many gene sets at once from the hit positions only, cheap
  enough to evaluate on every whole-pipeline permutation
- share single implementations of the RWR diffusion (`.nsea_diffuse()`) and
  TF-IDF specificity weighting (`.nsea_specific_weights()`); in `signed` mode
  positive and negative seeds propagate in a single linear solve. Robust
  node-name extraction guards against `Matrix` accessor quirks on serialized
  sparse matrices
- add regression tests for whole-pipeline reproducibility and p-value calibration
  on null data

# enrichit 0.2.2

- align `gsea_gson()` p-value filtering with historical clusterProfiler/DOSE
  behaviour: `pvalueCutoff` now requires both the raw and the adjusted p-value to
  pass (previously only the raw p-value was filtered), restoring significant
  pathway counts comparable to clusterProfiler <= 4.18.x (2026-08-14)
- clarify and harden the `seed` interface of GSEA (2026-08-14). `gsea()` treats
  `seed = TRUE` as a fixed default seed (consistent with the C\++ default) instead
  of silently coercing it to the integer `1`; `gsea_gson()` exposes an explicit
  `seed` argument; `seed = FALSE` (default) draws a fresh seed from R's RNG on
  each run, while a numeric seed (or `set.seed()` before the call) makes the
  result reproducible. Regression tests assert identical results across runs with
  a fixed seed

# enrichit 0.2.1

- fix `gsea()` to intersect gene sets with `names(geneList)` before applying
  `minGSSize` / `maxGSSize`, so the size filter constrains the actual overlap
  rather than the raw gene-set size; also guard `gsea_gson()` against `NA` pvalue
  rows leaking into the result table (2026-08-04, clusterProfiler#824)

# enrichit 0.2.0

- multi-omics integration at pathway and topology level (2026-06-24):
  `aggregate_enrichment()` for Late Fusion over several `enrichResult` /
  `gseaResult` objects; correlation-aware p-value aggregation
  (`method = "brown"`) and `method = "weighted_mean"` for signed statistics with
  layer-specific weights in `aggregate_omics()`; `mnsea()` / `mnsea_gson()` for
  multi-layer network-based enrichment, with `prepare_multilayer_network()`,
  `propagate_multilayer()` and `collapse_multilayer_scores()`; `mnseaResult` to
  store diffusion results, collapsed scores, layer weights and cached explanation
  tables; `get_mnsea_contribution()` and `extract_mnsea_subnetwork()` for
  explanation-ready extraction. Provides the engine layer for
  `clusterProfiler::mnseGO()`, `mnseKEGG()`, `mnseMKEGG()` and `mnseWP()`
- Early Integration helpers (2026-06-23): `aggregate_omics()`,
  `harmonize_ids()` and `select_features_for_ora()`, with a `conflict_policy`
  argument (`"keep_all"` / `"strict"` / `"penalty"`) for directional conflicts in
  signed statistics; early fusion before `ora()`, `gsea()` and `nsea()` through a
  decoupled aggregation layer
- `get_omics_contribution()` and `classify_omics_pattern()` for multi-omics
  contribution tracing
- Weighted Enrichment Analysis (2026-06-23): a `weight` argument in `ora()`,
  `ora_gson()`, `gsea()` and `gsea_gson()`; Weighted ORA uses Wallenius'
  noncentral hypergeometric distribution via `BiasedUrn`, Weighted GSEA fuses
  external weights with the ranked statistics
- Network-based Set Enrichment Analysis (2026-06-23): `nsea()` and
  `nsea_gson()` based on Random Walk with Restart, with `mode = "signed"` for
  bidirectional propagation; `prepare_network()` for parsing and normalizing edge
  lists or sparse matrices; fast RWR via `RcppEigen` sparse multiplication; a
  zero-dependency integration strategy for diffusion followed by multilevel GSEA.
  Provides the engine layer for `clusterProfiler::nseGO()`, `nseKEGG()`,
  `nseMKEGG()` and `nseWP()`
- align multilevel GSEA rank scaling with `fgsea::prepareStats()` to reduce drift
  against the long-used fgsea backend (2026-06-22): replace the fixed `* 1e6`
  scaling in `prepare_gsea_inputs()` with fgsea-style total-weight normalization
  and integer rounding, with a regression test against `fgsea::fgseaMultilevel()`
- exclude zero-overlap gene sets from `ora_gson()` before multiple-testing
  correction (2026-06-22): `Count = 0` rows stay out of `p.adjust` / `qvalue` so
  ORA results match historical `DOSE` / `clusterProfiler` behaviour, resolving
  inflated adjustment in `enricher()`, `enrichKEGG()` and `compareCluster()`
  (clusterProfiler#821, #819)

# enrichit 0.1.5

- add Bayesian term selection for ORA results (2026-06-16): `bayes_enrich()`
  estimates posterior probabilities of active explanatory terms from
  `enrichResult` objects, `bayes_summary()` returns them ordered by posterior.
  Candidate term spaces may come from top ranked rows, significant
  `as.data.frame()` rows, all `x@result` rows, or explicit term IDs. The output
  stays an `enrichResult` with posterior, posterior odds, Bayesian rank, active
  flag and covered-gene columns

# enrichit 0.1.4

- fix GSON helpers (2026-04-06): `TERM2NAME()` maps `gsid` to term names,
  `TERMID2EXTID()` returns gene vectors in the requested term order, and missing
  terms return `character(0)`. Resolves malformed `Description.*` columns in
  `clusterProfiler::groupGO()`

# enrichit 0.1.3

- fix `gsea_gson()` (2026-03-10): force sorting of `geneList` so the result object
  is consistent with the input, which resolves incorrect metric/colour alignment
  in `enrichplot::gseaplot()`

# enrichit 0.1.2

- improve robustness of `calculate_qvalue()` (2026-02-02): handle a missing
  `qvalue` package, retry with different parameters, and handle invalid p-values
  (`NA`, infinite, out of range)
- validate p-values in `ora_gson()` (2026-02-02), warning on values outside
  [0, 1]
- use `phyper` instead of summing `dhyper` for better performance and precision
  in the C\++ ORA p-value calculation (2026-02-02)

# enrichit 0.1.1

- `setReadable()` converts gene IDs to any type, not only SYMBOL (2026-01-21)
- add an `organism` slot to `compareClusterResult` (2026-01-20)

# enrichit 0.1.0

- fix bugs in `gsea_gson()` and `ora_gson()` (2026-01-11): fill missing columns
  (e.g. qvalues) with `NA`, and handle `NA` or duplicate gene-set IDs in result
  rownames
- `calculate_qvalue()` returns `NA` instead of `NULL` when qvalue estimation
  fails (2026-01-11)
- update `ora_gson()` output columns (2026-01-11)

# enrichit 0.0.9

- add leading edge analysis for GSEA (2026-01-10)

# enrichit 0.0.8

- fix p-value calculation in multilevel GSEA (2025-12-10), by learning the source
  code of 'fgsea'

# enrichit 0.0.7

- add `gseaScores()` to calculate GSEA scores for a single gene set (2025-12-07)

# enrichit 0.0.6

- add vignette (2025-12-07)

# enrichit 0.0.5

- implement the multilevel GSEA algorithm (2025-12-06)

# enrichit 0.0.4

- implement a simplified adaptive early-stopping GSEA algorithm (2025-12-05): run
  an initial batch of permutations, stop if p > threshold, otherwise increase the
  permutation count geometrically until the p-value stabilizes or the maximum is
  reached

# enrichit 0.0.3

- implement `ora_gson()` and `gsea_gson()` as replacements for
  `enricher_internal` and `GSEA_internal` (2025-12-05)
- move helper functions and class definitions from `DOSE` to `enrichit` (2025-12-05),
  to extend this package as the base package for the `clusterProfiler` family

# enrichit 0.0.2

- `gsea()`: Gene Set Enrichment Analysis using C\++ via Rcpp (2025-12-04)

# enrichit 0.0.1

- `ora()`: fast Over-Representation Analysis using C\++ via Rcpp (2025-12-03)
