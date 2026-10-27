# enrichit 0.2.7.003

- NSEA's whole-pipeline permutation null can no longer be poisoned by a skipped permutation. The
  null matrix was initialised to `0` and a degenerate permutation was skipped with `next`, leaving
  a zero in the null that then counted against every pathway while the p-value denominator stayed
  at `nPerm` -- a bias toward non-significance if the path were ever reachable. The matrix is now
  initialised to `NA`, filled sequentially, the denominator counts only permutations that ran
  (`nPerm` reports that count), and an all-degenerate loop is an error instead of a null of
  nothing. `sum(geneList[common])` is invariant under label permutation, so the path is not
  reachable through ordinary input; the fix turns "unreachable by luck" into "impossible by
  construction" (2026-10-27, Tue)
- ORA results keep the denominators the engine already computed: `DESize` (the GeneRatio
  denominator -- query genes inside the universe), `UniverseSize` (the BgRatio denominator --
  universe intersected with annotated genes) and `SetSize` are now columns of the result table.
  They were computed and then dropped, which is why four users in a row had to ask what the
  denominators of the `20/207` strings actually are (clusterProfiler#703, #523, #618, #723); all
  three ratios are now checkable from the table itself, as they long have been in GSEA results
  (2026-10-27, Tue)

- The multilevel engine's "too few same-side permutations" guard is a **rate**, not a count
  (`src/gsea_multilevel.cpp`). It read `modeFraction < 10` against the raw number of same-side
  permutations, which matches fgsea exactly at the default `nPermSimple = 1000` -- but with any
  smaller `nPermSimple` it turned *every* pathway `NA`: `nPermSimple = 10` returned an all-NA table
  where `nPermSimple = 1000` returned `p = 1e-10` on the same input (measured 2026-10-27). The guard
  is now `modeRate < 0.01`, scale-free, so lowering `nPermSimple` costs precision instead of
  destroying the result; at the default the two forms are identical, and a 200-set input produced
  byte-identical results across all nine reported columns (2026-10-27, Tue)
- `gsea()` rejects `eps` outside `[0, 1]` instead of passing it on: `eps` is the smallest p-value the
  multilevel method can estimate, so `eps = 2` put the floor above every possible p-value and strong
  pathways came back with `p = 2`, which then flowed into `p.adjust()` and `-log10()` in plots
  (fgsea clamps silently; erroring is more useful). The C++ entry point clamps as defence in depth,
  since it is reachable directly through Rcpp (2026-10-27, Tue)
- An **all-zero** ranked statistic is an error, not a table of `NaN`. `scale_fgsea_ranks()` computed
  `2^30 / 0 = Inf`, every rank became `NA`, and the multilevel path returned a full result frame with
  `ES = NaN` and no message. The all-zero case is rejected in `prepare_gsea_inputs()` and the scaler
  guards a direct caller; a ranking that is zero except for one gene remains valid (2026-10-27, Tue)
- A gene set covering the **entire ranked list** is excluded with a message instead of returning
  `ES = NaN`: with no "not in set" positions the running score never falls back, so the score is
  undefined. fgsea excludes such sets through its default `maxSize = length(stats) - 1`; enrichit
  returned the set as an all-NaN row and the accompanying warning blamed "unbalanced gene-level
  statistic values" -- the wrong cause. The `NA` warning now describes the real causes and mentions
  `reportNA = TRUE` (2026-10-27, Tue)
- Dropping gene sets whose p-value could not be estimated now **says how many went**: with the default
  `reportNA = FALSE` a silently shorter table read as "these were the only gene sets tested". The drop
  is extracted into `.gsea_drop_na_rows()` so the message is unit-testable without constructing a
  degenerate permutation null; the source-fingerprint test in `test-gsea.R` follows the extraction
  (2026-10-27, Tue)

- `get_enriched()` applies `qvalueCutoff` **row by row** instead of switching it off wholesale
  when any q-value is `NA`. The old guard (`if (!any(is.na(...)))`) meant a single term whose
  Storey q-value could not be estimated leaked every unfiltered row into the result with no
  message -- the user set `qvalueCutoff = 0.2` and got back terms with q = 0.9
  (clusterProfiler#710, #584, #614, #332). Now: a row with a q-value is kept only when it passes;
  a row with an `NA` q-value bypasses this cutoff (its `pvalue`/`p.adjust` filters still apply)
  with a warning naming the count; an all-`NA` column disables the cutoff with a warning pointing
  at `p.adjust`; a result with no `qvalue` column says that instead of blaming estimation. The
  empty-result case is skipped explicitly -- on a zero-row frame both all-`NA` tests are vacuously
  TRUE, so a result whose terms were all filtered out by `pvalueCutoff` would otherwise warn about
  q-value estimation on every `show()`, `dim()`, `head()`, `tail()` and `as.data.frame()`, since
  all of them route through here. The filter is reached from the `enrichResult` methods only;
  `gseaResult` is filtered at construction in `gsea_gson()`, and the `gseaResult` test covers the
  shared helper rather than a live path. The same row-wise rule replaces the wholesale guard in
  clusterProfiler's `enrichDAVID()`, which carried a copy of it (2026-10-27, Tue)

- `bayes_enrich()`: default `false_negative` is 0.9, not 0.1. The old value
  assumed an ORA hit list holds 90% of each active program's genes; a real KEGG
  pathway with 22 of 158 genes violates that for every candidate, which drove
  every posterior to 0 and made the sampler look broken. At the new rate the
  KEGG / Reactome / GO `geneList` examples give 4 / 18 / 7 active programs (was
  dead / 1 / 2); `prior = 0.3` now floods Reactome with 41, so the hit rate, not
  the prior, was the blind knob. Also warns when no candidate reaches
  `posterior_cutoff`, and a test pins `bayes_active = NA` outside the candidate
  window (2026-10-26)
- `filter_gene_sets()` is vectorized: one id mapping over the concatenated sets
  instead of `lapply(gene_sets, intersect, universe)`, which hashed the universe
  once per set. 14,209 human GO BP sets against 18,842 genes: 14.0 s -> 0.9 s,
  byte-identical output, pinned by a test. Used by ORA, GSEA, NSEA, DSEA and
  Bayes enrichment (2026-10-18)
- `ora()` returns `geneID` in the order the DE genes were supplied. It was built
  by iterating a `std::unordered_set`, so the same analysis returned the same
  genes in a different order in different processes; p-values and counts are
  unchanged (2026-10-18)
- `ora_gson()` intersects the sets with the background only when the background
  was actually restricted. The full annotated list is the identity but not free
  (20 s for those 14,209 sets), so the test is
  `length(extID) < length(annotatedGenes)`, not `universe != NULL`. Single
  ENTREZID `enrichGO()`: 19 s -> 4.2 s; 13-cluster SYMBOL `compareCluster()`:
  381 s -> 197 s, the rest being clusterProfiler's per-cluster keytype mapping
  (2026-10-18)
- `bayes_enrich()` and `bayes_summary()` reject `burnin + thin > n_iter` up
  front. `n_iter = 101, burnin = 100, thin = 2` passed the old checks but left no
  post-burn-in sample, so the sampler divided by zero and every posterior became
  `NaN` — surfacing much later as `missing value where TRUE/FALSE needed`. A
  `n_kept == 0` guard stays behind the check so a future change to the keeping
  rule cannot silently reintroduce it (2026-10-26)
- `get_enriched()` no longer returns all-`NA` rows. Subsetting a data frame with a
  logical that carries `NA` keeps the row but blanks every column
  (`d[c(NA, TRUE), ]` gives one all-`NA` row and one real row), so a term whose
  `pvalue` or `p.adjust` was `NA` came back as a corrupted row with `ID = NA`.
  The cutoff filters now index with `which()` (2026-10-26)
- add `dsea()` and `dsea_scores()`: design-aware GSEA (transform, rank by the
  requested contrast, run the package engine, return BH-adjusted p-values and the
  tested contrast) and per-sample ssGSEA scores with `set_sizes` (scored) and
  `raw_set_sizes` (input). Counts, gene-set identifiers and a full-rank,
  sample-aligned design are validated at the entry point (2026-10-18)
- `dsea()` does not guess the contrast: `contrast = NULL` takes the non-intercept
  column of a two-column design and requires an explicit choice for a wider one;
  intercept and cell-means columns (`~ 0 + condition`) are rejected, because
  ranking by either orders genes by mean expression. Warns when a
  `DESeqDataSet`'s stored design differs from `design`, since the transform reads
  the former and the ranking the latter (2026-10-18)
- add `as_ranked()` and `as_genes()` adapters for `DESeqResults` and fitted
  `DESeqDataSet` objects; `DESeq2` is suggested, not imported (2026-10-18)
- remove the golden-reference script and data for the retired `camera` / `roast` /
  `fry` engines, and the design note for the never-implemented `dsea_nb` set test
  (2026-10-18)

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
