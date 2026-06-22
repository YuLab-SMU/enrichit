# Benchmark script to separate backend drift from Monte Carlo variance
library(fgsea)
library(DOSE)
devtools::load_all(".")

DEFAULT_CFG <- list(
  label = "default",
  sampleSize = 101L,
  nPermSimple = 1000L,
  eps = 1e-50
)

PRECISE_CFG <- list(
  label = "precise",
  sampleSize = 501L,
  nPermSimple = 10000L,
  eps = 1e-50
)

BENCHMARK_SEEDS <- c(1L, 2L)
SUBSET_SEED <- 123L
N_PATHWAYS <- 100L
MIN_SIZE <- 15L
MAX_SIZE <- 500L

summarize_pair <- function(a, b) {
  stopifnot(length(a) == length(b))
  data.frame(
    n = length(a),
    mae = mean(abs(a - b)),
    mae_log10 = mean(abs(log10(a) - log10(b))),
    max_log10 = max(abs(log10(a) - log10(b))),
    spearman = cor(a, b, method = "spearman")
  )
}

compare_results <- function(lhs, rhs, lhs_name, rhs_name) {
  lhs_ids <- lhs$ID
  rhs_ids <- rhs$ID
  common <- intersect(lhs_ids, rhs_ids)

  lhs_sub <- lhs[match(common, lhs_ids), , drop = FALSE]
  rhs_sub <- rhs[match(common, rhs_ids), , drop = FALSE]

  data.frame(
    ID = common,
    lhs_p = lhs_sub$pvalue,
    rhs_p = rhs_sub$pvalue,
    lhs_es = lhs_sub$enrichmentScore,
    rhs_es = rhs_sub$enrichmentScore,
    ABS_DIFF = abs(lhs_sub$pvalue - rhs_sub$pvalue),
    LOG10_DIFF = abs(log10(lhs_sub$pvalue) - log10(rhs_sub$pvalue)),
    RATIO = lhs_sub$pvalue / rhs_sub$pvalue,
    lhs = lhs_name,
    rhs = rhs_name
  )
}

to_enrichit_like <- function(fgsea_res) {
  data.frame(
    ID = as.character(fgsea_res$pathway),
    pvalue = as.numeric(fgsea_res$pval),
    enrichmentScore = as.numeric(fgsea_res$ES),
    NES = as.numeric(fgsea_res$NES),
    stringsAsFactors = FALSE
  )
}

run_fgsea <- function(pathways_subset, stats, seed, cfg) {
  set.seed(seed)
  start_time <- Sys.time()
  res <- fgseaMultilevel(
    pathways_subset,
    stats,
    minSize = MIN_SIZE,
    maxSize = MAX_SIZE,
    eps = cfg$eps,
    sampleSize = cfg$sampleSize,
    nPermSimple = cfg$nPermSimple
  )
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  list(result = to_enrichit_like(as.data.frame(res)), time = elapsed)
}

run_enrichit <- function(pathways_subset, stats, seed, cfg) {
  start_time <- Sys.time()
  res <- gsea(
    stats,
    pathways_subset,
    method = "multilevel",
    minGSSize = MIN_SIZE,
    maxGSSize = MAX_SIZE,
    pvalThreshold = 1.0,
    eps = cfg$eps,
    sampleSize = cfg$sampleSize,
    nPermSimple = cfg$nPermSimple,
    seed = seed
  )
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  list(result = res, time = elapsed)
}

print_comparison <- function(df, title) {
  cat("\n", title, "\n", sep = "")
  cat("Top 10 largest absolute differences:\n")
  print(head(df[order(df$ABS_DIFF, decreasing = TRUE), ], 10))
  cat("\nTop 10 largest log10 differences:\n")
  print(head(df[order(df$LOG10_DIFF, decreasing = TRUE), ], 10))
  cat("\nP-value ratio summary (lhs/rhs):\n")
  print(summary(df$RATIO))
  cat("\nGeometric mean of p-value ratio:", exp(mean(log(df$RATIO))), "\n")
  cat("Median p-value ratio:", median(df$RATIO), "\n")
}

# Load the same data as in benchmark_vs_fgsea.R
data(geneList, package = "DOSE")
stats <- geneList
cat("Stats summary:\n")
print(summary(stats))
print(head(stats))
print(tail(stats))

x <- DOSE:::get_dose_data("HDO")
pathways <- split(x@gsid2gene$gene, x@gsid2gene$gsid)

set.seed(SUBSET_SEED)
subset_idx <- sample(length(pathways), min(N_PATHWAYS, length(pathways)))
pathways_subset <- pathways[subset_idx]

cat("Running on", length(pathways_subset), "pathways\n")
cat("Seeds:", paste(BENCHMARK_SEEDS, collapse = ", "), "\n")

run_suite <- function(cfg) {
  cat("\n=== Configuration:", cfg$label, "===\n")
  cat("sampleSize =", cfg$sampleSize, "; nPermSimple =", cfg$nPermSimple, "; eps =", cfg$eps, "\n")

  fgsea_runs <- lapply(BENCHMARK_SEEDS, function(seed) run_fgsea(pathways_subset, stats, seed, cfg))
  enrichit_runs <- lapply(BENCHMARK_SEEDS, function(seed) run_enrichit(pathways_subset, stats, seed, cfg))

  fgsea_ref <- fgsea_runs[[1]]$result
  enrichit_ref <- enrichit_runs[[1]]$result

  fgsea_vs_enrichit <- compare_results(fgsea_ref, enrichit_ref, "fgsea", "enrichit")
  fgsea_vs_fgsea <- compare_results(fgsea_ref, fgsea_runs[[2]]$result, "fgsea(seed1)", "fgsea(seed2)")
  enrichit_vs_enrichit <- compare_results(enrichit_ref, enrichit_runs[[2]]$result, "enrichit(seed1)", "enrichit(seed2)")

  cat("\nTime comparison:\n")
  cat("  fgsea seed1:   ", round(fgsea_runs[[1]]$time, 4), "s\n")
  cat("  fgsea seed2:   ", round(fgsea_runs[[2]]$time, 4), "s\n")
  cat("  enrichit seed1:", round(enrichit_runs[[1]]$time, 4), "s\n")
  cat("  enrichit seed2:", round(enrichit_runs[[2]]$time, 4), "s\n")
  cat("  Ratio (enrichit/fgsea, seed1):", round(enrichit_runs[[1]]$time / fgsea_runs[[1]]$time, 2), "x\n")

  summary_table <- rbind(
    fgsea_vs_enrichit = summarize_pair(fgsea_vs_enrichit$lhs_p, fgsea_vs_enrichit$rhs_p),
    fgsea_vs_fgsea = summarize_pair(fgsea_vs_fgsea$lhs_p, fgsea_vs_fgsea$rhs_p),
    enrichit_vs_enrichit = summarize_pair(enrichit_vs_enrichit$lhs_p, enrichit_vs_enrichit$rhs_p)
  )
  print(summary_table)

  cat("\nES correlation (fgsea vs enrichit, seed1):", round(cor(fgsea_vs_enrichit$lhs_es, fgsea_vs_enrichit$rhs_es), 6), "\n")
  cat("If `fgsea vs enrichit` is not worse than `fgsea vs fgsea`, most visible drift is Monte Carlo noise, not backend bias.\n")

  print_comparison(fgsea_vs_enrichit, "fgsea vs enrichit (same seed)")
}

run_suite(DEFAULT_CFG)
run_suite(PRECISE_CFG)
