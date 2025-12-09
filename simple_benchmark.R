# Simplified benchmark to understand p-value differences
library(fgsea)
library(enrichit)

# Load the same data as in benchmark_vs_fgsea.R
data(geneList, package = "DOSE")
stats <- geneList

x <- DOSE:::get_dose_data("HDO")
pathways <- split(x@gsid2gene$gene, x@gsid2gene$gsid)

# Run on a subset for faster testing
set.seed(123)
subset_idx <- sample(length(pathways), min(100, length(pathways)))
pathways_subset <- pathways[subset_idx]

cat("Running on", length(pathways_subset), "pathways\n")

# Run fgsea
cat("Running fgseaMultilevel...\n")
fgsea_res <- fgseaMultilevel(pathways_subset, stats, minSize=15, maxSize=500, eps=1e-50, sampleSize=101)

# Run enrichit
cat("Running enrichit gsea (multilevel)...\n")
enrichit_res <- gsea(stats, pathways_subset, method = "multilevel", minPerm = 101,
                     minGSSize=15, maxGSSize=500, pvalThreshold = 1.0, eps = 1e-50)

# Match by pathway name
common <- intersect(fgsea_res$pathway, enrichit_res$ID)
fgsea_sub <- fgsea_res[match(common, fgsea_res$pathway), ]
enrichit_sub <- enrichit_res[match(common, enrichit_res$ID), ]

pval_fgsea <- fgsea_sub$pval
pval_enrichit <- enrichit_sub$pvalue

# Calculate statistics
cor_pval <- cor(pval_fgsea, pval_enrichit, method = "spearman")
mae <- mean(abs(pval_fgsea - pval_enrichit))
rmse <- sqrt(mean((pval_fgsea - pval_enrichit)^2))
mae_log <- mean(abs(log10(pval_fgsea + 1e-300) - log10(pval_enrichit + 1e-300)))

cat("\n=== Results ===\n")
cat("Correlation (Spearman):", round(cor_pval, 6), "\n")
cat("MAE:", signif(mae, 6), "; RMSE:", signif(rmse, 6), "\n")
cat("MAE on log10 scale:", round(mae_log, 6), "\n")

# Check which gives smaller p-values
enrichit_smaller <- sum(pval_enrichit < pval_fgsea)
fgsea_smaller <- sum(pval_fgsea < pval_enrichit)
equal <- sum(abs(pval_fgsea - pval_enrichit) < 1e-10)

cat("\nP-value comparison:\n")
cat("enrichit gives smaller p-values for", enrichit_smaller, "out of", length(common), "pathways\n")
cat("fgsea gives smaller p-values for", fgsea_smaller, "out of", length(common), "pathways\n")
cat("Equal p-values for", equal, "pathways\n")

# Look at largest differences
df <- data.frame(ID = common, FGSEA = pval_fgsea, ENRICHIT = pval_enrichit,
                 ABS_DIFF = abs(pval_fgsea - pval_enrichit),
                 LOG10_DIFF = abs(log10(pval_fgsea + 1e-300) - log10(pval_enrichit + 1e-300)),
                 RATIO = pval_enrichit / pval_fgsea)

cat("\nTop 10 largest absolute differences:\n")
print(head(df[order(df$ABS_DIFF, decreasing = TRUE), ], 10))

cat("\nTop 10 largest log10 differences:\n")
print(head(df[order(df$LOG10_DIFF, decreasing = TRUE), ], 10))

cat("\nP-value ratio summary (enrichit/fgsea):\n")
print(summary(df$RATIO))

# Check if enrichit consistently gives smaller p-values
cat("\nGeometric mean of p-value ratio:", exp(mean(log(df$RATIO))), "\n")
cat("Median p-value ratio:", median(df$RATIO), "\n")

# Plot if possible
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Log-log plot
  df$log10_fgsea <- log10(df$FGSEA + 1e-300)
  df$log10_enrichit <- log10(df$ENRICHIT + 1e-300)

  p <- ggplot(df, aes(x = log10_fgsea, y = log10_enrichit)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    ggtitle("P-value comparison: enrichit vs fgsea (log10 scale)") +
    xlab("log10(fgsea p-value)") +
    ylab("log10(enrichit p-value)") +
    theme_minimal()

  print(p)

  # Ratio distribution
  p2 <- ggplot(df, aes(x = RATIO)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    scale_x_log10() +
    ggtitle("Distribution of p-value ratio (enrichit/fgsea)") +
    xlab("enrichit p-value / fgsea p-value (log10 scale)") +
    theme_minimal()

  print(p2)
}