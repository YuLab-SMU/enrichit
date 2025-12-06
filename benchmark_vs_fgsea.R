
if (!requireNamespace("fgsea", quietly = TRUE)) {
  message("fgsea not installed. Skipping benchmark.")
  q()
}
if (!requireNamespace("enrichit", quietly = TRUE)) {
  message("enrichit not installed.")
  q()
}

library(fgsea)
library(enrichit)
library(ggplot2)

data(geneList, package="DOSE")
stats <- geneList

x = DOSE:::get_dose_data('HDO')
pathways <- split(x@gsid2gene$gene, x@gsid2gene$gsid)


message("Running fgseaMultilevel...")
start_time <- Sys.time()
fgsea_res <- fgseaMultilevel(pathways, stats, minSize=15, maxSize=500, eps=1e-50)
fgsea_time <- Sys.time() - start_time
message("fgsea time: ", fgsea_time)

message("Running enrichit gsea (multilevel)...")
start_time <- Sys.time()
# Note: minPerm in enrichit maps to sampleSize in fgsea roughly?
# refine: fgsea default sampleSize=101. enrichit gsea default minPerm=1000.
# To be fair, set minPerm=101? Or keep defaults?
# Let's try to match params. fgseaMultilevel uses sampleSize=101.
# I'll use minPerm=101 to match.
enrichit_res <- gsea(stats, pathways, method = "multilevel", minPerm = 101, minGSSize=15, maxGSSize=500, pvalThreshold = 1.0, eps = 1e-50) 
enrichit_time <- Sys.time() - start_time
message("enrichit time: ", enrichit_time)

# Compare P-values
# Join results
common <- intersect(fgsea_res$pathway, enrichit_res$ID)
fgsea_sub <- fgsea_res[match(common, fgsea_res$pathway), ]
enrichit_sub <- enrichit_res[match(common, enrichit_res$ID), ]

pval_fgsea <- fgsea_sub$pval
pval_enrichit <- enrichit_sub$pvalue

cor_pval <- cor(pval_fgsea, pval_enrichit, method = "spearman")
message("Correlation of p-values (Spearman): ", cor_pval)

# Print comparison
df <- data.frame(ID = common, FGSEA = pval_fgsea, ENRICHIT = pval_enrichit)
print(head(df))

if (cor_pval > 0.9) {
    message("SUCCESS: P-values are highly correlated.")
} else {
    message("WARNING: Low correlation.")
}
