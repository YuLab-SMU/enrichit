library(testthat)
library(enrichit)

# Define a helper function to compute expected p-values using phyper
expected_pvalue <- function(set_size, overlap, de_total, universe_size) {
  # hypergeometric: probability of >= overlap successes
  # q = overlap - 1, lower.tail = FALSE
  phyper(overlap - 1, set_size, universe_size - set_size, de_total, lower.tail = FALSE)
}

test_that("ORA function returns correct p-values and format", {
  # Example data from documentation
  de_genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
  all_genes <- paste0("Gene", 1:1000)
  gs1 <- paste0("Gene", 1:50)
  gs2 <- paste0("Gene", 51:150)
  gs3 <- paste0("Gene", 151:300)
  gene_sets <- list(Pathway1 = gs1, Pathway2 = gs2, Pathway3 = gs3)

  result <- ora(gene = de_genes, gene_sets = gene_sets, universe = all_genes)

  # Verify columns exist
  expected_cols <- c("ID", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment", "pvalue", "geneID", "Count")
  expect_true(all(expected_cols %in% colnames(result)))

  # Manual calculations
  manual_p1 <- expected_pvalue(set_size = 50, overlap = 5, de_total = 5, universe_size = 1000)
  manual_p2 <- expected_pvalue(set_size = 100, overlap = 0, de_total = 5, universe_size = 1000)
  manual_p3 <- expected_pvalue(set_size = 150, overlap = 0, de_total = 5, universe_size = 1000)

  # Compare with ORA output (allow tiny numerical tolerance)
  expect_equal(result$pvalue[result$ID == "Pathway1"], manual_p1, tolerance = 1e-10)
  expect_equal(result$pvalue[result$ID == "Pathway2"], manual_p2, tolerance = 1e-10)
  expect_equal(result$pvalue[result$ID == "Pathway3"], manual_p3, tolerance = 1e-10)
  
  # Check content of new columns for Pathway1
  p1_res <- result[result$ID == "Pathway1", ]
  expect_equal(p1_res$Count, 5)
  expect_equal(p1_res$GeneRatio, "5/5")
  expect_equal(p1_res$BgRatio, "50/1000")
  expect_equal(p1_res$RichFactor, 5/50)
  expect_equal(p1_res$FoldEnrichment, (5/5) / (50/1000))
  
  # Check that all expected genes are present in geneID
  genes_in_id <- strsplit(p1_res$geneID, "/")[[1]]
  expect_true(all(c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5") %in% genes_in_id))
  expect_equal(length(genes_in_id), 5)
})

test_that("ora_gson excludes zero-overlap sets from multiple testing correction", {
  gsid2gene <- data.frame(
    gsid = c("set1", "set1", "set2", "set2", "set3", "set3"),
    gene = c("Gene1", "Gene2", "Gene2", "Gene3", "Gene4", "Gene5"),
    stringsAsFactors = FALSE
  )
  gsid2name <- data.frame(
    gsid = c("set1", "set2", "set3"),
    name = c("Set 1", "Set 2", "Set 3"),
    stringsAsFactors = FALSE
  )
  gson_obj <- gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = "test",
    gsname = "test",
    version = "test",
    accessed_date = as.character(Sys.Date()),
    keytype = "UNKNOWN"
  )

  res <- ora_gson(
    gene = c("Gene1", "Gene2"),
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5"),
    minGSSize = 1,
    maxGSSize = 10,
    qvalueCutoff = 1,
    gson = gson_obj
  )

  expect_s4_class(res, "enrichResult")
  expect_equal(res@result$ID, c("set1", "set2"))
  expect_true(all(res@result$Count > 0))

  raw_p <- c(
    expected_pvalue(set_size = 2, overlap = 2, de_total = 2, universe_size = 5),
    expected_pvalue(set_size = 2, overlap = 1, de_total = 2, universe_size = 5)
  )
  expect_equal(res@result$pvalue, raw_p, tolerance = 1e-10)
  expect_equal(res@result$p.adjust, p.adjust(raw_p, method = "BH"), tolerance = 1e-10)
})

test_that("weighted ORA runs on a small universe", {
  skip_if_not_installed("BiasedUrn")

  de_genes <- c("Gene1", "Gene2", "Gene3")
  all_genes <- paste0("Gene", 1:30)
  gene_sets <- list(
    Pathway1 = paste0("Gene", 1:6),
    Pathway2 = paste0("Gene", 10:18)
  )
  weight <- setNames(rep(1, 30), all_genes)
  weight[c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6")] <- 3

  unweighted <- ora(gene = de_genes, gene_sets = gene_sets, universe = all_genes)
  weighted <- ora(gene = de_genes, gene_sets = gene_sets, universe = all_genes, weight = weight)

  expect_true(is.data.frame(weighted))
  expect_true(all(c("ID", "pvalue", "Count") %in% colnames(weighted)))
  expect_false(isTRUE(all.equal(weighted$pvalue, unweighted$pvalue)))
})
