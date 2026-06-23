library(testthat)
library(enrichit)

make_enrich_result_for_test <- function(result_df, pvalue_cutoff = 0.05) {
  all_genes <- unique(unlist(strsplit(result_df$geneID, "/", fixed = TRUE), use.names = FALSE))
  gene_sets <- setNames(
    strsplit(result_df$geneID, "/", fixed = TRUE),
    result_df$ID
  )

  new("enrichResult",
    result = result_df,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    organism = "human",
    ontology = "TEST",
    gene = all_genes,
    keytype = "SYMBOL",
    universe = all_genes,
    gene2Symbol = character(0),
    geneSets = gene_sets,
    readable = FALSE,
    termsim = matrix(0, nrow = 0, ncol = 0),
    method = "ORA",
    dr = list()
  )
}

test_that("aggregate_enrichment combines shared pathways with Fisher method", {
  res1 <- data.frame(
    ID = c("Path1", "Path2"),
    Description = c("Path 1", "Path 2"),
    pvalue = c(0.001, 0.20),
    geneID = c("GeneA/GeneB", "GeneC"),
    Count = c(2, 1),
    stringsAsFactors = FALSE
  )
  res2 <- data.frame(
    ID = c("Path1", "Path2"),
    Description = c("Path 1", "Path 2"),
    pvalue = c(0.01, 0.03),
    core_enrichment = c("GeneB/GeneD", "GeneC/GeneE"),
    stringsAsFactors = FALSE
  )

  agg <- aggregate_enrichment(list(rna = res1, prot = res2), method = "fisher")

  expect_s4_class(agg, "enrichResult")
  expect_equal(agg@result$ID, c("Path1", "Path2"))

  expected_p <- c(
    stats::pchisq(-2 * sum(log(c(0.001, 0.01))), df = 4, lower.tail = FALSE),
    stats::pchisq(-2 * sum(log(c(0.20, 0.03))), df = 4, lower.tail = FALSE)
  )
  expect_equal(agg@result$pvalue, expected_p, tolerance = 1e-12)
  expect_equal(agg@result$Count, c(3L, 2L))
  expect_equal(agg@result$geneID[agg@result$ID == "Path1"], "GeneA/GeneB/GeneD")
  expect_equal(agg@result$geneID[agg@result$ID == "Path2"], "GeneC/GeneE")
})

test_that("aggregate_enrichment reads raw enrichResult tables instead of filtered as.data.frame output", {
  res1_df <- data.frame(
    ID = c("Path1", "Path2"),
    Description = c("Path 1", "Path 2"),
    pvalue = c(0.001, 0.20),
    p.adjust = c(0.002, 0.20),
    qvalue = c(0.002, 0.20),
    geneID = c("GeneA/GeneB", "GeneC"),
    Count = c(2, 1),
    stringsAsFactors = FALSE
  )
  res2_df <- data.frame(
    ID = c("Path1", "Path2"),
    Description = c("Path 1", "Path 2"),
    pvalue = c(0.02, 0.03),
    p.adjust = c(0.02, 0.03),
    qvalue = c(0.02, 0.03),
    geneID = c("GeneB/GeneD", "GeneC/GeneE"),
    Count = c(2, 2),
    stringsAsFactors = FALSE
  )

  res1 <- make_enrich_result_for_test(res1_df, pvalue_cutoff = 0.05)
  res2 <- make_enrich_result_for_test(res2_df, pvalue_cutoff = 0.05)

  expect_equal(as.data.frame(res1)$ID, "Path1")

  agg <- aggregate_enrichment(list(res1, res2), method = "fisher")

  expect_s4_class(agg, "enrichResult")
  expect_equal(sort(agg@result$ID), c("Path1", "Path2"))
  expect_equal(agg@ontology, "TEST")
  expect_equal(agg@organism, "human")
})
