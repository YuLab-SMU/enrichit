library(testthat)
library(enrichit)

make_enrich_result_accessor_test <- function() {
  result_df <- data.frame(
    ID = c("Path1", "Path2"),
    Description = c("Path 1", "Path 2"),
    GeneRatio = c("2/3", "1/3"),
    BgRatio = c("2/10", "4/10"),
    RichFactor = c(1, 0.25),
    FoldEnrichment = c(5, 1.25),
    pvalue = c(0.001, 0.02),
    p.adjust = c(0.002, 0.03),
    qvalue = c(0.002, 0.03),
    geneID = c("GeneA/GeneB", "GeneC"),
    Count = c(2, 1),
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- result_df$ID

  new("enrichResult",
    result = result_df,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    organism = "human",
    ontology = "TEST",
    gene = c("GeneA", "GeneB", "GeneC"),
    keytype = "SYMBOL",
    universe = paste0("Gene", LETTERS[1:10]),
    gene2Symbol = character(0),
    geneSets = list(Path1 = c("GeneA", "GeneB"), Path2 = "GeneC"),
    readable = FALSE,
    termsim = matrix(0, nrow = 0, ncol = 0),
    method = "ORA",
    dr = list()
  )
}

make_gsea_result_accessor_test <- function() {
  result_df <- data.frame(
    ID = c("Top", "Bottom"),
    Description = c("Top pathway", "Bottom pathway"),
    setSize = c(2L, 2L),
    enrichmentScore = c(0.7, -0.6),
    NES = c(1.8, -1.5),
    pvalue = c(0.01, 0.03),
    p.adjust = c(0.02, 0.03),
    qvalue = c(0.02, 0.03),
    rank = c(2L, 5L),
    leading_edge = c("tags=100%, list=40%, signal=70%", "tags=100%, list=40%, signal=70%"),
    core_enrichment = c("GeneA/GeneB", "GeneD/GeneE"),
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- result_df$ID

  new("gseaResult",
    result = result_df,
    organism = "human",
    setType = "TEST",
    geneSets = list(Top = c("GeneA", "GeneB"), Bottom = c("GeneD", "GeneE")),
    geneList = c(GeneA = 3, GeneB = 2, GeneC = 1, GeneD = -2, GeneE = -3),
    keytype = "SYMBOL",
    permScores = matrix(0, nrow = 0, ncol = 0),
    params = list(pvalueCutoff = 1, pAdjustMethod = "BH"),
    gene2Symbol = character(0),
    readable = FALSE,
    termsim = matrix(0, nrow = 0, ncol = 0),
    method = "GSEA",
    dr = list()
  )
}

test_that("geneID and geneInCategory return expected feature mappings", {
  enrich_res <- make_enrich_result_accessor_test()
  gsea_res <- make_gsea_result_accessor_test()

  expect_equal(geneID(enrich_res), c("GeneA/GeneB", "GeneC"))
  expect_equal(geneID(gsea_res), c("GeneA/GeneB", "GeneD/GeneE"))

  enrich_cat <- geneInCategory(enrich_res)
  gsea_cat <- geneInCategory(gsea_res)

  expect_equal(enrich_cat$Path1, c("GeneA", "GeneB"))
  expect_equal(gsea_cat$Bottom, c("GeneD", "GeneE"))
})

test_that("gsfilter filters enrichResult by Count and GSSize", {
  enrich_res <- make_enrich_result_accessor_test()

  by_count <- gsfilter(enrich_res, by = "Count", min = 2)
  expect_s4_class(by_count, "enrichResult")
  expect_equal(by_count@result$ID, "Path1")

  by_size <- gsfilter(enrich_res, by = "GSSize", max = 2)
  expect_equal(by_size@result$ID, "Path1")
})

test_that("summary delegates to as.data.frame for enrichResult and gseaResult", {
  enrich_res <- make_enrich_result_accessor_test()
  gsea_res <- make_gsea_result_accessor_test()

  expect_warning(
    expect_equal(summary(enrich_res), as.data.frame(enrich_res)),
    "deprecated"
  )
  expect_warning(
    expect_equal(summary(gsea_res), as.data.frame(gsea_res)),
    "deprecated"
  )
})
