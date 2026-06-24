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

  agg <- suppressWarnings(
    aggregate_enrichment(list(rna = res1, prot = res2), method = "fisher")
  )

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

  agg <- suppressWarnings(
    aggregate_enrichment(list(res1, res2), method = "fisher")
  )

  expect_s4_class(agg, "enrichResult")
  expect_equal(sort(agg@result$ID), c("Path1", "Path2"))
  expect_equal(agg@ontology, "TEST")
  expect_equal(agg@organism, "human")
})

test_that("aggregate_omics supports Brown method with explicit covariance", {
  p_mat <- matrix(
    c(0.01, 0.02,
      0.20, 0.05),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2"), c("rna", "prot"))
  )
  cov_mat <- matrix(c(4, 1, 1, 4), nrow = 2)

  res <- aggregate_omics(
    p_mat,
    method = "brown",
    input = "pvalue",
    cov_matrix = cov_mat
  )

  expect_s3_class(res, "omics_aggregated")
  expect_equal(names(res$pvalue), c("Gene1", "Gene2"))

  x1 <- -2 * sum(log(c(0.01, 0.02)))
  x2 <- -2 * sum(log(c(0.20, 0.05)))
  expected <- stats::pchisq(c(x1, x2) / 1.25, df = 3.2, lower.tail = FALSE)

  expect_equal(unname(res$pvalue), expected, tolerance = 1e-12)
  expect_equal(
    unname(res$score),
    -log10(expected),
    tolerance = 1e-12
  )
})

test_that("aggregate_omics handles signed conflicts and weighted mean", {
  score_mat <- matrix(
    c(1, 3,
      2, -4,
      -2, -6),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2", "Gene3"), c("rna", "prot"))
  )

  strict_res <- aggregate_omics(
    score_mat,
    method = "mean",
    input = "signed_score",
    conflict_policy = "strict"
  )
  expect_equal(strict_res$feature_id, c("Gene1", "Gene3"))
  expect_equal(unname(strict_res$score), c(2, -4))

  penalty_res <- aggregate_omics(
    score_mat,
    method = "mean",
    input = "signed_score",
    conflict_policy = "penalty"
  )
  expect_equal(unname(penalty_res$score), c(2, -0.5, -4))

  weighted_res <- aggregate_omics(
    score_mat,
    method = "weighted_mean",
    input = "signed_score",
    weights = c(1, 3)
  )
  expect_equal(unname(weighted_res$score), c(2.5, -2.5, -5), tolerance = 1e-12)
})

test_that("harmonize_ids and select_features_for_ora keep the workflow lightweight", {
  agg <- structure(
    list(
      input_type = "pvalue",
      feature_type = "protein",
      feature_id = c("P1", "P2", "P3"),
      score = c(P1 = 4, P2 = 2, P3 = 1),
      pvalue = c(P1 = 1e-4, P2 = 1e-2, P3 = 0.20)
    ),
    class = "omics_aggregated"
  )

  mapping <- data.frame(
    source_id = c("P1", "P2", "P3"),
    target_id = c("G1", "G1", "G2"),
    stringsAsFactors = FALSE
  )

  harm <- harmonize_ids(agg, mapping, collapse = "min_p")

  expect_s3_class(harm, "omics_aggregated")
  expect_equal(harm$feature_type, "gene")
  expect_equal(harm$feature_id, c("G1", "G2"))
  expect_equal(unname(harm$pvalue), c(1e-4, 0.20))
  expect_equal(unname(harm$score), c(4, 1))

  ora_input <- select_features_for_ora(harm, cutoff = 0.05, by = "pvalue")
  expect_equal(ora_input$gene, "G1")
  expect_equal(ora_input$universe, c("G1", "G2"))
})

test_that("get_omics_contribution extracts pathway-level source statistics", {
  agg <- structure(
    list(
      input_type = "pvalue",
      feature_type = "gene",
      feature_id = c("GeneA", "GeneB", "GeneC"),
      score = c(GeneA = 4, GeneB = 2, GeneC = 1),
      pvalue = c(GeneA = 1e-4, GeneB = 1e-2, GeneC = 0.2),
      original_matrix = matrix(
        c(1e-4, 1e-2,
          1e-2, 0.2,
          0.2, 0.3),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(c("GeneA", "GeneB", "GeneC"), c("rna", "prot"))
      )
    ),
    class = "omics_aggregated"
  )

  res_df <- data.frame(
    ID = c("Path1", "Path2"),
    Description = c("Path 1", "Path 2"),
    pvalue = c(0.001, 0.2),
    p.adjust = c(0.002, 0.2),
    qvalue = c(0.002, 0.2),
    geneID = c("GeneA/GeneB", "GeneC"),
    Count = c(2, 1),
    stringsAsFactors = FALSE
  )
  res <- make_enrich_result_for_test(res_df, pvalue_cutoff = 1)

  contrib <- get_omics_contribution(res, agg, pathway_id = "Path1")

  expect_true(is.data.frame(contrib))
  expect_equal(contrib$Feature, c("GeneA", "GeneB"))
  expect_true(all(c("rna", "prot", "Aggregated_Score", "Aggregated_Pvalue") %in% colnames(contrib)))
  expect_equal(contrib$Aggregated_Score, c(4, 2))
})

test_that("classify_omics_pattern labels merged pathways by single-omics support", {
  merged_df <- data.frame(
    ID = c("PathEnhanced", "PathRNA", "PathShared", "PathNS"),
    Description = c("Enhanced", "RNA only", "Shared", "Not sig"),
    pvalue = c(0.01, 0.01, 0.01, 0.2),
    p.adjust = c(0.01, 0.01, 0.01, 0.2),
    qvalue = c(0.01, 0.01, 0.01, 0.2),
    geneID = c("GeneA", "GeneB", "GeneC", "GeneD"),
    Count = c(1, 1, 1, 1),
    stringsAsFactors = FALSE
  )
  rna_df <- data.frame(
    ID = c("PathRNA", "PathShared"),
    Description = c("RNA only", "Shared"),
    pvalue = c(0.01, 0.02),
    p.adjust = c(0.01, 0.02),
    qvalue = c(0.01, 0.02),
    geneID = c("GeneB", "GeneC"),
    Count = c(1, 1),
    stringsAsFactors = FALSE
  )
  prot_df <- data.frame(
    ID = c("PathShared"),
    Description = c("Shared"),
    pvalue = c(0.03),
    p.adjust = c(0.03),
    qvalue = c(0.03),
    geneID = c("GeneC"),
    Count = c(1),
    stringsAsFactors = FALSE
  )

  merged <- make_enrich_result_for_test(merged_df, pvalue_cutoff = 1)
  rna <- make_enrich_result_for_test(rna_df, pvalue_cutoff = 1)
  prot <- make_enrich_result_for_test(prot_df, pvalue_cutoff = 1)

  classified <- classify_omics_pattern(
    merged_res = merged,
    single_res = list(rna = rna, prot = prot),
    p_cutoff = 0.05,
    by = "p.adjust"
  )

  out <- classified@result$Omics_Pattern
  names(out) <- classified@result$ID

  expect_equal(unname(out["PathEnhanced"]), "Enhanced (1+1>2)")
  expect_equal(unname(out["PathRNA"]), "rna-Specific")
  expect_equal(unname(out["PathShared"]), "Shared")
  expect_equal(unname(out["PathNS"]), "Not Significant")
})
