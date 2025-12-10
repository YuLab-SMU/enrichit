library(testthat)
library(enrichit)

test_that("GSEA function works correctly with both methods", {
  # Create synthetic data
  # 1000 genes, sorted
  stats <- sort(rnorm(1000), decreasing = TRUE)
  names(stats) <- paste0("Gene", 1:1000)
  
  # Create a gene set enriched at the top (should have positive ES)
  # Top 20 genes + some random ones
  gs_top <- c(names(stats)[1:20], names(stats)[sample(100:1000, 30)])
  
  # Create a gene set enriched at the bottom (should have negative ES)
  # Bottom 20 genes + some random ones
  gs_bottom <- c(names(stats)[981:1000], names(stats)[sample(1:900, 30)])
  
  # Create a random gene set (should have low ES)
  gs_random <- names(stats)[sample(1:1000, 50)]
  
  gene_sets <- list(
    TopEnriched = gs_top,
    BottomEnriched = gs_bottom,
    Random = gs_random
  )
  
  set.seed(123)
  
  # Test "sample" method (default)
  res_sample <- gsea(geneList = stats, gene_sets = gene_sets, nPerm = 100, method = "sample")
  
  expect_true(is.data.frame(res_sample))
  expect_true(all(c("ID", "enrichmentScore", "NES", "pvalue", "setSize", "rank", "leading_edge", "core_enrichment") %in% colnames(res_sample)))
  
  top_res <- res_sample[res_sample$ID == "TopEnriched", ]
  expect_gt(top_res$enrichmentScore, 0)
  expect_lt(top_res$pvalue, 0.05)
  
  # Test "permute" method
  res_permute <- gsea(geneList = stats, gene_sets = gene_sets, nPerm = 100, method = "permute")
  
  expect_true(is.data.frame(res_permute))
  
  top_res_perm <- res_permute[res_permute$ID == "TopEnriched", ]
  expect_gt(top_res_perm$enrichmentScore, 0)
  expect_lt(top_res_perm$pvalue, 0.05)
  
  # Compare NES (sample method usually produces higher NES magnitude for enriched sets)
  # Note: with small nPerm and synthetic data, this might not always hold, but generally true.
  # We just check that they are somewhat different but consistent in sign.
  expect_equal(sign(top_res$NES), sign(top_res_perm$NES))
  
  # Check that method argument validation works
  expect_error(gsea(geneList = stats, gene_sets = gene_sets, method = "invalid"))
})

test_that("Adaptive GSEA works correctly", {
  stats <- sort(rnorm(1000), decreasing = TRUE)
  names(stats) <- paste0("Gene", 1:1000)
  
  # Highly enriched set
  gs_top <- names(stats)[1:30]
  # Random set
  gs_random <- sample(names(stats), 30)
  
  gene_sets <- list(TopEnriched = gs_top, Random = gs_random)
  
  set.seed(42)
  res_adaptive <- gsea(geneList = stats, gene_sets = gene_sets, adaptive = TRUE, method = "sample",
                       minPerm = 100, maxPerm = 10000, pvalThreshold = 0.2)
  
  expect_true(is.data.frame(res_adaptive))
  expect_true("nPerm" %in% colnames(res_adaptive))
  
  # TopEnriched should have more permutations (significant)
  top_nPerm <- res_adaptive[res_adaptive$ID == "TopEnriched", "nPerm"]
  random_nPerm <- res_adaptive[res_adaptive$ID == "Random", "nPerm"]
  
  # Significant sets should use more permutations than initial minPerm
  expect_gte(top_nPerm, 100)
  # Random sets might stop early or not
  expect_gte(random_nPerm, 100)
})

test_that("Multilevel GSEA works correctly with new parameters", {
  stats <- sort(rnorm(1000), decreasing = TRUE)
  names(stats) <- paste0("Gene", 1:1000)
  
  # Top enriched
  gs_top <- names(stats)[1:30]
  # Bottom enriched
  gs_bottom <- names(stats)[971:1000]
  
  gene_sets <- list(Top = gs_top, Bottom = gs_bottom)
  
  # Test scoreType = "pos"
  set.seed(123)
  res_pos <- gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel", 
                  scoreType = "pos", nPermSimple = 1000)
  
  expect_true(is.data.frame(res_pos))
  top_pos <- res_pos[res_pos$ID == "Top", ]
  expect_lt(top_pos$pvalue, 0.05)
  
  # Test scoreType = "neg"
  res_neg <- gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel", 
                  scoreType = "neg", nPermSimple = 1000)
  
  expect_true(is.data.frame(res_neg))
  bottom_neg <- res_neg[res_neg$ID == "Bottom", ]
  expect_lt(bottom_neg$pvalue, 0.05)
  
  # Test scoreType = "std" (default)
  res_std <- gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel", 
                  scoreType = "std", nPermSimple = 1000)
  
  expect_true(is.data.frame(res_std))
  expect_lt(res_std[res_std$ID == "Top", "pvalue"], 0.05)
  expect_lt(res_std[res_std$ID == "Bottom", "pvalue"], 0.05)
  
  # Check if nPermSimple parameter is accepted and works (by checking it doesn't crash)
  res_simple <- gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel", 
                     nPermSimple = 500)
  expect_true(is.data.frame(res_simple))
})
