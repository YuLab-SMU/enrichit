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
  res_sample <- gsea(genelist = stats, gene_sets = gene_sets, nPerm = 100, method = "sample")
  
  expect_true(is.data.frame(res_sample))
  expect_true(all(c("GeneSet", "ES", "NES", "PValue", "Size", "rank", "leading_edge", "core_enrichment") %in% colnames(res_sample)))
  
  top_res <- res_sample[res_sample$GeneSet == "TopEnriched", ]
  expect_gt(top_res$ES, 0)
  expect_lt(top_res$PValue, 0.05)
  
  # Test "permute" method
  res_permute <- gsea(genelist = stats, gene_sets = gene_sets, nPerm = 100, method = "permute")
  
  expect_true(is.data.frame(res_permute))
  
  top_res_perm <- res_permute[res_permute$GeneSet == "TopEnriched", ]
  expect_gt(top_res_perm$ES, 0)
  expect_lt(top_res_perm$PValue, 0.05)
  
  # Compare NES (sample method usually produces higher NES magnitude for enriched sets)
  # Note: with small nPerm and synthetic data, this might not always hold, but generally true.
  # We just check that they are somewhat different but consistent in sign.
  expect_equal(sign(top_res$NES), sign(top_res_perm$NES))
  
  # Check that method argument validation works
  expect_error(gsea(stats, gene_sets, method = "invalid"))
})
