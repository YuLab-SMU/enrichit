library(testthat)
library(enrichit)

# Define a helper function to compute expected p-values using phyper
expected_pvalue <- function(set_size, overlap, de_total, universe_size) {
  # hypergeometric: probability of >= overlap successes
  # q = overlap - 1, lower.tail = FALSE
  phyper(overlap - 1, set_size, universe_size - set_size, de_total, lower.tail = FALSE)
}

test_that("ORA function returns correct p-values", {
  # Example data from documentation
  de_genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
  all_genes <- paste0("Gene", 1:1000)
  gs1 <- paste0("Gene", 1:50)
  gs2 <- paste0("Gene", 51:150)
  gs3 <- paste0("Gene", 151:300)
  gene_sets <- list(Pathway1 = gs1, Pathway2 = gs2, Pathway3 = gs3)

  result <- ora(gene = de_genes, gene_sets = gene_sets, universe = all_genes)

  # Verify columns exist
  expect_true(all(c("GeneSet", "SetSize", "DEInSet", "DESize", "PValue") %in% colnames(result)))

  # Manual calculations
  manual_p1 <- expected_pvalue(set_size = 50, overlap = 5, de_total = 5, universe_size = 1000)
  manual_p2 <- expected_pvalue(set_size = 100, overlap = 0, de_total = 5, universe_size = 1000)
  manual_p3 <- expected_pvalue(set_size = 150, overlap = 0, de_total = 5, universe_size = 1000)

  # Compare with ORA output (allow tiny numerical tolerance)
  expect_equal(result$PValue[result$GeneSet == "Pathway1"], manual_p1, tolerance = 1e-10)
  expect_equal(result$PValue[result$GeneSet == "Pathway2"], manual_p2, tolerance = 1e-10)
  expect_equal(result$PValue[result$GeneSet == "Pathway3"], manual_p3, tolerance = 1e-10)
})
