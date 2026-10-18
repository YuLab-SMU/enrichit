
test_that("filter_gene_sets is identical to per-set intersect", {
    set.seed(1)
    universe <- paste0("G", 1:200)
    sets <- list(
        a = c("G1", "G2", "G2", "G500"),
        b = character(),
        c = c("G3", "G200", NA),
        d = c("G7", "G7", "G8")
    )
    old <- lapply(sets, intersect, universe)
    new <- enrichit:::filter_gene_sets(sets, universe)
    expect_identical(new, old)
    expect_identical(names(new), names(sets))
    expect_identical(enrichit:::filter_gene_sets(list(), universe), list())
})
