## Tests for the DseaData container built by as_dsea().
##
## These tests pin down REBUILD_SPEC section 2 ("fix 3"): the library size used
## to build E must be mean(colSums(ct)) * norm.factors, NOT colSums(ct) *
## norm.factors. The latter counts sequencing depth twice and leaves an
## unmodellable -log2(colSums_j) per-sample offset in E.

## ---------------------------------------------------------------------------
## Fixtures: a small, fast, fully-fitted DESeqDataSet.
## ---------------------------------------------------------------------------

make_fitted_dds <- function() {
    testthat::skip_if_not_installed("DESeq2")
    set.seed(42)
    dds <- DESeq2::makeExampleDESeqDataSet(n = 200, m = 8)
    suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
}

make_unfitted_dds <- function() {
    testthat::skip_if_not_installed("DESeq2")
    set.seed(42)
    DESeq2::makeExampleDESeqDataSet(n = 100, m = 6)
}

test_that("as_dsea() is an S3 generic with a DESeqDataSet method", {
    expect_true(is.function(as_dsea))
    expect_true(is.function(getS3method("as_dsea", "DESeqDataSet", optional = TRUE)))
})

test_that("as_dsea() returns a DseaData list with E, weights and design", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)

    expect_s3_class(dd, "DseaData")
    expect_true(is.list(dd))
    ## Order-independent on purpose: sort() is locale-sensitive (in the C locale
    ## that R CMD check uses "E" < "design", under collation it is the reverse), so
    ## comparing against a hand-sorted literal is exactly the kind of assertion that
    ## passes in devtools::test() and fails in check (REBUILD_ACCEPTANCE.md D7).
    expect_setequal(names(dd), c("E", "weights", "design"))

    ct <- DESeq2::counts(dds)
    expect_identical(dim(dd$E), dim(ct))
    expect_identical(dim(dd$weights), dim(ct))
    expect_identical(rownames(dd$E), rownames(ct))
    expect_identical(colnames(dd$E), colnames(ct))
    ## design is samples x coefficients
    expect_identical(nrow(dd$design), ncol(ct))
})


## ---------------------------------------------------------------------------
## The E formula, checked bit-exactly against an independently computed
## reference. Tolerance 1e-12 (the two expressions are algebraically identical
## but not written identically, so allow rounding, not a different formula).
## ---------------------------------------------------------------------------

test_that("E uses mean(colSums(ct)) * norm.factors as the library size", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    ct <- DESeq2::counts(dds)
    sf <- DESeq2::sizeFactors(dds)

    ## Reference written out longhand, independently of as_dsea()'s internals.
    lib <- mean(colSums(ct)) * sf
    denom <- matrix(lib, nrow(ct), ncol(ct), byrow = TRUE) + 1
    E_ref <- log2((ct + 0.5) / denom * 1e6)

    dd <- as_dsea(dds)
    expect_equal(dd$E, E_ref, tolerance = 1e-12)
})


## ---------------------------------------------------------------------------
## Regression: the double-counted denominator must be caught.
## ---------------------------------------------------------------------------

test_that("E does NOT use the double-counted colSums(ct) * norm.factors", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    ct <- DESeq2::counts(dds)
    sf <- DESeq2::sizeFactors(dds)

    ## The pre-fix (wrong) denominator.
    lib_bad <- colSums(ct) * sf
    denom_bad <- matrix(lib_bad, nrow(ct), ncol(ct), byrow = TRUE) + 1
    E_bad <- log2((ct + 0.5) / denom_bad * 1e6)

    dd <- as_dsea(dds)

    ## The two formulas must NOT agree. Their difference is exactly
    ## log2(colSums_j / mean(colSums)) -- a per-sample offset (up to the small
    ## +1 in the denominator, hence an approximate rather than exact match),
    ## so it is constant down every column and non-zero for every sample.
    expect_false(isTRUE(all.equal(dd$E, E_bad, tolerance = 1e-6)))

    offset <- colMeans(dd$E - E_bad)
    expect_equal(offset, log2(colSums(ct) / mean(colSums(ct))), tolerance = 1e-4)
    expect_true(all(abs(offset) > 1e-3))

    ## And the per-column difference is constant (it is a sample-level effect,
    ## not a gene-level one) -- this is what makes it unmodellable by a design
    ## that is constant across genes.
    spread <- apply(dd$E - E_bad, 2, function(z) max(z) - min(z))
    expect_true(all(spread < 1e-3))
})


test_that("mean(colSums) convention removes the per-sample depth offset from E", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    ct <- DESeq2::counts(dds)
    dd <- as_dsea(dds)

    ## With the corrected denominator, colMeans(E) should still track the
    ## (negated) depth, but the residual per-sample offset left after removing
    ## the design should be much smaller than with the double-counted version.
    res_sd <- function(E) {
        X <- dd$design
        R <- t(lm.fit(X, t(E))$residuals)
        stats::sd(colMeans(R))
    }

    sf <- DESeq2::sizeFactors(dds)
    lib_bad <- colSums(ct) * sf
    E_bad <- log2((ct + 0.5) /
                      (matrix(lib_bad, nrow(ct), ncol(ct), byrow = TRUE) + 1) * 1e6)

    expect_lt(res_sd(dd$E), res_sd(E_bad))
})


## ---------------------------------------------------------------------------
## norm.factors / mu.floor / design
## ---------------------------------------------------------------------------

test_that("norm.factors = NULL uses raw mean depth only", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    ct <- DESeq2::counts(dds)

    dd <- as_dsea(dds, norm.factors = NULL)
    lib <- rep(mean(colSums(ct)), ncol(ct))
    E_ref <- log2((ct + 0.5) /
                      (matrix(lib, nrow(ct), ncol(ct), byrow = TRUE) + 1) * 1e6)
    expect_equal(dd$E, E_ref, tolerance = 1e-12)
})


test_that("norm.factors of all ones equals the size-factor-free depth", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    dd <- as_dsea(dds, norm.factors = rep(1, ncol(DESeq2::counts(dds))))
    dd_null <- as_dsea(dds, norm.factors = NULL)
    expect_equal(dd$E, dd_null$E, tolerance = 1e-12)
})


test_that("mu.floor is validated", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "a", NULL)) {
        expect_error(as_dsea(dds, mu.floor = bad), "mu.floor")
    }
    expect_s3_class(as_dsea(dds, mu.floor = 0.5), "DseaData")
})


test_that("weights follow the delta-method formula", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()

    mu <- SummarizedExperiment::assay(dds, "mu")
    alpha <- SummarizedExperiment::mcols(dds)$dispersion
    mt <- pmax(mu, 1)
    w_ref <- (mt + 0.5)^2 * log(2)^2 / (mt + alpha * mt^2)
    bad <- !is.finite(w_ref)
    w_ref[bad] <- min(w_ref[!bad])

    dd <- as_dsea(dds)
    expect_equal(dd$weights, w_ref, tolerance = 1e-12)
    expect_true(all(is.finite(dd$weights)))
    expect_true(all(dd$weights > 0))
})


test_that("mu.floor changes the weights as expected", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    w1 <- as_dsea(dds, mu.floor = 1)$weights
    w10 <- as_dsea(dds, mu.floor = 10)$weights
    ## A larger floor can only raise pmax(mu, floor), which raises the numerator
    ## relative to the alpha term for low-mu genes; the matrices must differ.
    expect_false(isTRUE(all.equal(w1, w10)))
})


test_that("as_dsea() errors clearly when DESeq() has not been run", {
    skip_if_not_installed("DESeq2")
    dds <- make_unfitted_dds()
    expect_error(as_dsea(dds), "DESeq")
})


test_that("design comes from the model matrix attribute when present", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    mm <- attr(dds, "modelMatrix")
    skip_if(is.null(mm))
    expect_equal(unname(dd$design), unname(as.matrix(mm)))
    expect_false(is.null(colnames(dd$design)))
})
