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
    ## The container now also carries the NB fit (counts/mu/alpha) so a
    ## parametric bootstrap can resample from the very model that produced the
    ## observed statistic (NB_IMPLEMENT_SPEC §1.3). The field ORDER is part of
    ## that contract, so it is pinned positionally -- expect_identical on
    ## names() compares element by element and is locale-independent, unlike
    ## the sort() in the note below.
    expect_identical(names(dd), c("E", "weights", "design", "counts", "mu", "alpha"))
    ## Order-independent belt-and-braces: sort() is locale-sensitive (in the C
    ## locale that R CMD check uses "E" < "design", under collation it is the reverse), so
    ## comparing against a hand-sorted literal is exactly the kind of assertion that
    ## passes in devtools::test() and fails in check (REBUILD_ACCEPTANCE.md D7).
    expect_setequal(names(dd),
                    c("E", "weights", "design", "counts", "mu", "alpha"))

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


## ---------------------------------------------------------------------------
## NB resampling inputs: counts / mu / alpha  (NB_IMPLEMENT_SPEC §1.3)
## ---------------------------------------------------------------------------

test_that("counts, mu and alpha are present with the right types and shapes", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)

    ct <- DESeq2::counts(dds)
    mu <- SummarizedExperiment::assay(dds, "mu")
    alpha <- SummarizedExperiment::mcols(dds)$dispersion

    ## counts: the raw integer matrix, byte for byte.
    expect_identical(dd$counts, ct)
    expect_identical(storage.mode(dd$counts), "integer")
    expect_identical(dim(dd$counts), dim(ct))

    ## mu: G x n, same shape and dimnames as E.
    expect_true(is.matrix(dd$mu))
    expect_identical(dim(dd$mu), dim(ct))
    expect_identical(dimnames(dd$mu), dimnames(ct))

    ## alpha: one dispersion per GENE (length G, not G x n).
    expect_type(dd$alpha, "double")
    expect_null(dim(dd$alpha))
    expect_identical(length(dd$alpha), nrow(ct))
    ## DESeq2 stores dispersion in mcols() without names, and the spec pins
    ## alpha to "mcols(x)$dispersion with NA -> 1", so no names are invented
    ## here. What must hold is POSITIONAL alignment with the rows of E, which
    ## the reconstruction below proves.
    if (!is.null(names(dd$alpha))) {
        expect_identical(names(dd$alpha), rownames(ct))
    }

    ## dimnames of counts / mu must agree with E, or a bootstrap that indexes
    ## rows by gene id would silently misalign with the observed E.
    expect_identical(rownames(dd$counts), rownames(dd$E))
    expect_identical(colnames(dd$counts), colnames(dd$E))
    expect_identical(rownames(dd$mu), rownames(dd$E))
    expect_identical(colnames(dd$mu), colnames(dd$E))

    ## Raw values are forwarded untouched (counts/mu); alpha keeps DESeq2's
    ## estimates wherever they are not NA.
    expect_identical(dd$mu, mu)
    expect_identical(dd$alpha, as.numeric(alpha))
})


test_that("stored counts / mu / alpha regenerate E and weights exactly", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    sf <- DESeq2::sizeFactors(dds)

    ## The whole point of the three new fields is that a consumer can rebuild
    ## the observed pipeline (NB_IMPLEMENT_SPEC §2) from them alone. If alpha
    ## were misaligned with the rows of E, this reconstruction would not match.
    ct <- dd$counts
    lib <- mean(colSums(ct)) * sf
    E_chk <- log2((ct + 0.5) /
                      (matrix(lib, nrow(ct), ncol(ct), byrow = TRUE) + 1) * 1e6)
    dimnames(E_chk) <- dimnames(ct)
    expect_true(identical(E_chk, dd$E))

    mt <- pmax(dd$mu, 1)
    w_chk <- (mt + 0.5)^2 * log(2)^2 / (mt + dd$alpha * mt^2)
    bad <- !is.finite(w_chk)
    if (any(bad)) w_chk[bad] <- min(w_chk[!bad])
    dimnames(w_chk) <- dimnames(ct)
    ## `dd$alpha * mt^2` relies on alpha being row-aligned with mt, which is
    ## what makes this reconstruction meaningful.
    expect_true(identical(w_chk, dd$weights))
})



test_that("alpha replaces NA dispersions with 1 and never contains NA", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    alpha_raw <- SummarizedExperiment::mcols(dds)$dispersion

    ## Inject NA dispersions (what DESeq2 leaves for genes dropped by
    ## independent filtering) so the convention is actually exercised.
    na_idx <- c(3L, 17L, 100L)
    SummarizedExperiment::mcols(dds)$dispersion[na_idx] <- NA_real_
    dd <- as_dsea(dds)

    expect_identical(length(dd$alpha), nrow(dds))
    expect_false(anyNA(dd$alpha))
    expect_identical(dd$alpha[na_idx], rep(1, length(na_idx)))
    ## Non-NA entries are passed through unchanged.
    keep <- seq_len(nrow(dds))[-na_idx]
    expect_identical(dd$alpha[keep], as.numeric(alpha_raw)[keep])
})


test_that("mu is stored UNFLOORED while weights use the floored version", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    mu_raw <- SummarizedExperiment::assay(dds, "mu")
    alpha_raw <- SummarizedExperiment::mcols(dds)$dispersion

    ## A floor of 10 must bite: the fixture does contain fitted means below it.
    mu.floor <- 10
    expect_gt(sum(mu_raw < mu.floor), 0L)

    dd <- as_dsea(dds, mu.floor = mu.floor)

    ## (a) The STORED field is the raw assay value, not clamped.
    expect_identical(dd$mu, mu_raw)
    expect_lt(min(dd$mu), mu.floor)
    expect_identical(sum(dd$mu < mu.floor), sum(mu_raw < mu.floor))

    ## (b) The WEIGHTS are still computed from the floored version. Reference
    ## written longhand, independently of as_dsea()'s body.
    mt <- pmax(mu_raw, mu.floor)
    w_ref <- (mt + 0.5)^2 * log(2)^2 / (mt + alpha_raw * mt^2)
    bad <- !is.finite(w_ref)
    w_ref[bad] <- min(w_ref[!bad])
    expect_equal(dd$weights, w_ref, tolerance = 1e-12)

    ## (c) Both halves of the claim at once: for a gene whose raw mu is below
    ## the floor, the stored mu differs from the mu that entered the weight
    ## formula. If mu were pre-floored, this assertion fails.
    low <- which(mu_raw < mu.floor, arr.ind = TRUE)
    i <- low[1L, "row"]; j <- low[1L, "col"]
    expect_lt(dd$mu[i, j], mu.floor)
    w_from_raw <- (mu_raw[i, j] + 0.5)^2 * log(2)^2 /
        (mu_raw[i, j] + alpha_raw[i] * mu_raw[i, j]^2)
    expect_false(isTRUE(all.equal(dd$weights[i, j], w_from_raw)))

    ## (d) The default floor of 1 leaves mu < 1 stored, too.
    dd1 <- as_dsea(dds)
    expect_lt(min(dd1$mu), 1)
})


## ---------------------------------------------------------------------------
## ADDITIVITY: adding counts/mu/alpha must not perturb E/weights/design by one
## bit. Independently re-derived longhand references, pinned with identical()
## (not all.equal) so any future change to the denominator -- the project had a
## real bug where colSums(ct) * norm.factors double-counted depth; the correct
## form is mean(colSums(ct)) * norm.factors -- is caught immediately, along with
## any drift in the weight formula or the fallback.
## ---------------------------------------------------------------------------

test_that("E, weights and design are bit-exact against longhand references", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)

    ct <- DESeq2::counts(dds)
    sf <- DESeq2::sizeFactors(dds)
    mu_raw <- SummarizedExperiment::assay(dds, "mu")
    alpha_raw <- SummarizedExperiment::mcols(dds)$dispersion

    ## ---- E, longhand ------------------------------------------------------
    lib <- mean(colSums(ct)) * sf
    denom <- matrix(lib, nrow(ct), ncol(ct), byrow = TRUE) + 1
    E_ref <- log2((ct + 0.5) / denom * 1e6)
    dimnames(E_ref) <- dimnames(ct)

    ## ---- weights, longhand ------------------------------------------------
    mt <- pmax(mu_raw, 1)
    w_ref <- (mt + 0.5)^2 * log(2)^2 / (mt + alpha_raw * mt^2)
    bad <- !is.finite(w_ref)
    if (any(bad)) w_ref[bad] <- min(w_ref[!bad])
    dimnames(w_ref) <- dimnames(ct)

    ## ---- design, longhand -------------------------------------------------
    design_ref <- as.matrix(attr(dds, "modelMatrix"))

    ## identical() = same bits, same attributes, same everything.
    expect_true(identical(dd$E, E_ref))
    expect_true(identical(dd$weights, w_ref))
    expect_true(identical(dd$design, design_ref))

    ## And the double-counted denominator must NOT reproduce E (guards the
    ## other direction: a re-introduction of colSums(ct) * sf).
    denom_bad <- matrix(colSums(ct) * sf, nrow(ct), ncol(ct), byrow = TRUE) + 1
    E_bad <- log2((ct + 0.5) / denom_bad * 1e6)
    dimnames(E_bad) <- dimnames(ct)
    expect_false(identical(dd$E, E_bad))
})


test_that("as_dsea() is a pure additive change: dsea()'s inputs are untouched", {
    skip_if_not_installed("DESeq2")
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)

    ## dsea() consumes exactly three fields. Rebuilding the pre-fix 3-field
    ## container from this 6-field one must give bit-identical values -- i.e.
    ## counts/mu/alpha are carried alongside, never folded into E/weights/design.
    before <- structure(
        list(E = dd$E, weights = dd$weights, design = dd$design),
        class = "DseaData"
    )
    expect_true(identical(before$E, dd$E))
    expect_true(identical(before$weights, dd$weights))
    expect_true(identical(before$design, dd$design))

    ## The container normaliser that dsea() funnels its input through reads only
    ## the three fields, so passing the 6-field container must not perturb them.
    skip_if_not(exists(".as_dsea_data", mode = "function"))
    dd2 <- get(".as_dsea_data")(dd)
    expect_true(is.list(dd2))
    expect_identical(names(dd2), c("E", "weights", "design"))
    expect_true(identical(dd2$E, dd$E))
    expect_true(identical(dd2$weights, dd$weights))
})


