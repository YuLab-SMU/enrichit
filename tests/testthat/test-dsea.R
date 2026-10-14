## Tests for dsea(): the three native set-level tests, their argument surface
## and the golden reference that pins their numerical output.
##
## REBUILD_SPEC section 5.3.  Everything statistical lives in src/set_test.cpp;
## these tests exercise it through the public dsea() entry point only.
##
## Two things this file deliberately does NOT do:
##   * compare against limma.  limma was removed from the package on purpose
##     (REBUILD_SPEC section 1), so parity is checked outside the package (see
##     dsea-benchmark/REBUILD_ACCEPTANCE.md) and inside the package by the
##     golden reference below.
##   * hard-code a formula as the expected value without saying which formula it
##     is.  The E-matrix test writes its reference out longhand and says, in a
##     comment, which of the two candidate denominators it is.

## ---------------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------------

## A small, fast, fully-fitted DESeqDataSet.  200 genes keeps DESeq() at ~1 s.
make_fitted_dds <- function(n = 200, m = 8) {
    testthat::skip_if_not_installed("DESeq2")
    set.seed(42)
    dds <- DESeq2::makeExampleDESeqDataSet(n = n, m = m)
    suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
}

## Gene sets of increasing size, contiguous slices of the gene order so that
## building them consumes no RNG.  Only roast draws random numbers.
make_gene_sets <- function(genes, sizes = c(5, 20, 50)) {
    sets <- lapply(sizes, function(k) utils::head(genes, k))
    stats::setNames(sets, sprintf("set%02d", seq_along(sets)))
}

## The contrast is discovered from the design, never hard-coded:
## makeExampleDESeqDataSet() names its second column "condition_B_vs_A" while
## the airway example in the documentation names it "dex_trt_vs_untrt"
## (REBUILD_SPEC trap 8.4).  Both are spelled by DESeq2 / model.matrix().
last_contrast <- function(dd) colnames(dd$design)[ncol(dd$design)]

## tests/refdata/dsea_golden.rds, found whether the suite is run from the
## package root (R CMD check / devtools::test), from tests/testthat, or from
## tests/ (testthat::test_path resolves against the running test file).
golden_path <- function() {
    cands <- c(
        tryCatch(testthat::test_path("refdata", "dsea_golden.rds"), error = function(e) ""),
        "refdata/dsea_golden.rds",          # cwd = tests/testthat
        "../refdata/dsea_golden.rds",       # cwd = tests
        "tests/refdata/dsea_golden.rds",    # cwd = package root
        system.file("refdata", "dsea_golden.rds", package = "enrichit")
    )
    hit <- cands[nzchar(cands) & file.exists(cands)]
    if (length(hit) == 0L) {
        stop("dsea_golden.rds not found; regenerate it with ",
             "Rscript inst/scripts/make_dsea_golden.R", call. = FALSE)
    }
    hit[1]
}

## ---------------------------------------------------------------------------
## as_dsea(): the E matrix
## ---------------------------------------------------------------------------

test_that("as_dsea() E matches an independently written longhand reference", {
    dds <- make_fitted_dds()
    ct <- DESeq2::counts(dds)
    sf <- DESeq2::sizeFactors(dds)
    dd <- as_dsea(dds)

    ## Reference written one SAMPLE AT A TIME, with scalar arithmetic.  This is
    ## algebraically the same expression as as_dsea()'s matrix form but shares
    ## no code path with it: no matrix recycling, no byrow = TRUE, no
    ## vectorised log2.  Denominator = mean(colSums(ct)) * norm.factors.
    E_ref <- vapply(seq_len(ncol(ct)), function(j) {
        lib_j <- mean(colSums(ct)) * sf[j]
        log2((ct[, j] + 0.5) / (lib_j + 1) * 1e6)
    }, numeric(nrow(ct)))
    dim(E_ref) <- dim(ct)
    dimnames(E_ref) <- dimnames(ct)

    expect_equal(dd$E, E_ref, tolerance = 1e-12)
    ## The tolerance above is the contract; this is the stronger form, and it
    ## holds exactly (max|diff| == 0 on the 200-gene fixture and on the full
    ## 63677-gene airway dataset).
    expect_true(max(abs(dd$E - E_ref)) == 0)
    expect_true(identical(unname(dd$E), unname(E_ref)))
    expect_identical(dimnames(dd$E), dimnames(E_ref))
})

test_that("as_dsea() E does NOT use the double-counted colSums * norm.factors", {
    dds <- make_fitted_dds()
    ct <- DESeq2::counts(dds)
    sf <- DESeq2::sizeFactors(dds)
    dd <- as_dsea(dds)

    ## The pre-fix-3 denominator: norm.factors are depth-like (median-of-ratios,
    ## geometric mean ~ 1), so multiplying by the raw per-sample depth counts
    ## the depth twice and leaves an unmodellable -log2(colSums_j) offset in E.
    E_bad <- vapply(seq_len(ncol(ct)), function(j) {
        log2((ct[, j] + 0.5) / (colSums(ct)[j] * sf[j] + 1) * 1e6)
    }, numeric(nrow(ct)))
    dim(E_bad) <- dim(ct)
    dimnames(E_bad) <- dimnames(ct)

    ## The two formulas must NOT agree.  If they ever do, the denominator has
    ## silently reverted to the double-counted form.
    expect_false(isTRUE(all.equal(dd$E, E_bad, tolerance = 1e-6)))

    ## Their difference is exactly log2(colSums_j / mean(colSums)) -- constant
    ## down every column (a sample-level effect, not a gene-level one).  The
    ## match is approximate only because of the +1 in the denominator.
    offset <- colMeans(dd$E - E_bad)
    expect_equal(offset, log2(colSums(ct) / mean(colSums(ct))), tolerance = 1e-4)
    expect_true(all(abs(offset) > 1e-3))
    spread <- apply(dd$E - E_bad, 2, function(z) max(z) - min(z))
    expect_true(all(spread < 1e-6))
})

## ---------------------------------------------------------------------------
## Rank-deficient designs
## ---------------------------------------------------------------------------

test_that("a rank-deficient design errors with 'not of full rank' in all three methods", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20))

    ## Column 3 is exactly 2 * column 2, so rank([1, x, 2x]) = 2 < p = 3.
    ## R's qr() reports the deficiency through $rank instead of raising, so the
    ## engine has to check it explicitly (REBUILD_SPEC section 3.1).
    x <- as.numeric(dd$design[, 2])
    bad <- cbind(1, x, 2 * x)
    colnames(bad) <- c("(Intercept)", "b", "b2")
    expect_identical(qr(bad)$rank, 2L)

    ## Unweighted and weighted: both branches of lm_effects() must rank-check.
    ## nrot is passed only where it is a valid argument; camera and fry reject it
    ## outright, which is a separate contract tested further down.
    for (w in list(NULL, dd$weights)) {
        for (m in c("camera", "roast", "fry")) {
            dots <- if (m == "roast") list(nrot = 9) else list()
            expect_error(
                do.call(dsea, c(list(list(E = dd$E, weights = w, design = bad),
                                     gs, contrast = 3L, method = m), dots)),
                "not of full rank",
                info = sprintf("method = %s, weights = %s", m,
                               if (is.null(w)) "NULL" else "G x n matrix")
            )
        }
    }

    ## And it is not an artefact of which column is asked for.
    for (ci in 1:3) {
        expect_error(
            dsea(list(E = dd$E, weights = dd$weights, design = bad), gs,
                 contrast = ci, method = "camera"),
            "not of full rank",
            info = sprintf("contrast column %d", ci)
        )
    }
})

## ---------------------------------------------------------------------------
## Per-method return contracts
## ---------------------------------------------------------------------------

test_that("camera returns its documented contract", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20, 50))
    ## sort = FALSE so that the rows stay in the caller's order and the
    ## assertions below can be written against the input, not against a ranking.
    res <- dsea(dds, gs, contrast = last_contrast(dd), method = "camera",
                sort = FALSE)

    expect_s3_class(res, "data.frame")
    expect_true(all(c("set", "NGenes", "Direction", "PValue") %in% colnames(res)))
    expect_identical(nrow(res), length(gs))
    expect_identical(res$set, names(gs))
    expect_identical(as.integer(res$NGenes), c(5L, 20L, 50L))
    expect_true(all(res$Direction %in% c("Down", "Up")))
    expect_true(all(res$PValue >= 0 & res$PValue <= 1))
    expect_true(all(is.finite(res$PValue)))

    ## The p-value is the two-sided one, and the Direction must agree with the
    ## sign of the statistic that carries it.
    expect_true(all((res$Direction == "Up") == (res$Statistic > 0)))
})

test_that("fry returns its documented contract", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20, 50))
    res <- dsea(dds, gs, contrast = last_contrast(dd), method = "fry",
                sort = FALSE)

    expect_true(all(c("set", "NGenes", "Direction", "PValue") %in% colnames(res)))
    expect_identical(res$set, names(gs))
    expect_identical(as.integer(res$NGenes), c(5L, 20L, 50L))
    expect_true(all(res$Direction %in% c("Down", "Up")))
    expect_true(all(res$PValue >= 0 & res$PValue <= 1))

    ## fry's statistic is the moderated set t, whose sign is the direction.
    expect_true(all((res$Direction == "Up") == (res$Statistic > 0)))

    ## Two-sided and closed form: p = 2 * pt(-|t|, df.residual).  This pins the
    ## relationship between the returned statistic and the returned p-value, so
    ## a change in either one that breaks it is caught here rather than by the
    ## golden file alone.
    dfr <- ncol(dd$E) - ncol(dd$design)
    expect_equal(res$PValue,
                 2 * stats::pt(-abs(res$Statistic), df = dfr),
                 info = "fry PValue must equal 2*pt(-|t|, df.residual)")
})

test_that("roast returns its documented contract", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20, 50))
    res <- dsea(dds, gs, contrast = last_contrast(dd), method = "roast",
                nrot = 99, sort = FALSE)

    expect_true(all(c("set", "NGenes", "Direction", "PValue") %in% colnames(res)))
    expect_identical(res$set, names(gs))
    expect_identical(as.integer(res$NGenes), c(5L, 20L, 50L))
    expect_true(all(res$Direction %in% c("Down", "Up")))
    expect_true(all(res$PValue >= 0 & res$PValue <= 1))

    ## With the default set.statistic = "mean" and nrot rotations, the UpOrDown
    ## p-value is a rotation count on the grid {1, ..., nrot + 1}/(nrot + 1),
    ## so (nrot + 1) * p must be an integer.  (The engine does not apply the
    ## mid-p adjustment unless midp = TRUE is asked for.)
    nrot <- 99
    expect_true(all(abs(res$PValue * (nrot + 1) -
                            round(res$PValue * (nrot + 1))) < 1e-9))
})

test_that("roast direction follows the winning one-sided tail, as limma does", {
    ## REBUILD_ACCEPTANCE.md D4. limma:::mroast.default labels a set
    ##   Direction[Up] <- "Up"   where   Up <- pv[,"Up"] < pv[,"Down"]
    ## i.e. "which one-sided rotation tail won", NOT the sign of the active
    ## proportion difference. The two rules part company on a near-null set,
    ## where both active proportions are 0 so the proportion difference is 0
    ## and a sign-of-statistic rule degenerates. An earlier version of this
    ## test asserted "Down" for exactly that degenerate case; verified against
    ## limma with synchronised seeds, the reference answer is "Up" (the set's
    ## observed moderated t is positive), and we now match it.
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    one <- list(single = utils::head(rownames(dd$E), 1))
    res <- dsea(dds, one, contrast = last_contrast(dd), method = "roast", nrot = 99)
    expect_true(res$Direction %in% c("Down", "Up"))
    expect_true(all(res$PValue >= 0 & res$PValue <= 1))

    ## The direction must agree with the sign of the set's observed moderated
    ## t -- that is what decides which tail wins. Computed here without limma:
    eff <- st_effects_cpp(as.matrix(dd$E), as.matrix(dd$weights),
                          as.matrix(dd$design),
                          as.integer(ncol(dd$design)))
    dfr <- ncol(eff) - 1L
    s2 <- rowMeans(eff[, -1L, drop = FALSE]^2)
    sv <- st_squeeze_var_cpp(s2, as.double(dfr))
    obs_t <- eff[1L, 1L] / sqrt(sv$var.post[1L])
    expect_identical(res$Direction, if (obs_t > 0) "Up" else "Down")

    ## Statistic still reports the active-proportion difference, which for a
    ## near-null set is 0 while the direction is still meaningful.
    expect_true(is.finite(res$Statistic))
})

test_that("fry computes a real p-value for a single-gene set", {
    ## limma:::.fryEffects computes t.stat and PValue for EVERY set; a
    ## single-gene set's colMeans is just that gene's own effects row.  Only
    ## limma's PValue.Mixed column is gated on NGenes > 1, and this engine does
    ## not return that column.  (REBUILD_SPEC section 3.3 was corrected on
    ## 2026-10-14 after reading the limma source; an earlier version of this
    ## test asserted NA and was wrong.)
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    one <- list(single = utils::head(rownames(dd$E), 1))
    res <- dsea(dds, one, contrast = last_contrast(dd), method = "fry")
    expect_identical(as.integer(res$NGenes), 1L)
    expect_false(is.na(res$PValue))
    expect_true(res$PValue >= 0 && res$PValue <= 1)
    expect_true(res$Direction %in% c("Down", "Up"))
})

test_that("camera Direction matches the sign of the set's mean effect", {
    ## Regression test for the sign inversion (REBUILD_ACCEPTANCE.md D3).
    ## Needs no limma: spike one set up and one down, then check the labels.
    ## A set that is genuinely up-regulated must come back "Up".
    skip_if_not_installed("DESeq2")
    set.seed(7)
    dds <- DESeq2::makeExampleDESeqDataSet(n = 300, m = 8)
    ct <- DESeq2::counts(dds); rn <- rownames(ct)
    trt <- dds$condition == "B"
    ct[rn[1:30], trt] <- as.integer(round(ct[rn[1:30], trt] * 8))    # UP
    ct[rn[60:89], trt] <- as.integer(round(ct[rn[60:89], trt] / 4))  # DOWN
    dds2 <- DESeq2::DESeqDataSetFromMatrix(
        ct, as.data.frame(SummarizedExperiment::colData(dds)), ~condition)
    dds2 <- suppressMessages(DESeq2::DESeq(dds2, quiet = TRUE))
    res <- dsea(dds2, list(UPSET = rn[1:30], DOWNSET = rn[60:89]),
                contrast = "condition_B_vs_A", method = "camera")
    expect_identical(res$Direction[res$set == "UPSET"], "Up")
    expect_identical(res$Direction[res$set == "DOWNSET"], "Down")
    expect_true(res$Statistic[res$set == "UPSET"] > 0)
    expect_true(res$Statistic[res$set == "DOWNSET"] < 0)
})

## ---------------------------------------------------------------------------
## Method-argument coverage
## ---------------------------------------------------------------------------

test_that("roast set.statistic variants all run and reach different p-values", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20, 50))
    ct <- last_contrast(dd)

    pvs <- vapply(c("mean", "floormean", "mean50", "msq"), function(s) {
        set.seed(42)
        res <- dsea(dds, gs, contrast = ct, method = "roast",
                    set.statistic = s, nrot = 99)
        expect_true(all(res$PValue >= 0 & res$PValue <= 1))
        expect_true(all(res$Direction %in% c("Down", "Up")))
        paste(sprintf("%.12g", res$PValue), collapse = "|")
    }, character(1))

    ## They are four different statistics, so they must not all coincide.
    expect_true(length(unique(pvs)) > 1L)
    ## match.arg() reports the offending argument as 'arg', so match on its
    ## message rather than on the parameter name.
    expect_error(dsea(dds, gs, contrast = ct, method = "roast",
                      set.statistic = "nope", nrot = 9), "should be one of")
})

test_that("fry standardize variants all run and change the variance moderation", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20, 50))
    ct <- last_contrast(dd)

    pvs <- vapply(c("posterior.sd", "none", "residual.sd", "p2"), function(s) {
        set.seed(42)
        res <- dsea(dds, gs, contrast = ct, method = "fry", standardize = s)
        expect_true(all(res$PValue >= 0 & res$PValue <= 1))
        paste(sprintf("%.12g", res$PValue), collapse = "|")
    }, character(1))

    ## "posterior.sd" borrows strength across genes, the others do not; the four
    ## must not collapse onto one vector.
    expect_true(length(unique(pvs)) > 1L)
    expect_error(dsea(dds, gs, contrast = ct, method = "fry",
                      standardize = "nope"), "should be one of")
})

test_that("camera inter.gene.cor: a fixed value and NA (estimated) both work", {
    ## dsea()'s own documentation: "inter.gene.cor ... Default 0.01. Pass NA to
    ## estimate it from the data."  REBUILD_SPEC section 4 says the same:
    ## "inter.gene.cor default 0.01; NA -> estimated".
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20, 50))
    ct <- last_contrast(dd)

    fixed <- dsea(dds, gs, contrast = ct, method = "camera", inter.gene.cor = 0.01)
    expect_true(all(fixed$PValue >= 0 & fixed$PValue <= 1))

    ## NA must estimate the inter-gene correlation instead of using a constant.
    est <- dsea(dds, gs, contrast = ct, method = "camera", inter.gene.cor = NA)
    expect_true(all(est$PValue >= 0 & est$PValue <= 1))
    expect_identical(est$set, fixed$set)

    ## The two are different tests: with the correlation estimated per set the
    ## variance inflation factor changes, so at least one p-value must move.
    expect_false(isTRUE(all.equal(est$PValue, fixed$PValue)))
})

test_that("camera allow.neg.cor is validated", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20))
    ct <- last_contrast(dd)
    expect_error(dsea(dds, gs, contrast = ct, method = "camera",
                      allow.neg.cor = "yes"), "allow.neg.cor")
})

## ---------------------------------------------------------------------------
## Entry-point behaviour
## ---------------------------------------------------------------------------

test_that("dsea() returns the gene set names as row names", {
    ## This was a real defect: the result frame's row names were being cleared,
    ## which destroys the only key that ties a row back to its gene set and
    ## breaks every caller that matches arms by name (the benchmark parity
    ## scripts do exactly that).  limma's camera/fry/mroast keep them, so this
    ## package must too.
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20, 50, 100))
    ct <- last_contrast(dd)

    for (m in c("camera", "roast", "fry")) {
        set.seed(42)
        res <- dsea(dds, gs, contrast = ct, method = m, sort = FALSE)
        expect_identical(rownames(res), res$set,
                         info = sprintf("method = %s must keep the set names as row names", m))
        expect_identical(rownames(res), names(gs))
    }
})

test_that("engine accepts only 'native'", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20))
    ct <- last_contrast(dd)

    for (e in c("limma", "other", "")) {
        expect_error(dsea(dds, gs, contrast = ct, method = "camera", engine = e),
                     info = sprintf("engine = '%s'", e))
    }
    expect_error(dsea(dds, gs, contrast = ct, method = "camera", engine = 1),
                 "'engine' must be a single character string")

    ## The default is "native", and it must not be an error.
    expect_s3_class(dsea(dds, gs, contrast = ct, method = "camera"), "data.frame")
})

test_that("contrast can be a name, a column index, a one-non-zero vector, or a negative index", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20))
    p <- ncol(dd$design)

    ref <- dsea(dds, gs, contrast = colnames(dd$design)[p], method = "camera")
    expect_identical(ref$PValue,
                     dsea(dds, gs, contrast = p, method = "camera")$PValue)
    expect_identical(ref$PValue,
                     dsea(dds, gs, contrast = c(0, 1), method = "camera")$PValue)
    ## -p asks for the same column read the other way round.
    expect_identical(ref$PValue,
                     dsea(dds, gs, contrast = -p, method = "camera")$PValue)

    expect_error(dsea(dds, gs, contrast = "no_such_column", method = "camera"),
                 "not a column of the design matrix")
    expect_error(dsea(dds, gs, contrast = 99, method = "camera"), "out of range")
    expect_error(dsea(dds, gs, contrast = c(1, 1), method = "camera"),
                 "non-zero entry")
})

test_that("dsea() accepts a bare list(E, weights, design) and a DseaData object", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20))
    ct <- last_contrast(dd)

    from_dds <- dsea(dds, gs, contrast = ct, method = "camera")
    from_obj <- dsea(dd, gs, contrast = ct, method = "camera")
    from_list <- dsea(list(E = dd$E, weights = dd$weights, design = dd$design),
                      gs, contrast = ct, method = "camera")
    expect_identical(from_dds$PValue, from_obj$PValue)
    expect_identical(from_dds$PValue, from_list$PValue)

    ## weights = NULL is the unweighted path and must still work.
    unw <- dsea(list(E = dd$E, weights = NULL, design = dd$design),
                gs, contrast = ct, method = "camera")
    expect_true(all(unw$PValue >= 0 & unw$PValue <= 1))
    expect_false(isTRUE(all.equal(unw$PValue, from_dds$PValue)))
})

test_that("gene.weights must be a named numeric vector, or R NULL", {
    ## numeric(0) is the classic trap: it is numeric and has length 0, and it
    ## reaches the C++ boundary as an empty vector rather than as "absent".
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20))
    ct <- last_contrast(dd)
    ids <- rownames(dd$E)

    expect_error(dsea(dds, gs, contrast = ct, method = "fry",
                      gene.weights = numeric(0)), "gene.weights")
    expect_error(dsea(dds, gs, contrast = ct, method = "roast",
                      gene.weights = numeric(0), nrot = 9), "gene.weights")
    expect_error(dsea(dds, gs, contrast = ct, method = "fry",
                      gene.weights = runif(length(ids))), "gene.weights")
    expect_error(dsea(dds, gs, contrast = ct, method = "fry",
                      gene.weights = setNames(runif(10), ids[1:10])),
                 "does not cover all genes")

    ## A properly named vector is accepted and changes the result.
    gw <- stats::setNames(runif(length(ids)), ids)
    set.seed(42)
    with_gw <- dsea(dds, gs, contrast = ct, method = "roast",
                    gene.weights = gw, nrot = 99)
    set.seed(42)
    without <- dsea(dds, gs, contrast = ct, method = "roast", nrot = 99)
    expect_true(all(with_gw$PValue >= 0 & with_gw$PValue <= 1))
    expect_false(isTRUE(all.equal(with_gw$PValue, without$PValue)))

    ## A constant vector carries no information and is documented to be ignored.
    set.seed(42)
    const <- dsea(dds, gs, contrast = ct, method = "roast",
                  gene.weights = stats::setNames(rep(1, length(ids)), ids), nrot = 99)
    set.seed(42)
    expect_identical(const$PValue, without$PValue)
})

test_that("sort = TRUE/FALSE controls the row order", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20, 50, 100))
    ct <- last_contrast(dd)

    unsorted <- dsea(dds, gs, contrast = ct, method = "camera", sort = FALSE)
    expect_identical(unsorted$set, names(gs))
    expect_identical(rownames(unsorted), names(gs))

    sorted <- dsea(dds, gs, contrast = ct, method = "camera", sort = TRUE)
    ## dsea() orders by PValue, then -NGenes, then set name.
    o <- order(unsorted$PValue, -unsorted$NGenes, unsorted$set)
    expect_identical(sorted$set, unsorted$set[o])
    expect_identical(sorted$PValue, unsorted$PValue[o])
    expect_true(all(diff(sorted$PValue) >= 0))
    expect_setequal(sorted$set, names(gs))
    expect_identical(rownames(sorted), sorted$set)
})

test_that("roast midp is a documented no-op, not a silently dropped argument", {
    ## limma:::mroast.default applies midp ONLY to the multiple-testing columns
    ## (TwoSidedP2 and MixedP2, then tab$FDR <- pmax(...)). dsea() returns
    ## unadjusted p-values and no FDR, so midp cannot change any returned value.
    ## Pinned so that, if an FDR column is ever added, this test fails and the
    ## claim in ?dsea gets revisited instead of rotting.
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(20, 50))
    ct <- last_contrast(dd)
    set.seed(11)
    on <- suppressWarnings(dsea(dds, gs, contrast = ct, method = "roast",
                               nrot = 199, midp = TRUE))
    set.seed(11)
    off <- suppressWarnings(dsea(dds, gs, contrast = ct, method = "roast",
                                 nrot = 199, midp = FALSE))
    expect_identical(on$PValue, off$PValue)
    expect_false("FDR" %in% colnames(on))
    ## A nonsense midp must be rejected rather than swallowed.
    expect_error(dsea(dds, gs, contrast = ct, method = "roast", midp = "yes"),
                 "midp")
    expect_error(dsea(dds, gs, contrast = ct, method = "roast", midp = NA),
                 "midp")
})

test_that("unknown per-method arguments are rejected", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20))
    ct <- last_contrast(dd)

    expect_error(dsea(dds, gs, contrast = ct, method = "camera", nrot = 9),
                 "unused argument")
    expect_error(dsea(dds, gs, contrast = ct, method = "fry", inter.gene.cor = 0.5),
                 "unused argument")
    expect_error(dsea(dds, gs, contrast = ct, method = "roast", standardize = "none"),
                 "unused argument")
})

test_that("dsea() validates its inputs", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20))
    ct <- last_contrast(dd)

    expect_error(dsea(dds, list(Z = "NOT_A_GENE"), contrast = ct, method = "camera"),
                 "none of the gene set members")
    expect_error(dsea(list(E = dd$E, weights = dd$weights,
                           design = dd$design[1:3, ]), gs, contrast = ct,
                     method = "camera"), "must equal ncol")

    ## An E with NA (or non-finite) entries must be refused: the set tests
    ## require a complete matrix, and a silent NA would propagate into every
    ## p-value.
    E_na <- dd$E
    E_na[1, 1] <- NA_real_
    expect_error(dsea(list(E = E_na, weights = NULL, design = dd$design),
                      list(A = rownames(dd$E)[1:20]), contrast = ct,
                      method = "camera"), "NA")
    E_inf <- dd$E
    E_inf[1, 1] <- Inf
    expect_error(dsea(list(E = E_inf, weights = NULL, design = dd$design),
                      list(A = rownames(dd$E)[1:20]), contrast = ct,
                      method = "camera"), "non-finite")
})

test_that("p-values are returned unadjusted and Method labels the engine", {
    dds <- make_fitted_dds()
    dd <- as_dsea(dds)
    gs <- make_gene_sets(rownames(dd$E), c(5, 20, 50))
    ct <- last_contrast(dd)

    for (m in c("camera", "roast", "fry")) {
        set.seed(42)
        res <- dsea(dds, gs, contrast = ct, method = m)
        expect_identical(unique(res$Method), m)
        ## dsea() returns p-values only: no FDR / p.adjust / qvalue column is
        ## added, because the right multiple-testing family depends on how the
        ## caller collected the sets.
        expect_false(any(c("FDR", "p.adjust", "qvalue") %in% colnames(res)))
        expect_true(all(res$PValue >= 0 & res$PValue <= 1))
    }
})

## ---------------------------------------------------------------------------
## Golden reference
## ---------------------------------------------------------------------------

test_that("the golden reference p-values do not drift", {
    ## limma is gone from the package, so this .rds is the only guard against a
    ## silent engine change.  Regenerate it ONLY with
    ##     Rscript inst/scripts/make_dsea_golden.R
    ## and only when the engine is final -- see the header of that script.
    testthat::skip_if_not_installed("DESeq2")
    gpath <- golden_path()
    golden <- readRDS(gpath)

    ## Reproduce the exact conditions recorded in the file's provenance block.
    expect_identical(golden$provenance$seed, 42L)
    set.seed(golden$provenance$seed)
    dds <- DESeq2::makeExampleDESeqDataSet(n = golden$provenance$n_genes,
                                           m = golden$provenance$m_samples)
    dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
    dd <- as_dsea(dds)
    genes <- rownames(dd$E)
    contrast <- colnames(dd$design)[ncol(dd$design)]
    expect_identical(contrast, golden$provenance$contrast)

    gene_sets <- lapply(golden$n_genes_per_set, function(k) utils::head(genes, k))
    names(gene_sets) <- golden$set_names

    for (m in c("camera", "roast", "fry")) {
        set.seed(golden$provenance$seed)
        res <- dsea(dds, gene_sets, contrast = contrast, method = m, sort = FALSE)
        expect_identical(res$set, golden$set_names)
        expect_identical(as.integer(res$NGenes), as.integer(golden$n_genes_per_set))
        ## TOLERANCE, not bit-exactness. The golden was recorded on one CPU/BLAS;
        ## R CMD check runs on ubuntu/macos/windows with different floating-point
        ## summation order, so the last 1-2 ulps legitimately differ. tolerance = 0
        ## (testthat's exact comparison) failed CI on exactly this with diffs of
        ## ~2e-15 that are NOT regressions. 1e-10 is far above hardware noise yet
        ## still catches any real algorithmic change: the four bugs this golden
        ## guards against (QR basis sign, depth double-count, Direction inversion,
        ## fry single-gene NA) all moved p-values by >= 1e-4 or changed a label,
        ## i.e. >= 6 orders of magnitude beyond this bound.
        expect_equal(res$PValue, golden$pvalue[[m]], tolerance = 1e-10,
                     info = sprintf("%s p-values drifted from the golden reference", m))
        ## Direction is locked too: a p-value-only golden cannot catch a sign
        ## inversion, because TwoSided (and fry's 2*pt(-|t|)) are sign-blind.
        ## That is exactly how the camera Direction bug (REBUILD_ACCEPTANCE D3)
        ## survived the first round of checks.
        if (!is.null(golden$direction)) {
            expect_identical(res$Direction, golden$direction[[m]],
                             info = sprintf("%s Direction drifted", m))
        }
    }
})
