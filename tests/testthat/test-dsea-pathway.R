## Tests for the pathway-score route (PATHWAY_DSEA_SPEC.md).
##
## A10/A11 use limma as a REFERENCE ONLY.  The package itself never calls it; if
## these two golden tests fail, our base-R re-implementations have drifted and
## the numbers below are the ones to compare against.

skip_if_no <- function(p) testthat::skip_if_not_installed(p)

## ---- fixtures --------------------------------------------------------------

make_counts <- function(n = 12, G = 1200, seed = 1) {
    set.seed(seed)
    mu <- 2^runif(G, 3, 10)
    m <- matrix(rnbinom(G * n, mu = rep(mu, n), size = 10), G, n)
    dimnames(m) <- list(paste0("g", seq_len(G)), paste0("s", seq_len(n)))
    m
}

make_dds <- function(cnt, n_treated = NULL) {
    skip_if_no("DESeq2")
    suppressMessages(library(DESeq2))
    n <- ncol(cnt)
    if (is.null(n_treated)) n_treated <- n %/% 2
    cd <- data.frame(cond = factor(c(rep("trt", n_treated), rep("untrt", n - n_treated))),
                     row.names = colnames(cnt))
    suppressMessages(DESeq(DESeqDataSetFromMatrix(cnt, cd, ~cond), quiet = TRUE))
}

make_sets <- function(G, n_sets = 5, size = 100, seed = 2) {
    set.seed(seed)
    out <- lapply(seq_len(n_sets), function(i) paste0("g", sample.int(G, size)))
    names(out) <- paste0("SET", seq_len(n_sets))
    out
}

two_group_design <- function(n, n_treated = NULL) {
    if (is.null(n_treated)) n_treated <- n %/% 2
    g <- c(rep(1, n_treated), rep(0, n - n_treated))
    X <- cbind(intercept = 1, cond = g)
    rownames(X) <- paste0("s", seq_len(n))
    X
}

## ---- A4 / A9: design validation -------------------------------------------

test_that("a rank-deficient design is rejected before any fitting", {
    n <- 8
    ## dex nested in cell: the airway trap, in the numeric form model.matrix
    ## would actually produce
    X <- cbind(intercept = 1, cell_b = rep(c(0, 1), each = 4), dex = rep(c(1, 0), each = 4))
    expect_equal(qr(X)$rank, 2L)
    expect_error(enrichit:::.dsea_check_design(X, paste0("s", seq_len(n))), "not full rank")
})

test_that("a non-numeric design is rejected with usable advice", {
    X <- data.frame(cond = rep(c("a", "b"), 4))
    expect_error(enrichit:::.dsea_check_design(X, paste0("s", 1:8)), "model.matrix")
})

test_that("design row names must match the sample names", {
    X <- two_group_design(8)
    rownames(X) <- paste0("other", 1:8)
    expect_error(enrichit:::.dsea_check_design(X, paste0("s", 1:8)), "do not match")
})

test_that("a design of the wrong length is rejected", {
    expect_error(enrichit:::.dsea_check_design(two_group_design(8), paste0("s", 1:7)),
                 "rows but")
})

## ---- A7: contrast resolution ----------------------------------------------

test_that("contrast resolves by name and by index to the same column", {
    X <- two_group_design(8)
    samples <- rownames(X)
    expect_equal(enrichit:::.dsea_resolve_one("cond", X),
                 enrichit:::.dsea_resolve_one(2L, X))
    expect_equal(enrichit:::.dsea_resolve_one("cond", X), 2L)
    expect_error(enrichit:::.dsea_resolve_one("nope", X), "not found")
    expect_error(enrichit:::.dsea_resolve_one(99L, X), "1..2")
    ## NULL tests every non-intercept column
    expect_equal(enrichit:::.dsea_resolve_contrasts(NULL, X), 2L)
})

## ---- input validation ------------------------------------------------------

test_that("counts input is validated", {
    x <- make_counts(n = 6, G = 50)
    expect_error(enrichit:::.dsea_counts(x * 2.5), "log-scale expression")
    x2 <- x; x2[1, 1] <- NA
    expect_error(enrichit:::.dsea_counts(x2), "missing counts")
    x3 <- x; x3[1, 1] <- -1
    expect_error(enrichit:::.dsea_counts(x3), "negative counts")
    x4 <- x; rownames(x4) <- NULL
    expect_error(enrichit:::.dsea_counts(x4), "row names")
    expect_error(enrichit:::.dsea_counts(42), "must be a DESeqDataSet")
})

test_that("gene sets are expanded, and unusable ones are rejected", {
    gs <- enrichit:::.dsea_expand_sets(list(A = c("g1", "g2"), B = list("g3", c("g4", "g5"))))
    expect_identical(gs$A, c("g1", "g2"))
    expect_identical(gs$B, c("g3", "g4", "g5"))
    expect_error(enrichit:::.dsea_expand_sets(list(A = character(0))), "at least one")
    expect_error(enrichit:::.dsea_expand_sets(list(c("g1", "g2"))), "named list")
    ## ssGSEA weights a duplicated member twice, so this is an error, not a
    ## silent de-duplication: the caller almost certainly did not mean it.
    expect_error(enrichit:::.dsea_expand_sets(list(A = c("g1", "g2", "g1"))),
                 "duplicate identifiers")
})

## ---- A1 / A3 / A5 / A6: scoring -------------------------------------------

test_that("dsea_scores returns a well-formed, reproducible score matrix", {
    skip_if_no("GSVA")
    cnt <- make_counts()
    dds <- make_dds(cnt)
    sets <- make_sets(nrow(cnt))
    s <- dsea_scores(dds, sets)

    ## A1: shape and names
    expect_identical(dim(s$scores), c(length(sets), ncol(cnt)))
    expect_identical(rownames(s$scores), names(sets))
    expect_identical(colnames(s$scores), colnames(cnt))
    expect_s3_class(s, "DseaScores")

    ## A3: deterministic
    expect_identical(dsea_scores(dds, sets)$scores, s$scores)

    ## A5: provenance is complete and serialisable
    expect_true(all(c("transform", "transform_args", "score_method", "n_input_sets",
                      "n_tested_sets", "dropped_fraction", "vst_version",
                      "gsva_version") %in% names(s$meta)))
    expect_false(s$meta$transform_args$blind)   # blind = FALSE is deliberate
    skip_if_no("jsonlite")
    expect_silent(jsonlite::toJSON(s$meta))

    ## A6: dropped sets are counted, not silently lost
    tiny <- make_sets(nrow(cnt), n_sets = 3, size = 3)
    s2 <- suppressMessages(dsea_scores(dds, tiny))
    expect_equal(s2$meta$n_tested_sets, 0L)
})

test_that("sets below min_size are dropped and reported", {
    skip_if_no("GSVA")
    cnt <- make_counts()
    dds <- make_dds(cnt)
    sets <- c(make_sets(nrow(cnt), n_sets = 3, size = 120),
              list(SMALL = paste0("g", 1:5)))
    s <- suppressMessages(dsea_scores(dds, sets, min_size = 10L))
    expect_equal(s$meta$n_input_sets, 4L)
    expect_equal(s$meta$n_tested_sets, 3L)
    expect_gt(s$meta$dropped_fraction, 0)
    expect_false("SMALL" %in% rownames(s$scores))
})

## ---- A10: golden tests against limma --------------------------------------

test_that("A10: fitFDist re-implementation matches limma exactly", {
    skip_if_no("limma")
    set.seed(11)
    for (dof in c(3, 6, 25)) {
        m <- 400
        s2 <- rchisq(m, df = dof) / dof * exp(rnorm(m, 0, 0.5))
        mine <- enrichit:::.dsea_fit_prior(s2, dof)
        theirs <- limma:::fitFDist(s2, dof)
        expect_equal(mine$df2, theirs$df2, tolerance = 1e-12)
        expect_equal(mine$scale, theirs$scale, tolerance = 1e-12)
    }
})

test_that("A10: OLS coefficients match limma::lmFit exactly", {
    skip_if_no("limma")
    set.seed(12)
    n <- 10; m <- 50
    X <- two_group_design(n)
    Z <- matrix(rnorm(m * n), m, n)
    rownames(Z) <- paste0("p", seq_len(m))
    mine <- enrichit:::.dsea_fit(Z, X, 2L, "none")
    theirs <- limma::lmFit(Z, X)
    expect_equal(unname(mine$logFC), unname(theirs$coefficients[, 2]), tolerance = 1e-12)
    ## and the variance multiplier is [(X'X)^-1]_ci,ci, not [R^-1]_ci,ci
    v <- solve(crossprod(X))[2, 2]
    expect_equal(unname(mine$se^2 / mine$moderation$se2), rep(v, m), tolerance = 1e-12)
})

test_that("A2: a standardised t is O(1) under the null", {
    set.seed(13)
    n <- 20; m <- 500
    X <- two_group_design(n)
    Z <- matrix(rnorm(m * n), m, n)          # no signal anywhere
    rownames(Z) <- paste0("p", seq_len(m))
    fit <- enrichit:::.dsea_fit(Z, X, 2L, "eb")
    expect_gt(median(abs(fit$t)), 0.5)
    expect_lt(median(abs(fit$t)), 2)
    ## standardised noise gives logFC ~ sqrt(v), v = (X'X)^-1[2,2] = 0.2 here
    expect_lt(median(abs(fit$logFC)), 0.6)
})

## ---- A11: consensus correlation cross-check -------------------------------

test_that("A11: consensus correlation tracks limma::duplicateCorrelation", {
    skip_if_no("limma")
    set.seed(14)
    n <- 12; m <- 120
    X <- two_group_design(n)
    shared <- matrix(rnorm(m * n), m, n)      # one strong common factor
    Z <- shared + matrix(rnorm(m * n) * 1.5, m, n)
    rho_ours <- enrichit:::.dsea_consensus_correlation(Z, X, NULL)
    rho_limma <- limma::duplicateCorrelation(Z, X)$consensus.correlation
    expect_true(is.finite(rho_ours))
    expect_lte(rho_ours, 0.99)
    expect_lt(abs(rho_ours - rho_limma), 0.05)
})

test_that("consensus correlation is zero when pathways are independent", {
    set.seed(15)
    n <- 12; m <- 300
    X <- two_group_design(n)
    Z <- matrix(rnorm(m * n), m, n)
    expect_lt(enrichit:::.dsea_consensus_correlation(Z, X, NULL), 0.3)
})

## ---- dsea() end to end ----------------------------------------------------

test_that("dsea returns contrasts, scores and provenance", {
    skip_if_no("GSVA")
    n <- 12; G <- 1500
    set.seed(21)
    cnt <- matrix(rnbinom(G * n, mu = rep(2^runif(G, 6.5, 9), n), size = 20), G, n,
                  dimnames = list(paste0("g", seq_len(G)), paste0("s", seq_len(n))))
    de <- paste0("g", 1:200)
    mult <- matrix(1, G, n); mult[seq_along(de), 1:6] <- 4   # treated = columns 1:6
    cnt <- round(cnt * mult)
    dds <- make_dds(cnt)
    ## A realistic set COLLECTION, not three sets.  ssGSEA is a within-sample rank
    ## statistic, so with only a handful of disjoint sets, enrichment in one
    ## mechanically depresses the others and every pathway ends up about equally
    ## "significant".  That is the composition-bias failure mode recorded in
    ## PATHWAY_DSEA_SPEC.md section 8, and it is why the benchmark must score
    ## thousands of overlapping sets rather than a toy list.
    sets <- c(list(SIG = de), lapply(seq_len(60), function(i) paste0("g", sample(600:G, 120))))
    names(sets)[-1] <- paste0("BG", seq_len(60))

    res <- dsea(dds, sets, design = two_group_design(ncol(cnt)), contrast = "cond")

    expect_s3_class(res, "DseaResult")
    expect_length(res$contrasts, 1L)
    r <- res$contrasts$cond
    expect_true(all(c("logFC","se","t","p.value","adj.p.value","direction") %in% names(r)))
    expect_identical(names(r$logFC), rownames(res$scores))
    expect_true(all(r$p.value >= 0 & r$p.value <= 1))
    expect_identical(levels(r$direction), c("up","down","ns"))
    expect_true(res$moderation$method == "eb")
    expect_true(res$correlation$consensus >= 0)

    ## the planted set must rank first and move in the planted direction
    expect_equal(unname(which.min(r$p.value)), which(names(r$p.value) == "SIG"))
    expect_gt(r$logFC[["SIG"]], 0)
    expect_identical(as.character(r$direction[["SIG"]]), "up")
    ## A2 in practice: nothing exploded in magnitude
    expect_lt(median(abs(r$t)), 1e3)
})

test_that("dsea reports every non-intercept column when contrast is NULL", {
    skip_if_no("GSVA")
    cnt <- make_counts(n = 12, G = 1500, seed = 22)
    dds <- make_dds(cnt)
    sets <- list(S1 = paste0("g", 1:120), S2 = paste0("g", 200:320))
    X <- two_group_design(ncol(cnt))
    X <- cbind(X, batch = rep(c(0, 1), 6))
    res <- dsea(dds, sets, design = X, contrast = NULL)
    expect_identical(names(res$contrasts), c("cond", "batch"))
})

test_that("dsea standardisation removes between-pathway scale differences", {
    Z <- rbind(pathwayA = c(1, 2, 3, 4, 5, 6, 7, 8),
               pathwayB = c(100, 200, 300, 400, 500, 600, 700, 800))
    Zs <- enrichit:::.dsea_standardize(Z)
    expect_equal(unname(Zs[1, ]), unname(Zs[2, ]), tolerance = 1e-12)
    expect_identical(attr(Zs, "dropped_constant"), character(0))
    ## a constant row is dropped and named, not silently divided by zero
    Zc <- rbind(Z, flat = rep(2, 8))
    expect_equal(attr(enrichit:::.dsea_standardize(Zc), "dropped_constant"), "flat")
})

test_that("duplicate members are rejected rather than silently de-duplicated", {
    skip_if_no("GSVA")
    cnt <- make_counts(n = 8, G = 800, seed = 23)
    dds <- make_dds(cnt)
    expect_error(dsea_scores(dds, list(A = c(paste0("g", 1:50), paste0("g", 1:20)))),
                 class = "error")
})