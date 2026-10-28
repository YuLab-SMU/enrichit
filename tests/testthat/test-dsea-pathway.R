## Tests for design-aware GSEA and the ssGSEA scorer (PATHWAY_DSEA_SPEC.md).

skip_if_no <- function(p) testthat::skip_if_not_installed(p)

make_counts <- function(n = 8, G = 1200, seed = 1) {
    set.seed(seed)
    mu <- 2^runif(G, 3, 10)
    m <- matrix(rnbinom(G * n, mu = rep(mu, n), size = 20), G, n)
    dimnames(m) <- list(paste0("g", seq_len(G)), paste0("s", seq_len(n)))
    m
}

two_group <- function(n = 8, n_trt = 4) {
    g <- c(rep(1, n_trt), rep(0, n - n_trt))
    X <- cbind(intercept = 1, cond = g)
    rownames(X) <- paste0("s", seq_len(n))
    X
}

make_dds <- function(cnt) {
    skip_if_no("DESeq2")
    n <- ncol(cnt)
    cd <- data.frame(cond = factor(c(rep("trt", n %/% 2), rep("untrt", n - n %/% 2))),
                     row.names = colnames(cnt))
    suppressMessages(DESeq2::DESeq(
        DESeq2::DESeqDataSetFromMatrix(cnt, cd, ~cond), quiet = TRUE))
}

## ---- ssGSEA ----------------------------------------------------------------

test_that("ssgsea_scores returns a well formed, reproducible matrix", {
    E <- make_counts(n = 6, G = 600)
    sets <- list(A = paste0("g", 1:80), B = paste0("g", 200:300), C = paste0("g", 400:500))
    S <- ssgsea_scores(E, sets)
    expect_identical(dim(S), c(3L, 6L))
    expect_identical(rownames(S), names(sets))
    expect_identical(colnames(S), colnames(E))
    expect_true(all(is.finite(S)))
    expect_identical(ssgsea_scores(E, sets), S)   # deterministic
})

test_that("a set with no member in a sample scores zero, not NA", {
    E <- make_counts(n = 4, G = 400)
    sets <- list(half = paste0("g", 1:100))
    S <- ssgsea_scores(E, sets, min_size = 10, max_size = 500)
    expect_false(anyNA(S))
    ## a set whose members are absent from the matrix is filtered out, not scored
    expect_error(ssgsea_scores(E, list(none = c("zzz1", "zzz2")), min_size = 1),
                 "no gene set has between")
})

test_that("ssGSEA responds to the shift the set actually carries", {
    set.seed(5); n <- 8; G <- 500
    E <- matrix(rnorm(G * n), G, n, dimnames = list(paste0("g", 1:G), paste0("s", 1:n)))
    sig <- paste0("g", 1:60)
    E[sig, 1:4] <- E[sig, 1:4] + 3
    S <- ssgsea_scores(E, list(SIG = sig, BG = paste0("g", 300:400)))
    ## the shifted set must score higher in the treated samples than the control
    expect_gt(mean(S["SIG", 1:4]), mean(S["SIG", 5:8]))
    expect_gt(mean(S["SIG", 1:4] - S["SIG", 5:8]), mean(S["BG", 1:4] - S["BG", 5:8]))
})

test_that("gene set input is validated", {
    E <- make_counts(n = 4, G = 300)
    expect_error(ssgsea_scores(E, list(A = c("g1", "g1", "g2"))), "duplicate identifiers")
    expect_error(ssgsea_scores(E, list(A = character(0))), "at least one")
    expect_error(ssgsea_scores(E, list(c("g1", "g2"))), "named list")
    expect_error(ssgsea_scores(unname(E), list(A = paste0("g", 1:20))), "row names")
})

## ---- design validation -----------------------------------------------------

test_that("a rank-deficient design is rejected before any fitting", {
    X <- cbind(intercept = 1, cell_b = rep(c(0, 1), each = 4), dex = rep(c(1, 0), each = 4))
    expect_equal(qr(X)$rank, 2L)
    expect_error(enrichit:::.dsea_check_design(X, paste0("s", 1:8)), "not full rank")
})

test_that("a non-numeric or misaligned design is rejected with usable advice", {
    expect_error(enrichit:::.dsea_check_design(data.frame(cond = rep(c("a", "b"), 4)),
                                              paste0("s", 1:8)), "model.matrix")
    X <- two_group(8)
    rownames(X) <- paste0("other", 1:8)
    expect_error(enrichit:::.dsea_check_design(X, paste0("s", 1:8)), "do not match")
    expect_error(enrichit:::.dsea_check_design(two_group(8), paste0("s", 1:7)), "rows but")
})

test_that("contrast resolves by name and by index to the same column", {
    X <- two_group(8)
    expect_equal(enrichit:::.dsea_resolve_one("cond", X), enrichit:::.dsea_resolve_one(2L, X))
    expect_error(enrichit:::.dsea_resolve_one("nope", X), "not found")
    expect_error(enrichit:::.dsea_resolve_one(9L, X), "1..2")
})

## ---- the ranking statistic -------------------------------------------------

test_that("the two-group ranking is exactly the group mean difference", {
    set.seed(3); n <- 8; G <- 200
    E <- matrix(rnorm(G * n), G, n, dimnames = list(paste0("g", 1:G), paste0("s", 1:n)))
    X <- two_group(n)
    r <- enrichit:::.dsea_rank_by_contrast(E, X, 2L)
    expect_equal(unname(r), unname(rowMeans(E[, 1:4]) - rowMeans(E[, 5:8])))
})

test_that("a nuisance term is removed before ranking, not ignored", {
    ## The design must be UNBALANCED for this test to mean anything: with a
    ## fully crossed 2x2 the contrast is orthogonal to the nuisance term and the
    ## fitted coefficient equals the raw difference exactly, so a broken
    ## adjustment would pass.  Here batch is unbalanced across condition.
    set.seed(4); n <- 8; G <- 300
    E <- matrix(rnorm(G * n), G, n, dimnames = list(paste0("g", 1:G), paste0("s", 1:n)))
    E[, 6:8] <- E[, 6:8] + 2                       # a batch offset
    cond  <- c(1, 1, 1, 0, 0, 0, 0, 0)
    batch <- c(0, 0, 0, 0, 0, 1, 1, 1)
    X <- cbind(intercept = 1, cond = cond, batch = batch)
    expect_equal(qr(X)$rank, 3L)                   # still full rank
    fitted <- enrichit:::.dsea_rank_by_contrast(E, X, 2L)
    raw <- rowMeans(E[, cond == 1, drop = FALSE]) - rowMeans(E[, cond == 0, drop = FALSE])
    ## unbalanced, so the adjustment must move the estimate
    expect_gt(max(abs(fitted - raw)), 0.5)
    ## and it must be the OLS coefficient
    ref <- solve(crossprod(X), crossprod(X, t(E)))[2, ]
    expect_equal(unname(fitted), unname(ref), tolerance = 1e-10)
})

## ---- dsea() end to end -----------------------------------------------------

test_that("dsea finds a planted set and returns a GSEA table", {
    skip_if_no("DESeq2")
    set.seed(11); n <- 8; G <- 1500
    cnt <- matrix(rnbinom(G * n, mu = rep(2^runif(G, 6, 9), n), size = 20), G, n,
                  dimnames = list(paste0("g", 1:G), paste0("s", 1:n)))
    sig <- paste0("g", 1:120)
    mult <- matrix(1, G, n); mult[seq_along(sig), 1:4] <- 4
    cnt <- round(cnt * mult)
    dds <- make_dds(cnt)
    ## a realistic collection, not three sets
    sets <- c(list(SIG = sig),
              lapply(seq_len(40), function(i) paste0("g", sample(400:G, 80))))
    names(sets)[-1] <- paste0("BG", seq_len(40))

    ## the planted set is strong enough to reach the eps floor, which gsea()
    ## announces; the assertion is about the returned table
    res <- suppressWarnings(suppressMessages(
        dsea(dds, sets, design = two_group(n), contrast = "cond", seed = 7, verbose = FALSE)))
    expect_s3_class(res, "data.frame")
    ## gsea() returns ID/NES/pvalue/setSize/... and no multiplicity column;
    ## dsea() adds `padj` so that the front end hands back something directly
    ## usable, and records which contrast was tested.
    expect_true(all(c("ID", "NES", "pvalue", "padj", "contrast") %in% names(res)))
    expect_identical(unique(res$contrast), "cond")
    ## the planted set must rank first by p-value, and survive correction
    expect_equal(as.character(res$ID[which.min(res$pvalue)]), "SIG")
    expect_lt(res$padj[res$ID == "SIG"], 0.05)
    ## and its NES must have the sign of the planted direction
    expect_gt(res$NES[res$ID == "SIG"], 0)
})

test_that("dsea rejects count-like inputs that are not counts", {
    bad <- matrix(runif(100), 10, 10, dimnames = list(paste0("g", 1:10), paste0("s", 1:10)))
    expect_error(enrichit:::.dsea_counts(bad), "log-scale expression")
})

## ---- dsea_scores() ---------------------------------------------------------

test_that("dsea_scores records its provenance and drops sets it cannot score", {
    skip_if_no("DESeq2")
    cnt <- make_counts(n = 8, G = 1500, seed = 21)
    dds <- make_dds(cnt)
    sets <- c(list(A = paste0("g", 1:100), B = paste0("g", 200:320)),
              list(TINY = paste0("g", 1:3)))
    sc <- suppressMessages(dsea_scores(dds, sets, min_size = 10L))
    expect_s3_class(sc, "DseaScores")
    expect_identical(dim(sc$scores), c(2L, 8L))
    expect_equal(sc$meta$n_input_sets, 3L)
    expect_equal(sc$meta$n_scored_sets, 2L)
    expect_gt(sc$meta$dropped_fraction, 0)
    expect_false(sc$meta$transform_args$blind)     # blind = FALSE is deliberate
    expect_true(all(c("transform", "transform_args", "scorer", "exponent",
                      "set_sizes") %in% names(sc$meta)))
    skip_if_no("jsonlite")
    expect_silent(jsonlite::toJSON(sc$meta))
})

test_that("bare counts matrices work on both paths", {
    skip_if_no("DESeq2")
    cnt <- make_counts(n = 8, G = 1200, seed = 31)
    sc <- dsea_scores(cnt, list(A = paste0("g", 1:80), B = paste0("g", 200:300)),
                      transform = "logcpm")
    expect_identical(dim(sc$scores), c(2L, 8L))
    expect_true(all(is.finite(sc$scores)))
})

## ---- contrast selection is a correctness question --------------------------
##
## The default used to be `contrast = 1L`, which for the standard
## `model.matrix(~ condition)` design is the INTERCEPT column.  Ranking by the
## intercept ranks by overall mean expression, so GSEA reports highly expressed
## pathways instead of differentially enriched ones, and nothing errors.  These
## tests pin the fixed behaviour.

test_that("the default contrast is the condition column, never the intercept", {
    skip_if_no("DESeq2")
    cnt <- make_counts(n = 8, G = 1200, seed = 41)
    dds <- make_dds(cnt)
    X <- two_group(8)                      # (Intercept), cond
    sets <- list(A = paste0("g", 1:80), B = paste0("g", 200:300))
    res <- dsea(dds, sets, design = X, seed = 1, verbose = FALSE)
    expect_identical(unique(res$contrast), "cond")
    expect_false(identical(unique(res$contrast), "(Intercept)"))
})

test_that("an intercept column cannot be used as the contrast", {
    cnt <- make_counts(n = 8, G = 1200, seed = 42)
    X <- two_group(8)
    expect_error(enrichit:::.dsea_check_contrast(X, 1L), "intercept")
    expect_error(enrichit:::.dsea_check_contrast(X, "(Intercept)"), "intercept")
})

test_that("a design with more than two columns needs an explicit contrast", {
    skip_if_no("DESeq2")
    cnt <- make_counts(n = 8, G = 1200, seed = 43)
    dds <- make_dds(cnt)
    g <- c(rep(1, 4), rep(0, 4))
    X <- cbind(intercept = 1, cond = g, batch = c(0, 0, 1, 1, 0, 0, 1, 1))
    expect_equal(qr(X)$rank, 3L)
    sets <- list(A = paste0("g", 1:80))
    expect_error(dsea(dds, sets, design = X, seed = 1, verbose = FALSE), "contrast")
})

test_that("a no-intercept design is rejected: its columns are cell means", {
    cnt <- make_counts(n = 8, G = 1200, seed = 44)
    ## ~0 + cond: two indicator columns, no intercept.  Each coefficient is a
    ## group mean, not a difference, so ranking by one is not a differential
    ## ranking at all.
    X0 <- cbind(condB = c(1, 1, 1, 1, 0, 0, 0, 0),
                condA = c(0, 0, 0, 0, 1, 1, 1, 1))
    rownames(X0) <- paste0("s", 1:8)
    expect_error(enrichit:::.dsea_check_contrast(X0, 2L), "cell mean|indicator|group mean")
})

test_that("an untreated covariate column with a slope is still allowed", {
    ## A continuous column is a real coefficient, not a group mean: it must not
    ## be rejected by the indicator guard.
    X <- cbind(intercept = 1, dose = seq(0, 3.5, length.out = 8))
    expect_silent(enrichit:::.dsea_check_contrast(X, 2L))
})

test_that("dsea warns when the fitted dds design differs from the passed design", {
    skip_if_no("DESeq2")
    cnt <- make_counts(n = 8, G = 1200, seed = 45)
    dds <- make_dds(cnt)                       # fitted with ~ cond
    cd <- data.frame(cond = dds$cond,
                     batch = factor(c(1, 2, 1, 2, 1, 2, 1, 2)),
                     row.names = colnames(cnt))
    X <- model.matrix(~ batch + cond, cd)      # richer than design(dds)
    sets <- list(A = paste0("g", 1:80))
    ci <- grep("cond", colnames(X))
    expect_warning(dsea(dds, sets, design = X, contrast = ci, seed = 1, verbose = FALSE),
                   "design")
})

## ---- gene set validation is the same on both entry points -----------------

test_that("dsea rejects the gene sets dsea_scores rejects", {
    skip_if_no("DESeq2")
    cnt <- make_counts(n = 8, G = 1200, seed = 46)
    dds <- make_dds(cnt)
    X <- two_group(8)
    expect_error(dsea(dds, list(A = c("g1", "g1", "g2")), design = X,
                      contrast = 2L, seed = 1, verbose = FALSE),
                 "duplicate identifiers")
    expect_error(dsea(dds, list(A = character(0)), design = X,
                      contrast = 2L, seed = 1, verbose = FALSE),
                 "at least one")
})

## ---- meta reports the size actually scored ---------------------------------

test_that("set_sizes records the overlap with the data, not the raw set size", {
    skip_if_no("DESeq2")
    cnt <- make_counts(n = 8, G = 1200, seed = 47)
    dds <- make_dds(cnt)
    ## 40 identifiers, only 20 of which exist in the counts matrix
    sets <- list(PARTLY = c(paste0("g", 1:20), paste0("absent", 1:20)))
    sc <- suppressMessages(dsea_scores(dds, sets, min_size = 10L))
    expect_equal(unname(sc$meta$set_sizes[["PARTLY"]]), 20L)
    expect_equal(unname(sc$meta$raw_set_sizes[["PARTLY"]]), 40L)
})
