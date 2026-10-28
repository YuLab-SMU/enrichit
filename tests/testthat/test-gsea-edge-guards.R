library(testthat)
library(enrichit)

## Regression tests for the four defects found in the 2026-10-27 review of the
## multilevel engine and its R layer:
##   D2  the "too few same-side permutations" guard was a count, so nPermSimple <= 10
##       turned every pathway NA
##   D1  eps was not clamped, so eps > 1 returned p = eps > 1
##   D4  an all-zero ranking produced all-NA ranks and a full table of NaN
##   D3  a gene set covering the entire ranked list returned NaN, where fgsea
##       excludes it

make_input <- function(n = 400, top_shift = 4) {
    ids <- sprintf("G%04d", seq_len(n))
    set.seed(42)
    st <- rnorm(n)
    names(st) <- ids
    st[1:40] <- st[1:40] + top_shift
    st[(n - 39):n] <- st[(n - 39):n] - top_shift
    list(gl = sort(st, decreasing = TRUE), ids = ids)
}

## gson is a Suggests; every test below builds a GSON, so guard centrally
skip_without_gson <- function() skip_if_not_installed("gson")

make_gson <- function(sets) {
    gson::gson(
        gsid2gene = data.frame(gsid = rep(names(sets), lengths(sets)),
                               gene = unlist(sets, use.names = FALSE)),
        gsid2name = data.frame(gsid = names(sets), name = names(sets))
    )
}

test_that("D2: a small nPermSimple no longer turns every pathway NA", {
    skip_without_gson()
    x <- make_input()
    g <- make_gson(list(top = x$ids[1:40], bot = x$ids[361:400]))

    res <- suppressWarnings(gsea_gson(x$gl, g, seed = 1, verbose = FALSE,
                                      pvalueCutoff = 1, reportNA = TRUE,
                                      nPermSimple = 10))
    expect_s4_class(res, "gseaResult")
    ## before the fix this returned NA for every row
    expect_lt(sum(is.na(res@result$pvalue)), nrow(res@result))
    expect_true(all(res@result$pvalue <= 1, na.rm = TRUE))
})

test_that("D2: the guard is a rate, so defaults are unchanged", {
    skip_without_gson()
    ## 10 per 1000 == 1%; the default nPermSimple = 1000 must therefore behave
    ## exactly as the old count-based threshold did. Both settings below sit on
    ## the same side of the boundary as before.
    x <- make_input()
    g <- make_gson(list(top = x$ids[1:40], bot = x$ids[361:400]))
    res <- suppressWarnings(gsea_gson(x$gl, g, seed = 7, verbose = FALSE,
                                      pvalueCutoff = 1, reportNA = TRUE))
    expect_identical(nrow(res@result), 2L)
    expect_false(anyNA(res@result$pvalue))
})

test_that("D2: dropping NA rows says how many went", {
    skip_without_gson()
    ## the message is reachable through gsea_gson() only for a degenerate null,
    ## which is hard to construct, so the helper is pinned directly
    df <- data.frame(ID = c("A", "B", "C"), pvalue = c(0.01, NA, 0.02),
                     stringsAsFactors = FALSE)

    expect_warning(
        kept <- .gsea_drop_na_rows(df, reportNA = FALSE),
        "1 of 3 gene set"
    )
    expect_identical(kept$ID, c("A", "C"))

    ## reportNA = TRUE keeps them, silently
    expect_silent(kept2 <- .gsea_drop_na_rows(df, reportNA = TRUE))
    expect_identical(nrow(kept2), 3L)

    ## nothing to drop, nothing to say
    expect_silent(.gsea_drop_na_rows(df[c(1, 3), ], reportNA = FALSE))
    ## empty and NULL inputs are passthroughs
    expect_silent(expect_identical(nrow(.gsea_drop_na_rows(df[0, ], FALSE)), 0L))
    expect_null(.gsea_drop_na_rows(NULL, reportNA = FALSE))
})

test_that("D1: eps outside [0, 1] is rejected instead of returning p > 1", {
    skip_without_gson()
    x <- make_input()
    g <- make_gson(list(top = x$ids[1:40]))

    ## list(), not c(): a mixed vector is coerced to character, so `eps = 2`,
    ## -0.1 and c(0.1, 0.2) would all arrive as strings and every iteration would
    ## exercise the same "not numeric" branch. That version of this test passed
    ## while the out-of-range case -- the one that produced p = eps > 1 -- went
    ## completely untested.
    for (bad in list(2, 1.5, -0.1, NA_real_, "0.1", c(0.1, 0.2))) {
        expect_error(
            suppressWarnings(gsea_gson(x$gl, g, verbose = FALSE, eps = bad)),
            "`eps` must be a single number in \\[0, 1\\]"
        )
    }

    ## the boundary values themselves are legal
    for (good in c(0, 1e-10, 0.5, 1)) {
        expect_s4_class(
            suppressWarnings(gsea_gson(x$gl, g, seed = 1, verbose = FALSE,
                                       pvalueCutoff = 1, eps = good)),
            "gseaResult"
        )
    }
})

test_that("D1: reported p-values stay within [0, 1]", {
    skip_without_gson()
    x <- make_input()
    g <- make_gson(list(top = x$ids[1:40], bot = x$ids[361:400]))
    res <- suppressWarnings(gsea_gson(x$gl, g, seed = 1, verbose = FALSE,
                                      pvalueCutoff = 1, reportNA = TRUE,
                                      eps = 1))
    expect_true(all(res@result$pvalue >= 0 & res@result$pvalue <= 1, na.rm = TRUE))
})

test_that("D4: an all-zero ranking is an error, not a table of NaN", {
    skip_without_gson()
    x <- make_input()
    g <- make_gson(list(top = x$ids[1:40]))
    zero <- stats::setNames(rep(0, length(x$ids)), x$ids)

    ## Match the input-validation message, not the scaler's: both mention "all
    ## statistics are zero", so the looser regexp passed even with the
    ## pre-validation removed -- the scaler's fallback fired instead. Only the
    ## input-validation wording distinguishes which layer caught it.
    expect_error(
        suppressWarnings(gsea_gson(zero, g, verbose = FALSE, pvalueCutoff = 1)),
        "ranking has no order"
    )
    ## the scaler guards a direct caller too
    expect_error(enrichit:::scale_fgsea_ranks(zero, 1), "total weight")
    ## a ranking that is zero except for one gene is still a valid ranking
    one <- zero
    one[1] <- 1
    expect_s4_class(
        suppressWarnings(gsea_gson(sort(one, decreasing = TRUE), g, seed = 1,
                                   verbose = FALSE, pvalueCutoff = 1)),
        "gseaResult"
    )
})

test_that("D3: a gene set covering the whole ranked list is excluded, not NaN", {
    skip_without_gson()
    x <- make_input()
    g <- make_gson(list(top = x$ids[1:40], all = x$ids))

    res <- suppressWarnings(gsea_gson(x$gl, g, seed = 1, verbose = FALSE,
                                      pvalueCutoff = 1, reportNA = TRUE))
    expect_identical(res@result$ID, "top")
    expect_false(any(is.nan(res@result$enrichmentScore)))

    ## a set covering n - 1 genes is still tested (fgsea's maxSize boundary)
    g2 <- make_gson(list(near = x$ids[1:(length(x$ids) - 1)]))
    res2 <- suppressWarnings(gsea_gson(x$gl, g2, seed = 1, verbose = FALSE,
                                       pvalueCutoff = 1, reportNA = TRUE,
                                       minGSSize = 1, maxGSSize = 1e6))
    expect_identical(res2@result$ID, "near")
})

test_that("D3: when only whole-list sets remain, the result is NULL", {
    skip_without_gson()
    x <- make_input()
    g <- make_gson(list(all = x$ids, all2 = x$ids))
    res <- suppressWarnings(gsea_gson(x$gl, g, seed = 1, verbose = FALSE,
                                      pvalueCutoff = 1, minGSSize = 1,
                                      maxGSSize = 1e6))
    expect_null(res)
})
