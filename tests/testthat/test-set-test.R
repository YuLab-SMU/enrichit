## set_test() / as_elist(): limma competitive gene-set tests on DESeq2 input

.make_dds_set <- function(seed = 123) {
    skip_if_not_installed("DESeq2")
    set.seed(seed)
    counts <- matrix(rnbinom(6 * 100, mu = 60, size = 50), nrow = 100)
    ## only genes 2:31 are up in trt; shifting ALL genes would be absorbed
    ## by the median-of-ratios size factors (see test-deseq2.R)
    counts[2:31, 4:6] <- counts[2:31, 4:6] * 3 + 200
    storage.mode(counts) <- "integer"
    rownames(counts) <- paste0("gene", seq_len(100))
    colnames(counts) <- paste0("s", 1:6)
    coldata <- data.frame(
        condition = factor(rep(c("ctrl", "trt"), each = 3), levels = c("ctrl", "trt")),
        row.names = paste0("s", 1:6)
    )
    ## gene1 all-zero -> mu = 0, dispersion NA
    counts[1, ] <- 0L
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts, colData = coldata, design = ~ condition
    )
    suppressWarnings(DESeq2::DESeq(dds, quiet = TRUE))
}

.make_sets <- function() {
    list(
        SPIKED = paste0("gene", 2:31),
        NULLSET = paste0("gene", 51:80)
    )
}

test_that("as_elist builds voom-style E with delta-method precision weights", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set()
    elist <- as_elist(dds)

    expect_s4_class(elist, "EList")
    expect_identical(dim(elist$E), c(100L, 6L))
    expect_identical(dimnames(elist$E), dimnames(elist$weights))
    ## DESeq2's stored model matrix uses resultsNames() conventions:
    ## 'Intercept' (no parentheses), 'condition_trt_vs_ctrl'
    expect_identical(
        colnames(elist$design),
        c("Intercept", "condition_trt_vs_ctrl")
    )

    ct <- DESeq2::counts(dds)
    lib <- colSums(ct) * DESeq2::sizeFactors(dds)
    E_expected <- log2((ct + 0.5) / (matrix(lib, 100, 6, byrow = TRUE) + 1) * 1e6)
    expect_equal(elist$E, E_expected, tolerance = 1e-12)

    ## scalar hand-check of w = (mu_f + 0.5)^2 log(2)^2 / (mu_f + alpha mu_f^2)
    i <- 2
    j <- 3
    m <- SummarizedExperiment::assay(dds, "mu")[i, j]
    a <- SummarizedExperiment::mcols(dds)$dispersion[i]
    mf <- max(m, 1)
    expect_equal(
        elist$weights[i, j],
        (mf + 0.5)^2 * log(2)^2 / (mf + a * mf^2),
        tolerance = 1e-10
    )

    ## camera requires strictly positive finite weights, incl. all-zero gene1
    expect_true(all(elist$weights > 0))
    expect_true(all(is.finite(elist$weights)))
})

test_that("precision weights increase with fitted abundance", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set()
    elist <- as_elist(dds)
    mu_g <- rowMeans(SummarizedExperiment::assay(dds, "mu"))
    w_g <- rowMeans(elist$weights)
    ## gene1 (all-zero) has NA fitted means
    ok <- is.finite(mu_g)
    expect_gt(cor(mu_g[ok], w_g[ok]), 0.5)
})

test_that("mu.floor and NA handling in the weight helper", {
    mu <- matrix(c(NA, 5, 0, 100, 2, 50), nrow = 3, byrow = TRUE)
    alpha <- c(0.1, 0.2, 0.3)
    expect_message(
        w <- enrichit:::.deseq2_precision_weights(mu, alpha, mu.floor = 1),
        "non-finite"
    )
    ln2sq <- log(2)^2
    ## mu = 5, alpha = 0.1
    expect_equal(w[1, 2], 5.5^2 * ln2sq / (5 + 0.1 * 25), tolerance = 1e-12)
    ## mu = 0 is floored to 1, alpha = 0.2
    expect_equal(w[2, 1], 1.5^2 * ln2sq / (1 + 0.2), tolerance = 1e-12)
    ## mu = 2, alpha = 0.3
    expect_equal(w[3, 1], 2.5^2 * ln2sq / (2 + 0.3 * 4), tolerance = 1e-12)
    ## NA fitted mean -> minimum finite weight
    expect_equal(w[1, 1], min(w[is.finite(w)]), tolerance = 1e-12)
    expect_true(all(w > 0 & is.finite(w)))
})

test_that("mu.floor applies at the dds level", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set()
    elist <- as_elist(dds, mu.floor = 10)
    elist1 <- as_elist(dds)

    ## observations with fitted means above the floor are unchanged
    mu <- SummarizedExperiment::assay(dds, "mu")
    min_finite_mu <- vapply(
        seq_len(nrow(mu)),
        function(i) {
            x <- mu[i, ]
            if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
        },
        numeric(1)
    )
    hi <- which(!is.na(min_finite_mu) & min_finite_mu > 10)
    expect_true(length(hi) > 50)
    expect_false("gene1" %in% rownames(mu)[hi])
    expect_equal(elist$weights[hi, ], elist1$weights[hi, ], tolerance = 1e-12)

    ## all-zero gene1 (NA mu) carries the minimum finite weight
    expect_equal(
        elist1$weights["gene1", 1],
        min(elist1$weights[is.finite(elist1$weights)]),
        tolerance = 1e-12
    )
    expect_error(as_elist(dds, mu.floor = -1), "mu.floor")
})

test_that("set_test(camera) detects the spiked set on a DESeqDataSet", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    res <- set_test(dds, .make_sets(), contrast = c("condition", "trt", "ctrl"))

    expect_s3_class(res, "data.frame")
    expect_equal(res["SPIKED", "NGenes"], 30)
    expect_equal(res["SPIKED", "Direction"], "Up")
    expect_lt(res["SPIKED", "PValue"], 0.05)
    expect_lt(res["SPIKED", "PValue"], res["NULLSET", "PValue"])
})

test_that("set_test accepts an EList and all contrast forms", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)
    sets <- .make_sets()

    base <- set_test(dds, sets, contrast = c("condition", "trt", "ctrl"))
    expect_equal(
        set_test(elist, sets, contrast = c("condition", "trt", "ctrl")),
        base
    )
    expect_equal(set_test(dds, sets, contrast = "condition_trt_vs_ctrl"), base)
    expect_equal(set_test(dds, sets, contrast = c(0, 1)), base)
    expect_equal(set_test(dds, sets), base)
})

test_that("roast and fry backends run on the EList path", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    sets <- .make_sets()

    rr <- set_test(dds, sets, contrast = "condition_trt_vs_ctrl", method = "roast")
    expect_equal(nrow(rr), 2)

    rf <- set_test(dds, sets, contrast = "condition_trt_vs_ctrl", method = "fry")
    expect_s3_class(rf, "data.frame")
    expect_equal(nrow(rf), 2)
})

test_that("delta weights agree with voom weights in ranking", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)

    ct <- DESeq2::counts(dds)
    v <- limma::voom(
        ct,
        design = elist$design,
        lib.size = colSums(ct) * DESeq2::sizeFactors(dds)
    )
    expect_gt(
        cor(as.numeric(elist$weights), as.numeric(v$weights), method = "spearman"),
        0.3
    )
    ## rank correlation is diluted by per-gene dispersion heterogeneity:
    ## the delta weights vary with alpha_i within an abundance cluster,
    ## voom's smooth trend cannot capture that -- both are valid weights

    sets <- .make_sets()
    gv <- limma::camera(v, sets, contrast = "condition_trt_vs_ctrl")
    gdelta <- set_test(elist, sets, contrast = "condition_trt_vs_ctrl")
    expect_lt(gv["SPIKED", "PValue"], 0.05)
    expect_lt(gdelta["SPIKED", "PValue"], 0.05)
})

test_that("informative errors and set filtering", {
    skip_if_not_installed("DESeq2")
    skip_if_not_installed("limma")

    set.seed(1)
    counts <- matrix(rnbinom(12, mu = 60, size = 50), nrow = 2)
    storage.mode(counts) <- "integer"
    coldata <- data.frame(
        condition = factor(rep(c("ctrl", "trt"), each = 3)),
        row.names = paste0("s", 1:6)
    )
    dds0 <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts, colData = coldata, design = ~ condition
    )
    expect_error(as_elist(dds0), "DESeq")
    expect_error(as_elist(1), "no as_elist")
    expect_error(set_test(1, list(A = "x")), "no as_elist")

    dds <- .make_dds_set()
    expect_error(
        set_test(dds, list(S = paste0("gene", 2:31)),
                 contrast = c("condition", "A", "B")),
        "maps to design column"
    )
    expect_error(
        set_test(dds, list(X = "nope1", Y = "nope2")),
        "none of the gene sets"
    )
    expect_warning(
        res <- set_test(
            dds, list(SPIKED = paste0("gene", 2:31), TINY = "gene1")
        ),
        "dropped"
    )
    expect_equal(nrow(res), 1)

    ## unnamed single vector is treated as one set
    res1 <- set_test(dds, paste0("gene", 2:31))
    expect_equal(nrow(res1), 1)
})
