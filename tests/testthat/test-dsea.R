## dsea() / as_elist(): limma competitive gene-set tests on DESeq2 input

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

test_that("dsea(camera) detects the spiked set on a DESeqDataSet", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    res <- dsea(dds, .make_sets(), contrast = c("condition", "trt", "ctrl"))

    expect_s3_class(res, "data.frame")
    expect_equal(res["SPIKED", "NGenes"], 30)
    expect_equal(res["SPIKED", "Direction"], "Up")
    expect_lt(res["SPIKED", "PValue"], 0.05)
    expect_lt(res["SPIKED", "PValue"], res["NULLSET", "PValue"])
})

test_that("dsea accepts an EList and all contrast forms", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)
    sets <- .make_sets()

    base <- dsea(dds, sets, contrast = c("condition", "trt", "ctrl"))
    expect_equal(
        dsea(elist, sets, contrast = c("condition", "trt", "ctrl")),
        base
    )
    expect_equal(dsea(dds, sets, contrast = "condition_trt_vs_ctrl"), base)
    ## numeric contrast vectors are a limma-backend feature
    expect_equal(
        dsea(dds, sets, contrast = c(0, 1), engine = "limma"),
        dsea(dds, sets, contrast = c("condition", "trt", "ctrl"),
                 engine = "limma")
    )
    expect_error(
        dsea(dds, sets, contrast = c(0, 1)),
        "engine = \"native\""
    )
    expect_equal(dsea(dds, sets), base)
})

test_that("roast and fry backends run on the EList path", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    sets <- .make_sets()

    rr <- dsea(dds, sets, contrast = "condition_trt_vs_ctrl", method = "roast")
    expect_equal(nrow(rr), 2)

    rf <- dsea(dds, sets, contrast = "condition_trt_vs_ctrl", method = "fry")
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
    gdelta <- dsea(elist, sets, contrast = "condition_trt_vs_ctrl")
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
    expect_error(dsea(1, list(A = "x")), "no as_elist")

    dds <- .make_dds_set()
    expect_error(
        dsea(dds, list(S = paste0("gene", 2:31)),
                 contrast = c("condition", "A", "B")),
        "maps to design column"
    )
    expect_error(
        dsea(dds, list(X = "nope1", Y = "nope2")),
        "none of the gene sets"
    )
    expect_warning(
        res <- dsea(
            dds, list(SPIKED = paste0("gene", 2:31), TINY = "gene1")
        ),
        "dropped"
    )
    expect_equal(nrow(res), 1)

    ## unnamed single vector is treated as one set
    res1 <- dsea(dds, paste0("gene", 2:31))
    expect_equal(nrow(res1), 1)
})

test_that("native camera matches limma (fixed, estimated, unweighted)", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)
    sets <- .make_sets()
    ct <- "condition_trt_vs_ctrl"

    nat <- dsea(elist, sets, contrast = ct)
    lim <- dsea(elist, sets, contrast = ct, engine = "limma")
    expect_identical(names(nat), names(lim))
    expect_equal(nat, lim, tolerance = 1e-8)

    nat_est <- dsea(elist, sets, contrast = ct, inter.gene.cor = NA)
    lim_est <- dsea(elist, sets, contrast = ct, engine = "limma",
                        inter.gene.cor = NA)
    expect_identical(names(nat_est), names(lim_est))
    expect_equal(nat_est, lim_est, tolerance = 1e-8)

    elist0 <- elist
    elist0$weights <- NULL
    nat_uw <- dsea(elist0, sets, contrast = ct)
    lim_uw <- dsea(elist0, sets, contrast = ct, engine = "limma")
    expect_equal(nat_uw, lim_uw, tolerance = 1e-8)
})

test_that("native fry matches limma (posterior.sd and p2)", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)
    sets <- .make_sets()
    ct <- "condition_trt_vs_ctrl"

    nat <- dsea(elist, sets, contrast = ct, method = "fry")
    lim <- dsea(elist, sets, contrast = ct, method = "fry",
                    engine = "limma")
    expect_identical(names(nat), names(lim))
    expect_equal(nat, lim, tolerance = 1e-8)

    nat_p2 <- dsea(elist, sets, contrast = ct, method = "fry",
                       standardize = "p2")
    lim_p2 <- dsea(elist, sets, contrast = ct, method = "fry",
                       engine = "limma", standardize = "p2")
    expect_equal(nat_p2, lim_p2, tolerance = 1e-8)
})

test_that("native roast reproduces limma bit-exactly under the same seed", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)
    sets <- .make_sets()
    ct <- "condition_trt_vs_ctrl"
    gw <- rep_len(c(1, -1), length.out = nrow(elist$E))

    for (ss in c("mean", "floormean", "mean50", "msq")) {
        set.seed(2026)
        nat <- dsea(elist, sets, contrast = ct, method = "roast",
                        nrot = 99, set.statistic = ss)
        set.seed(2026)
        lim <- dsea(elist, sets, contrast = ct, method = "roast",
                        nrot = 99, set.statistic = ss, engine = "limma")
        expect_identical(names(nat), names(lim), label = ss)
        expect_equal(nat, lim, tolerance = 0, label = ss)
    }

    set.seed(2026)
    nat_gw <- dsea(elist, sets, contrast = ct, method = "roast",
                       nrot = 99, gene.weights = gw)
    set.seed(2026)
    lim_gw <- dsea(elist, sets, contrast = ct, method = "roast",
                       nrot = 99, gene.weights = gw, engine = "limma")
    expect_equal(nat_gw, lim_gw, tolerance = 0)

    elist0 <- elist
    elist0$weights <- NULL
    set.seed(2026)
    nat_uw <- dsea(elist0, sets, contrast = ct, method = "roast",
                       nrot = 99)
    set.seed(2026)
    lim_uw <- dsea(elist0, sets, contrast = ct, method = "roast",
                       nrot = 99, engine = "limma")
    expect_equal(nat_uw, lim_uw, tolerance = 0)
})

test_that("native roast honors midp", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)
    sets <- .make_sets()
    ct <- "condition_trt_vs_ctrl"

    set.seed(7)
    nat <- dsea(elist, sets, contrast = ct, method = "roast",
                    nrot = 99, midp = FALSE)
    set.seed(7)
    lim <- dsea(elist, sets, contrast = ct, method = "roast",
                    nrot = 99, midp = FALSE, engine = "limma")
    expect_equal(nat, lim, tolerance = 0)
    expect_error(
        dsea(elist, sets, contrast = ct, method = "roast",
                 adjust.method = "holm"),
        "not supported"
    )
})

test_that("single-set native outputs match limma shapes", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)
    one <- paste0("gene", 2:31)
    ct <- "condition_trt_vs_ctrl"

    expect_false("FDR" %in% names(dsea(elist, one, contrast = ct)))
    expect_false("FDR" %in% names(dsea(elist, one, contrast = ct,
                                           engine = "limma")))
    rf <- dsea(elist, one, contrast = ct, method = "fry")
    expect_false("FDR" %in% names(rf))
    expect_true("PValue.Mixed" %in% names(rf))
    rf_lim <- dsea(elist, one, contrast = ct, method = "fry",
                       engine = "limma")
    expect_equal(rf, rf_lim, tolerance = 1e-8)

    set.seed(11)
    nat <- dsea(elist, one, contrast = ct, method = "roast", nrot = 99)
    set.seed(11)
    lim <- dsea(elist, one, contrast = ct, method = "roast", nrot = 99,
                    engine = "limma")
    expect_equal(nat, lim, tolerance = 0)
})

test_that("native engine rejects unsupported limma arguments", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)
    sets <- .make_sets()
    ct <- "condition_trt_vs_ctrl"

    expect_error(
        dsea(elist, sets, contrast = ct, use.ranks = TRUE),
        "not supported by dsea\\(method = \"camera\", engine = \"native\"\\)"
    )
    expect_error(
        dsea(elist, sets, contrast = ct, method = "roast",
                 trend.var = TRUE),
        "not supported by dsea\\(method = \"roast\", engine = \"native\"\\)"
    )
    expect_error(
        dsea(elist, sets, contrast = ct, method = "fry",
                 trend.var = TRUE),
        "not supported by dsea\\(method = \"fry\", engine = \"native\"\\)"
    )
    expect_error(
        dsea(elist, sets, contrast = ct, method = "camera",
                 engine = "native", TRUE),
        "must be named"
    )
    expect_error(
        dsea(elist, sets, contrast = "no_such_column"),
        "not found in the design"
    )
})

test_that("native roast msq reproduces limma across chunks with gene.weights", {
    skip_if_not_installed("limma")
    dds <- .make_dds_set(42)
    elist <- as_elist(dds)
    sets <- .make_sets()
    ct <- "condition_trt_vs_ctrl"
    gw <- rep_len(c(2, -1), nrow(elist$E))

    set.seed(2026)
    nat <- dsea(elist, sets, contrast = ct, method = "roast",
                    nrot = 2000, set.statistic = "msq", gene.weights = gw)
    set.seed(2026)
    lim <- dsea(elist, sets, contrast = ct, method = "roast",
                    nrot = 2000, set.statistic = "msq", gene.weights = gw,
                    engine = "limma")
    expect_equal(nat, lim, tolerance = 0)
})

test_that("native roast matches limma on a two-gene dataset", {
    skip_if_not_installed("limma")
    el <- structure(list(
        E = matrix(rnorm(2 * 8), 2, 8,
                   dimnames = list(c("g1", "g2"), paste0("s", 1:8))),
        weights = NULL,
        design = cbind(Intercept = 1, cond = c(0, 0, 0, 0, 1, 1, 1, 1))
    ), class = "EList")
    one <- list(S1 = c("g1", "g2"))

    set.seed(5)
    nat <- dsea(el, one, method = "roast", nrot = 99)
    set.seed(5)
    lim <- dsea(el, one, method = "roast", nrot = 99, engine = "limma")
    expect_equal(nat, lim, tolerance = 0)
})

test_that("native engine errors on numerically rank-deficient designs", {
    x <- c(1, -1, 1, -1, 1, -1, 1, -1)
    d <- cbind(Intercept = 1, x = x,
               x2 = x + c(1e-8, 0, 2e-8, 1e-8, 0, 3e-8, 0, 1e-8))
    el <- structure(list(
        E = matrix(rnorm(50 * 8), 50, 8,
                   dimnames = list(paste0("g", 1:50), paste0("s", 1:8))),
        weights = NULL, design = d
    ), class = "EList")
    sets <- list(S1 = paste0("g", 1:10), S2 = paste0("g", 11:25))

    expect_error(dsea(el, sets, contrast = 2), "not of full rank")
    expect_error(dsea(el, sets, contrast = 2, method = "roast",
                          nrot = 99), "not of full rank")
    expect_error(dsea(el, sets, contrast = 2, method = "fry"),
                 "not of full rank")
})
