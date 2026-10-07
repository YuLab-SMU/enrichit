## DESeq2 adapters: as_ranked() / as_genes()
##
## The mock objects are plain data.frames carrying the "DESeqResults"
## class so the S3 methods can be exercised without DESeq2 installed;
## the code paths under test only rely on `$`, `rownames()` and S3
## dispatch, which behave identically on the real S4 class.

make_mock_res <- function() {
    res <- data.frame(
        baseMean = c(100, 90, 80, 70, 60, 50, 40, 30),
        log2FoldChange = c(2, -3, 1, -1, 0.5, -0.5, 4, NA),
        lfcSE = c(0.2, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
        stat = c(8, -7, 4, -4, 2, -2, 10, NA),
        pvalue = c(1e-10, 1e-9, 0.01, 0.01, 0.2, 0.2, 1e-12, NA),
        padj = c(1e-9, 1e-8, 0.02, 0.02, 0.4, 0.4, 1e-11, NA),
        row.names = paste0("gene", 1:8)
    )
    class(res) <- c("DESeqResults", "data.frame")
    res
}

test_that("as_ranked.DESeqResults returns sorted, NA-free named vector", {
    res <- make_mock_res()
    expect_message(r <- as_ranked(res), "non-finite")
    expect_type(r, "double")
    expect_false(anyNA(r))
    expect_true(all(diff(r) <= 0))
    expect_equal(names(r), c("gene7", "gene1", "gene3", "gene5", "gene6", "gene4", "gene2"))
    expect_equal(unname(r), c(10, 8, 4, 2, -2, -4, -7))
})

test_that("as_ranked.DESeqResults honours statistic = 'log2FoldChange'", {
    res <- make_mock_res()
    r <- as_ranked(res, statistic = "log2FoldChange")
    expect_equal(names(r)[1], "gene7")
    expect_equal(max(r), 4)
    expect_true(all(diff(r) <= 0))
})

test_that("as_ranked.DESeqResults gives a helpful error when 'stat' is absent", {
    res <- make_mock_res()
    res$stat <- NULL
    expect_error(as_ranked(res), "lfcShrink")
    expect_type(as_ranked(res, statistic = "log2FoldChange"), "double")
})

test_that("as_ranked.DESeqResults requires gene identifiers", {
    l <- structure(list(stat = c(1, 2)), class = "DESeqResults")
    expect_error(as_ranked(l), "no row names")
})

test_that("as_ranked.DESeqResults errors when nothing is finite", {
    res <- make_mock_res()
    res$stat <- NA_real_
    expect_error(as_ranked(res), "no finite")
})

test_that("as_ranked.DESeqResults resolves duplicated ids by max |stat|", {
    ## data.frames forbid duplicate row names; a matrix mock allows them
    ## and exercises the same `$`/rownames() code path
    m <- matrix(
        c(8, -7, 3, 2, -3, 1),
        ncol = 2,
        dimnames = list(c("gene1", "gene1", "gene2"), c("stat", "log2FoldChange"))
    )
    class(m) <- c("DESeqResults", "matrix")
    expect_warning(r <- as_ranked(m), "duplicated")
    expect_false(anyDuplicated(names(r)) > 0)
    expect_equal(r[["gene1"]], 8)
})

test_that("as_genes.DESeqResults selects up/down/both with cutoffs", {
    res <- make_mock_res()

    both <- as_genes(res)
    expect_setequal(both, c("gene1", "gene2", "gene3", "gene4", "gene7"))

    expect_setequal(as_genes(res, direction = "up"), c("gene1", "gene3", "gene7"))
    expect_setequal(as_genes(res, direction = "down"), c("gene2", "gene4"))

    expect_setequal(as_genes(res, minAbsLFC = 1.5), c("gene1", "gene2", "gene7"))
    expect_setequal(as_genes(res, direction = "up", minAbsLFC = 3), "gene7")

    expect_setequal(as_genes(res, cutoff = 0.005), c("gene1", "gene2", "gene7"))

    expect_warning(
        empty <- as_genes(res, cutoff = 1e-50),
        "no gene has padj"
    )
    expect_length(empty, 0)
})

test_that("as_genes.DESeqResults honours strict direction signs at zero", {
    res <- make_mock_res()
    res$log2FoldChange[4] <- 0 # gene4: significant but exactly zero LFC
    expect_setequal(as_genes(res, direction = "down"), "gene2")
    expect_setequal(as_genes(res, direction = "up"), c("gene1", "gene3", "gene7"))
})

test_that("as_genes.DESeqResults reports NA padj exclusions", {
    res <- make_mock_res()
    expect_message(as_genes(res), "NA adjusted p-value")
    expect_silent(as_genes(res[1:7, ]))
})

## ---- integration with a real DESeq2 run ---------------------------------

test_that("as_ranked/as_genes work on a fitted DESeqDataSet", {
    skip_if_not_installed("DESeq2")
    set.seed(123)

    counts <- matrix(rnbinom(6 * 100, mu = 60, size = 50), nrow = 100)
    ## only genes 2:31 are up in trt; if ALL genes were shifted, the
    ## median-of-ratios size factors would absorb the global shift and
    ## no gene would be called DE
    counts[2:31, 4:6] <- counts[2:31, 4:6] * 3 + 200
    storage.mode(counts) <- "integer"
    rownames(counts) <- paste0("gene", seq_len(100))
    colnames(counts) <- paste0("s", 1:6)
    coldata <- data.frame(
        condition = factor(rep(c("ctrl", "trt"), each = 3), levels = c("ctrl", "trt")),
        row.names = paste0("s", 1:6)
    )
    ## gene1 has all-zero counts -> NA statistic in results()
    counts[1, ] <- 0L

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
    suppressWarnings(dds <- DESeq2::DESeq(dds, quiet = TRUE))

    expect_message(
        r <- as_ranked(dds),
        "non-finite"
    )
    expect_length(r, 99)
    expect_false(anyNA(r))
    expect_true(all(diff(r) <= 0))
    expect_false("gene1" %in% names(r))

    r_named <- as_ranked(dds, name = "condition_trt_vs_ctrl")
    expect_equal(names(r_named), names(r))

    r_rev <- as_ranked(dds, contrast = c("condition", "ctrl", "trt"))
    expect_equal(unname(sort(r_rev)), unname(sort(-as.numeric(r_named))), tolerance = 1e-6)

    r_lfc <- as_ranked(dds, statistic = "log2FoldChange")
    expect_false(anyNA(r_lfc))

    up <- as_genes(dds, direction = "up")
    down <- as_genes(dds, direction = "down")
    expect_type(up, "character")
    expect_type(down, "character")
    expect_gt(length(up), 0)
    expect_false("gene1" %in% up)
    expect_true(all(up %in% paste0("gene", 2:31)))

    g10 <- as_genes(dds, cutoff = 0.1)
    g05 <- as_genes(dds, cutoff = 0.05)
    expect_true(length(g10) >= length(g05))
    expect_true(all(g05 %in% g10))
})

test_that("as_ranked output feeds directly into gsea() and ora()", {
    skip_if_not_installed("DESeq2")
    set.seed(42)

    counts <- matrix(rnbinom(6 * 100, mu = 60, size = 50), nrow = 100)
    counts[2:31, 4:6] <- counts[2:31, 4:6] * 3 + 200
    storage.mode(counts) <- "integer"
    rownames(counts) <- paste0("gene", seq_len(100))
    colnames(counts) <- paste0("s", 1:6)
    coldata <- data.frame(
        condition = factor(rep(c("ctrl", "trt"), each = 3), levels = c("ctrl", "trt")),
        row.names = paste0("s", 1:6)
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
    suppressWarnings(dds <- DESeq2::DESeq(dds, quiet = TRUE))

    ranked <- as_ranked(dds)
    gene_sets <- list(
        SET_A = paste0("gene", 1:20),
        SET_B = paste0("gene", 51:70)
    )
    expect_s3_class(
        suppressWarnings(gsea(ranked, gene_sets, nPerm = 200, seed = 1, verbose = FALSE)),
        "data.frame"
    )

    sig <- as_genes(dds, direction = "up")
    expect_gt(length(sig), 0)
    expect_s3_class(
        ora(gene = sig, gene_sets = gene_sets, universe = rownames(dds)),
        "data.frame"
    )
})
