test_that("as_enrichResult builds enrichResult from canonical table", {
    df <- data.frame(
        ID = c("T1", "T2"),
        Description = c("path one", "path two"),
        GeneRatio = c("5/100", "3/100"),
        BgRatio = c("50/1000", "60/1000"),
        pvalue = c(1e-4, 0.01),
        geneID = c("g1/g2/g3/g4/g5", "g1/g2/g3"),
        Count = c(5, 3),
        stringsAsFactors = FALSE
    )
    x <- as_enrichResult(df)

    expect_s4_class(x, "enrichResult")
    expected <- c(
        "ID", "Description", "GeneRatio", "BgRatio", "pvalue",
        "p.adjust", "qvalue", "geneID", "Count"
    )
    expect_true(all(expected %in% colnames(x@result)))
    expect_equal(rownames(x@result), c("T1", "T2"))
    expect_false(anyNA(x@result$p.adjust))
    expect_false(anyNA(x@result$qvalue))
    expect_equal(x@result$GeneRatio, c("5/100", "3/100"))
    expect_equal(x@gene, c("g1", "g2", "g3", "g4", "g5"))
    expect_equal(names(x@geneSets), c("T1", "T2"))
    expect_equal(x@geneSets$T1, c("g1", "g2", "g3", "g4", "g5"))
    expect_equal(x@pvalueCutoff, 1)
    expect_equal(x@ontology, "UNKNOWN")
})

test_that("as_enrichResult derives ratios and statistics from geneSets", {
    df <- data.frame(
        ID = c("T1", "T2"),
        pvalue = c(1e-4, 0.01),
        geneID = c("g1/g2/g3/g4/g5", "g1/g2/g3"),
        stringsAsFactors = FALSE
    )
    geneSets <- list(
        T1 = paste0("g", 1:50),
        T2 = paste0("g", 51:110)
    )
    universe <- paste0("g", 1:1000)
    x <- as_enrichResult(df, geneSets = geneSets, gene = paste0("g", 1:100), universe = universe)

    expect_equal(x@result$GeneRatio, c("5/100", "3/100"))
    expect_equal(x@result$BgRatio, c("50/1000", "60/1000"))
    expect_equal(x@result$RichFactor, c(5 / 50, 3 / 60))
    expect_equal(x@result$FoldEnrichment, c(5 / 100 / (50 / 1000), 3 / 100 / (60 / 1000)))
    expect_true("zScore" %in% colnames(x@result))
    expect_equal(x@universe, universe)
    expect_equal(x@geneSets, lapply(geneSets, unique))
})

test_that("as_enrichResult accepts common column aliases and separators", {
    df <- data.frame(
        term_id = c("GO:1", "GO:2"),
        description = c("a", "b"),
        PValue = c(0.01, 0.2),
        padj = c(0.02, 0.3),
        Genes = c("g1; g2", "g3;g4"),
        stringsAsFactors = FALSE
    )
    x <- as_enrichResult(df)
    expect_equal(x@result$ID, c("GO:1", "GO:2"))
    expect_equal(x@result$geneID, c("g1/g2", "g3/g4"))
    expect_equal(x@result$Count, c(2, 2))
    expect_equal(x@result$GeneRatio, c("2/4", "2/4"))
    expect_equal(x@result$p.adjust, c(0.02, 0.3))
})

test_that("as_enrichResult rebuilds k/n from numeric GeneRatio when gene known", {
    df <- data.frame(
        ID = c("T1"), pvalue = c(0.01), Count = c(5),
        GeneRatio = c(0.05), stringsAsFactors = FALSE
    )
    x <- as_enrichResult(df, gene = paste0("g", 1:100))
    expect_equal(x@result$GeneRatio, "5/100")
})

test_that("as_gseaResult fills rank/leading_edge/core_enrichment", {
    set.seed(1)
    stats <- rnorm(30)
    names(stats) <- paste0("g", 1:30)
    stats <- sort(stats, decreasing = TRUE)

    gs <- list(
        P1 = paste0("g", 1:10),
        P2 = paste0("g", 15:25)
    )
    es <- c(0.5, -0.4)
    df <- data.frame(
        pathway = c("P1", "P2"),
        ES = es,
        NES = c(1.8, -1.5),
        pval = c(0.01, 0.05),
        padj = c(0.02, 0.05),
        size = c(10, 11),
        stringsAsFactors = FALSE
    )
    x <- as_gseaResult(df, geneList = stats, geneSets = gs)

    expect_s4_class(x, "gseaResult")
    expect_true(all(c("rank", "leading_edge", "core_enrichment", "setSize") %in% colnames(x@result)))
    expect_true(all(x@result$rank > 0))
    expect_true(all(nzchar(x@result$leading_edge)))
    expect_false(anyNA(x@result$qvalue))
    expect_equal(x@result$setSize, c(10, 11))
    expect_equal(x@params$exponent, 1)
    expect_equal(x@geneList, stats)
    expect_equal(names(x@geneSets), c("P1", "P2"))
    ## core_enrichment genes must be members of the corresponding gene set
    core <- strsplit(x@result$core_enrichment[1], "/", fixed = TRUE)[[1]]
    expect_true(all(core %in% gs$P1))
})

test_that("as_gseaResult prefers fgsea leadingEdge list column", {
    stats <- c(g1 = 3, g2 = 2, g3 = 1, g4 = -1, g5 = -2)
    df <- data.frame(
        pathway = c("P1"),
        ES = c(0.6),
        NES = c(1.5),
        pval = c(0.01),
        padj = c(0.02),
        size = c(3),
        stringsAsFactors = FALSE
    )
    df$leadingEdge <- list(c("g1", "g2"))
    expect_warning(
        x <- as_gseaResult(df, geneList = stats),
        "core_enrichment"
    )
    expect_equal(x@result$core_enrichment, "g1/g2")
    expect_equal(names(x@geneSets), "P1")
})

test_that("as_gseaResult validates geneList", {
    df <- data.frame(ID = "P1", ES = 0.5, pvalue = 0.01, stringsAsFactors = FALSE)
    expect_error(as_gseaResult(df, geneList = c(1, 2, 3)), "named")
    stats <- c(g1 = 1, g2 = 2)
    expect_error(
        suppressWarnings(
            as_gseaResult(
                data.frame(ID = "P1", pvalue = 0.01, stringsAsFactors = FALSE),
                geneList = stats,
                geneSets = list(P1 = "g1")
            )
        ),
        "enrichmentScore"
    )
    expect_error(
        as_gseaResult(data.frame(ID = "P1", ES = 0.5, pvalue = NA_real_), geneList = stats),
        "p-value"
    )
})

test_that("as_enrichResult guards invalid input", {
    expect_error(as_enrichResult(data.frame(pvalue = 0.01)), "ID")
    expect_error(
        as_enrichResult(data.frame(ID = "T1", pvalue = NA_real_, geneID = "g1")),
        "p-value"
    )
    expect_error(as_enrichResult(data.frame(ID = character(), pvalue = numeric())), "empty")

    df <- data.frame(
        ID = c("T1", "T1"), pvalue = c(0.01, 0.02),
        geneID = c("g1", "g2"), stringsAsFactors = FALSE
    )
    expect_warning(x <- as_enrichResult(df), "duplicated")
    expect_equal(rownames(x@result), c("T1", "T1.1"))

    df2 <- data.frame(
        ID = c("T1", "T2"), pvalue = c(-0.5, 2),
        geneID = c("g1", "g2"), stringsAsFactors = FALSE
    )
    expect_warning(x2 <- as_enrichResult(df2), "clamped")
    expect_equal(x2@result$pvalue, c(0, 1))
})

test_that("as_enrichResult derives RichFactor from BgRatio when geneSets absent", {
    df <- data.frame(
        ID = "T1", GeneRatio = "5/100", BgRatio = "50/1000",
        pvalue = 0.01, geneID = "g1/g2/g3/g4/g5", stringsAsFactors = FALSE
    )
    x <- suppressWarnings(
        as_enrichResult(df, gene = paste0("g", 1:100), universe = paste0("g", 1:1000))
    )
    expect_equal(x@result$RichFactor, 5 / 50)
    expect_equal(x@result$FoldEnrichment, (5 / 100) / (50 / 1000))
})

test_that("converted results do not claim a clusterProfiler citation", {
    ora <- suppressWarnings(as_enrichResult(data.frame(
        ID = "T1", pvalue = 0.01, geneID = "g1",
        stringsAsFactors = FALSE
    )))
    gsea <- suppressWarnings(as_gseaResult(
        data.frame(ID = "T1", ES = 0.5, pvalue = 0.01,
                   core_enrichment = "g1", stringsAsFactors = FALSE),
        geneList = c(g1 = 1),
        geneSets = list(T1 = "g1")
    ))

    ora_output <- capture.output(show(ora))
    gsea_output <- capture.output(show(gsea))
    expect_false(any(grepl("clusterProfiler", ora_output, fixed = TRUE)))
    expect_false(any(grepl("clusterProfiler", gsea_output, fixed = TRUE)))
})
