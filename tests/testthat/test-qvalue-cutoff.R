library(testthat)
library(enrichit)

## get_enriched() applies qvalueCutoff row by row. It used to be skipped
## wholesale when any qvalue was NA, so one unestimable qvalue leaked
## unfiltered rows into the result with no message
## (clusterProfiler#710, #584, #614, #332).

mk_result <- function(qvalue, qvalueCutoff = 0.02) {
  result_df <- data.frame(
    ID = c("A", "B", "C"),
    Description = c("a", "b", "c"),
    GeneRatio = c("2/3", "1/3", "1/3"),
    BgRatio = c("2/10", "4/10", "4/10"),
    RichFactor = c(1, 0.25, 0.25),
    FoldEnrichment = c(5, 1.25, 1.25),
    pvalue = c(0.001, 0.02, 0.03),
    p.adjust = c(0.002, 0.03, 0.03),
    qvalue = qvalue,
    geneID = c("G1/G2", "G3", "G4"),
    Count = c(2L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- result_df$ID
  new("enrichResult",
    result = result_df,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalueCutoff,
    organism = "human",
    ontology = "TEST",
    gene = c("G1", "G2", "G3", "G4"),
    keytype = "SYMBOL",
    universe = paste0("G", LETTERS[1:10]),
    gene2Symbol = character(0),
    geneSets = list(A = c("G1", "G2"), B = "G3", C = "G4"),
    readable = FALSE,
    termsim = matrix(0, nrow = 0, ncol = 0),
    method = "ORA",
    dr = list()
  )
}

test_that("a complete qvalue column is filtered row-wise as before", {
    r <- get_enriched(mk_result(c(0.002, 0.03, 0.03)))
    expect_identical(rownames(r@result), "A")
})

test_that("a single NA qvalue no longer disables the cutoff", {
    ## B (q = 0.03 > 0.02) must now be dropped even though C is NA;
    ## C bypasses on its own NA and stays, warned about.
    warnings <- capture_warnings(r <- get_enriched(mk_result(c(0.002, 0.03, NA))))
    expect_setequal(rownames(r@result), c("A", "C"))
    expect_true(any(grepl("bypass qvalueCutoff", warnings)))
    expect_true(any(grepl("1 of 3 terms", warnings)))
})

test_that("an all-NA qvalue column disables the cutoff with a warning", {
    warnings <- capture_warnings(r <- get_enriched(mk_result(c(NA, NA, NA))))
    expect_identical(nrow(r@result), 3L)
    expect_true(any(grepl("could not be estimated", warnings)))
    expect_true(any(grepl("Filter on p.adjust", warnings)))
})

test_that("qvalueCutoff = 1 keeps every row, NA or not", {
    warnings <- capture_warnings(
        r <- get_enriched(mk_result(c(0.002, 0.03, NA), qvalueCutoff = 1))
    )
    expect_identical(nrow(r@result), 3L)
})

test_that("the gseaResult path behaves the same", {
    make_gsea <- function(qvalue, qvalueCutoff = 0.02) {
        result_df <- data.frame(
            ID = c("Top", "Bottom"),
            Description = c("Top pathway", "Bottom pathway"),
            setSize = c(2L, 2L),
            enrichmentScore = c(0.7, -0.6),
            NES = c(1.8, -1.5),
            pvalue = c(0.01, 0.03),
            p.adjust = c(0.02, 0.03),
            qvalue = qvalue,
            rank = c(2L, 5L),
            leading_edge = rep("tags=100%, list=40%, signal=70%", 2),
            core_enrichment = c("GeneA/GeneB", "GeneD/GeneE"),
            stringsAsFactors = FALSE
        )
        rownames(result_df) <- result_df$ID
        new("gseaResult",
            result = result_df,
            organism = "human",
            setType = "TEST",
            geneSets = list(Top = c("GeneA", "GeneB"), Bottom = c("GeneD", "GeneE")),
            geneList = structure(c(2, 1, -1, -2), .Names = paste0("Gene", LETTERS[1:4])),
            keytype = "SYMBOL",
            permScores = matrix(numeric(0), nrow = 0, ncol = 0),
            params = list(),
            gene2Symbol = character(0),
            readable = FALSE,
            termsim = matrix(0, nrow = 0, ncol = 0),
            method = "GSEA",
            dr = list()
        )
    }

    ## one NA among two: the non-NA row that fails the cutoff is dropped,
    ## the NA row bypasses with a warning
    warnings <- capture_warnings(
        r <- get_enriched(make_gsea(c(0.02, NA), qvalueCutoff = 0.01))
    )
    expect_setequal(rownames(r@result), "Bottom")
    expect_true(any(grepl("bypass qvalueCutoff", warnings)))
})

test_that("an empty result does not blame q-value estimation", {
    ## every term filtered out by pvalueCutoff: the two all-NA checks are
    ## vacuously TRUE on an empty frame, so without an nrow guard this warned
    ## "qvalue could not be estimated for any term" on every show()/dim()/
    ## as.data.frame() of a no-hit result
    obj <- mk_result(c(0.002, 0.03, 0.03))
    obj@result <- obj@result[0, , drop = FALSE]
    warnings <- capture_warnings(r <- get_enriched(obj))
    expect_identical(nrow(r@result), 0L)
    expect_length(warnings, 0L)
})

test_that("a result without a qvalue column says so, instead of blaming estimation", {
    obj <- mk_result(c(0.002, 0.03, 0.03))
    obj@result$qvalue <- NULL
    warnings <- capture_warnings(r <- get_enriched(obj))
    expect_identical(nrow(r@result), 3L)
    expect_true(any(grepl("no `qvalue` column", warnings, fixed = TRUE)))
    expect_false(any(grepl("could not be estimated", warnings, fixed = TRUE)))
})
