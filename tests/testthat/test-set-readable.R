library(testthat)
library(enrichit)

## setReadable() and EXTID2NAME() are the only exported functions that had no
## test coverage at all (2026-10-27 review, L3). setReadable() is re-exported by
## clusterProfiler, so any regression here ships to its users. The invariants
## pinned below: it changes only the display columns, never the statistics, and
## it reports failure loudly instead of silently returning unreadable IDs.

skip <- function() {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("gson")
}

make_ora <- function() {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    set.seed(11)
    u <- as.character(head(keys(org.Hs.eg.db, "ENTREZID"), 400))
    gene <- u[1:60]
    gs <- list(P1 = u[1:25], P2 = u[30:60], P3 = u[100:140])
    g <- gson::gson(
        gsid2gene = data.frame(gsid = rep(names(gs), lengths(gs)),
                               gene = unlist(gs, use.names = FALSE)),
        gsid2name = data.frame(gsid = names(gs), name = c("P one", "P two", "P three"))
    )
    ora_gson(gene, pvalueCutoff = 1, universe = u, gson = g, qvalueCutoff = 1)
}

test_that("setReadable changes display, not statistics", {
    skip()
    res <- make_ora()
    expect_false(res@readable)

    conv <- setReadable(res, org.Hs.eg.db, keyType = "ENTREZID")
    expect_true(conv@readable)

    ## every mapped geneID entry became a symbol: no bare integer remains
    ids <- unlist(strsplit(as.character(conv@result$geneID), "/", fixed = TRUE))
    ids <- ids[nzchar(ids)]
    expect_true(all(!grepl("^\\d+$", ids)))

    ## the statistics are untouched -- this is display only
    expect_identical(conv@result$pvalue, res@result$pvalue)
    expect_identical(conv@result$p.adjust, res@result$p.adjust)
    expect_identical(conv@result$Count, res@result$Count)
    expect_identical(conv@result$ID, res@result$ID)

    ## @gene stays the *query* list in its input namespace (upstream contract:
    ## setReadable rewrites the result table and @gene2Symbol, not the query
    ## record), so the assertion is that it is unchanged, not readable
    expect_identical(conv@gene, res@gene)
})

test_that("setReadable is idempotent on the statistics across repeated calls", {
    skip()
    res <- make_ora()
    once <- setReadable(res, org.Hs.eg.db, keyType = "ENTREZID")
    twice <- setReadable(once, org.Hs.eg.db, keyType = "SYMBOL")
    expect_identical(twice@result$pvalue, once@result$pvalue)
    expect_identical(twice@result$geneID, once@result$geneID)
})

test_that("EXTID2NAME maps the requested ids and tolerates unmappable ones", {
    skip()
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    out <- EXTID2NAME(org.Hs.eg.db, c("7157", "672"), "ENTREZID", "SYMBOL")
    expect_identical(unname(out["7157"]), "TP53")
    expect_identical(unname(out["672"]), "BRCA1")

    ## an id with no symbol comes back unchanged rather than as NA
    out2 <- EXTID2NAME(org.Hs.eg.db, c("7157", "1"), "ENTREZID", "SYMBOL")
    expect_false(anyNA(out2))
    expect_true(all(nzchar(out2)))
})

test_that("setReadable rejects wrong input loudly", {
    skip()
    expect_error(setReadable(data.frame(a = 1), org.Hs.eg.db), "input should be")
})
