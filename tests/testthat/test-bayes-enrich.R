library(enrichit)

make_enrich_result <- function(gene_sets,
                               gene,
                               universe,
                               pvalueCutoff = 1,
                               qvalueCutoff = 1) {
    result <- data.frame(
        ID = names(gene_sets),
        Description = names(gene_sets),
        GeneRatio = paste0(
            vapply(gene_sets, function(x) length(intersect(x, gene)), integer(1)),
            "/",
            length(gene)
        ),
        BgRatio = vapply(gene_sets, function(x) paste0(length(x), "/", length(universe)), character(1)),
        RichFactor = 1,
        FoldEnrichment = 1,
        zScore = 0,
        pvalue = seq_along(gene_sets) / 10,
        p.adjust = seq_along(gene_sets) / 10,
        qvalue = seq_along(gene_sets) / 10,
        geneID = vapply(gene_sets, function(x) paste(intersect(x, gene), collapse = "/"), character(1)),
        Count = vapply(gene_sets, function(x) length(intersect(x, gene)), integer(1)),
        stringsAsFactors = FALSE
    )

    new(
        "enrichResult",
        result = result,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = "BH",
        qvalueCutoff = qvalueCutoff,
        organism = "UNKNOWN",
        ontology = "UNKNOWN",
        gene = as.character(gene),
        keytype = "UNKNOWN",
        universe = as.character(universe),
        geneSets = gene_sets,
        readable = FALSE
    )
}

test_that("bayes_enrich annotates enrichResult with posterior columns", {
    gene_sets <- list(
        T1 = c("g1", "g2", "g3", "g4"),
        T2 = c("g1", "g2", "g3"),
        T3 = c("g8", "g9", "g10")
    )
    gene <- c("g1", "g2", "g3")
    universe <- paste0("g", 1:10)
    x <- make_enrich_result(gene_sets, gene, universe)

    y <- bayes_enrich(
        x,
        n_terms = Inf,
        n_iter = 300,
        burnin = 100,
        seed = 1
    )

    res <- as.data.frame(y)
    expect_s4_class(y, "enrichResult")
    expect_true(all(c(
        "posterior", "posterior_odds", "bayes_rank",
        "bayes_active", "bayes_covered_gene", "bayes_covered_count"
    ) %in% names(res)))
    expect_true(all(res$posterior[!is.na(res$posterior)] >= 0))
    expect_true(all(res$posterior[!is.na(res$posterior)] <= 1))
    expect_true(res$posterior[res$ID == "T2"] > res$posterior[res$ID == "T3"])

    tab <- bayes_summary(y)
    expect_true(all(diff(tab$posterior) <= 0))
    expect_identical(names(tab), names(res))

    top1 <- bayes_summary(y, n = 1)
    expect_equal(nrow(top1), 1)

    active_tab <- bayes_summary(y, active = TRUE)
    expect_true(all(active_tab$bayes_active))
})

test_that("bayes_enrich candidate modes use the expected term space", {
    gene_sets <- list(
        T1 = c("g1", "g2", "g3"),
        T2 = c("g1", "g2"),
        T3 = c("g8", "g9")
    )
    x <- make_enrich_result(
        gene_sets,
        gene = c("g1", "g2", "g3"),
        universe = paste0("g", 1:10),
        pvalueCutoff = 0.2,
        qvalueCutoff = 0.2
    )

    selected_sig <- .select_bayes_terms(
        x,
        candidate = "significant",
        n_terms = Inf,
        by = "p.adjust"
    )
    selected_all <- .select_bayes_terms(
        x,
        candidate = "all",
        n_terms = Inf,
        by = "p.adjust"
    )

    expect_identical(selected_sig$term_ids, as.character(as.data.frame(x)$ID))
    expect_identical(selected_all$term_ids, as.character(x@result$ID))
    expect_lt(length(selected_sig$term_ids), length(selected_all$term_ids))
})

test_that("bayes_enrich validates inputs", {
    gene_sets <- list(T1 = c("g1", "g2"))
    x <- make_enrich_result(
        gene_sets,
        gene = c("g1"),
        universe = c("g1", "g2", "g3")
    )

    expect_error(
        bayes_enrich(x, prior = 1),
        "must be in"
    )
    expect_error(
        bayes_enrich(x, n_iter = 10, burnin = 10),
        "Require"
    )
    expect_error(
        bayes_summary(x),
        "Run bayes_enrich"
    )
})
