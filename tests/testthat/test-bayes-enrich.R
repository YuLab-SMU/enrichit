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

test_that("bayes_enrich warns when no term reaches the posterior cutoff", {
    # reproduce the original dead-default regime: gene sets large relative to
    # the hit list, scored with the old false negative rate
    gene_sets <- lapply(seq_len(6), function(i) {
        paste0("g", seq.int((i - 1) * 50 + 1, i * 50))
    })
    names(gene_sets) <- paste0("T", 1:6)
    # ten hits in each of the first three sets, none in the rest
    gene <- as.character(unlist(lapply(c(0, 50, 100), function(o) paste0("g", o + 1:10))))
    universe <- paste0("g", 1:2000)
    x <- make_enrich_result(gene_sets, gene, universe)

    expect_warning(
        y <- bayes_enrich(x, n_terms = Inf, false_negative = 0.1,
                          n_iter = 200, burnin = 50, seed = 3),
        "no candidate term reached the posterior cutoff"
    )
    expect_true(sum(as.data.frame(y)$bayes_active, na.rm = TRUE) == 0L)

    # the same input with the realistic rate must select terms and stay silent
    expect_silent(
        z <- bayes_enrich(x, n_terms = Inf, false_negative = 0.9,
                          n_iter = 200, burnin = 50, seed = 3)
    )
    expect_true(sum(as.data.frame(z)$bayes_active, na.rm = TRUE) > 0L)
})

test_that("unmodeled rows keep bayes_active = NA, not FALSE", {
    gene_sets <- lapply(seq_len(8), function(i) {
        paste0("g", seq.int((i - 1) * 10 + 1, i * 10))
    })
    names(gene_sets) <- paste0("T", 1:8)
    gene <- paste0("g", c(1:5, 11:15, 21:25))
    universe <- paste0("g", 1:200)
    # all eight terms are part of the reported result, but only five are
    # candidates, so three rows must come back without a posterior
    x <- make_enrich_result(gene_sets, gene, universe,
                            pvalueCutoff = 1, qvalueCutoff = 1)

    y <- suppressWarnings(bayes_enrich(x, candidate = "all", n_terms = 5,
                                       n_iter = 200, burnin = 50, seed = 2))
    res <- as.data.frame(y)

    unmodeled <- is.na(res$posterior)
    expect_true(any(unmodeled))                     # the test premise
    expect_true(all(is.na(res$bayes_active[unmodeled])))
    expect_true(all(!is.na(res$bayes_active[!unmodeled])))

    # bayes_summary(active = TRUE) must not surface unmodeled rows
    active <- bayes_summary(y, active = TRUE)
    expect_true(all(!is.na(active$bayes_active)))
    expect_true(all(active$posterior >= 0.5))
})

test_that("bayes_enrich rejects a burnin/thin combination that keeps no sample", {
    gene_sets <- list(T1 = c("g1", "g2", "g3"))
    x <- make_enrich_result(gene_sets, gene = c("g1", "g2"),
                            universe = paste0("g", 1:10))

    # passes the old checks (burnin < n_iter, thin >= 1) yet leaves no
    # post-burn-in sample, which used to divide by zero and return NaN
    expect_error(
        bayes_enrich(x, n_terms = Inf, n_iter = 101, burnin = 100, thin = 2),
        "exceeds n_iter"
    )
    # one iteration later and a sample is kept again
    expect_silent(
        bayes_enrich(x, n_terms = Inf, n_iter = 102, burnin = 100, thin = 2, seed = 1)
    )
})
