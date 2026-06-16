#' Bayesian term selection for enrichment results
#'
#' `bayes_enrich()` adds a model-based selection layer on top of ORA results.
#' It estimates the posterior probability that each candidate term is an
#' active biological program explaining the observed input genes.
#'
#' The implementation uses a lightweight Metropolis-Hastings sampler over
#' binary latent term states. Given active terms, each gene is modeled as
#' observed with probability `1 - false_negative` if covered by at least one
#' active term, and with probability `false_positive` otherwise. The prior
#' probability that a candidate term is active is `prior`.
#'
#' This is intended as a result-compression and interpretation layer, not as a
#' replacement for ORA p-values.
#'
#' @param x An `enrichResult` object, typically from `ora_gson()` or a package
#'   that builds on enrichit, such as clusterProfiler.
#' @param candidate Candidate terms to include. `"significant"` uses
#'   `as.data.frame(x)`, `"all"` uses `x@result`, `"top"` uses the top
#'   `n_terms` rows from `x@result` ordered by `by`; or provide a character
#'   vector of term IDs.
#' @param n_terms Maximum number of candidate terms when `candidate = "top"` or
#'   when more candidates are supplied than this value. Use `Inf` to disable.
#' @param by Column used to order candidate terms.
#' @param prior Prior probability that a term is active.
#' @param false_positive Probability of observing a gene not covered by active
#'   terms.
#' @param false_negative Probability of missing a gene covered by active terms.
#' @param n_iter Total number of MCMC iterations.
#' @param burnin Number of initial iterations discarded.
#' @param thin Keep one sample every `thin` iterations after burn-in.
#' @param posterior_cutoff Terms with posterior greater than or equal to this
#'   value are marked active.
#' @param seed Optional random seed.
#' @param verbose Print sampler progress.
#' @return The input `enrichResult` object with additional columns in
#'   `@result`: `posterior`, `posterior_odds`, `bayes_rank`,
#'   `bayes_active`, `bayes_covered_gene`, and `bayes_covered_count`.
#' @export
bayes_enrich <- function(x,
                         candidate = c("top", "significant", "all"),
                         n_terms = 200,
                         by = "p.adjust",
                         prior = 0.1,
                         false_positive = 0.01,
                         false_negative = 0.1,
                         n_iter = 5000,
                         burnin = 1000,
                         thin = 1,
                         posterior_cutoff = 0.5,
                         seed = NULL,
                         verbose = FALSE) {
    if (!is(x, "enrichResult")) {
        stop("x should be an enrichResult object.", call. = FALSE)
    }

    .validate_bayes_enrich_args(
        prior = prior,
        false_positive = false_positive,
        false_negative = false_negative,
        n_iter = n_iter,
        burnin = burnin,
        thin = thin,
        posterior_cutoff = posterior_cutoff
    )

    all_res <- x@result
    if (nrow(all_res) == 0L) {
        return(x)
    }

    selected <- .select_bayes_terms(
        x = x,
        candidate = candidate,
        n_terms = n_terms,
        by = by
    )
    term_ids <- selected$term_ids

    gene_sets <- x@geneSets[term_ids]
    gene_sets <- gene_sets[!vapply(gene_sets, is.null, logical(1))]
    gene_sets <- lapply(gene_sets, as.character)
    term_ids <- names(gene_sets)

    if (length(term_ids) == 0L) {
        stop("No candidate terms have gene sets in x@geneSets.", call. = FALSE)
    }

    observed <- unique(as.character(x@gene))
    observed <- observed[!is.na(observed) & observed != ""]
    if (length(observed) == 0L) {
        stop("No observed genes available in x@gene.", call. = FALSE)
    }

    universe <- unique(as.character(x@universe))
    universe <- universe[!is.na(universe) & universe != ""]
    if (length(universe) == 0L) {
        universe <- unique(c(observed, unlist(gene_sets, use.names = FALSE)))
    }

    model_input <- .build_bayes_enrich_input(
        gene_sets = gene_sets,
        observed = observed,
        universe = universe
    )

    if (length(model_input$terms) == 0L) {
        stop("No candidate term overlaps the analysis universe.", call. = FALSE)
    }

    posterior <- .sample_active_terms(
        term_genes = model_input$term_genes,
        observed = model_input$observed,
        prior = prior,
        false_positive = false_positive,
        false_negative = false_negative,
        n_iter = n_iter,
        burnin = burnin,
        thin = thin,
        seed = seed,
        verbose = verbose
    )

    active <- posterior >= posterior_cutoff
    odds <- posterior / pmax(1 - posterior, .Machine$double.eps)
    ranks <- rank(-posterior, ties.method = "first")
    covered <- .posterior_covered_genes(
        term_genes = model_input$term_genes,
        observed = model_input$observed,
        genes = model_input$genes
    )

    out <- data.frame(
        ID = model_input$terms,
        posterior = posterior,
        posterior_odds = odds,
        bayes_rank = ranks,
        bayes_active = active,
        bayes_covered_gene = covered$gene,
        bayes_covered_count = covered$count,
        stringsAsFactors = FALSE
    )

    res <- merge(all_res, out, by = "ID", all.x = TRUE, sort = FALSE)
    ord <- match(all_res$ID, res$ID)
    res <- res[ord, , drop = FALSE]

    x@result <- res
    x@method <- if (length(x@method)) {
        paste0(x@method, "; bayes_enrich")
    } else {
        "bayes_enrich"
    }
    x
}

#' Summarize Bayesian enrichment results
#'
#' Return a data frame sorted by posterior probability from a result processed
#' by `bayes_enrich()`. This is a convenience wrapper around sorting
#' `as.data.frame(x)` by decreasing `posterior`.
#'
#' @param x An `enrichResult` object processed by `bayes_enrich()`.
#' @param active Logical. If `TRUE`, keep only rows with `bayes_active = TRUE`.
#' @param n Number of rows to return. Use `Inf` to return all rows.
#' @return A data frame ordered by decreasing `posterior`.
#' @export
bayes_summary <- function(x, active = FALSE, n = Inf) {
    if (!is(x, "enrichResult")) {
        stop("x should be an enrichResult object.", call. = FALSE)
    }

    res <- as.data.frame(x)
    if (!"posterior" %in% names(res)) {
        stop("x does not contain a posterior column. Run bayes_enrich(x) first.",
             call. = FALSE)
    }

    if (isTRUE(active)) {
        if (!"bayes_active" %in% names(res)) {
            stop("x does not contain a bayes_active column. Run bayes_enrich(x) first.",
                 call. = FALSE)
        }
        res <- res[!is.na(res$bayes_active) & res$bayes_active, , drop = FALSE]
    }

    if ("p.adjust" %in% names(res)) {
        ord <- order(-res$posterior, res$p.adjust, na.last = TRUE)
    } else {
        ord <- order(-res$posterior, na.last = TRUE)
    }
    res <- res[ord, , drop = FALSE]

    if (is.finite(n) && nrow(res) > n) {
        res <- res[seq_len(n), , drop = FALSE]
    }
    rownames(res) <- NULL
    res
}

.validate_bayes_enrich_args <- function(prior,
                                        false_positive,
                                        false_negative,
                                        n_iter,
                                        burnin,
                                        thin,
                                        posterior_cutoff) {
    probs <- c(
        prior = prior,
        false_positive = false_positive,
        false_negative = false_negative,
        posterior_cutoff = posterior_cutoff
    )
    if (any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
        stop("prior, false_positive, false_negative, and posterior_cutoff must be in (0, 1).",
             call. = FALSE)
    }
    if (n_iter <= 1L || burnin < 0L || burnin >= n_iter || thin < 1L) {
        stop("Require n_iter > 1, 0 <= burnin < n_iter, and thin >= 1.",
             call. = FALSE)
    }
}

.select_bayes_terms <- function(x, candidate, n_terms, by) {
    candidate_choices <- c("top", "significant", "all")
    all_res <- x@result
    if (is.character(candidate) &&
        (length(candidate) > 1L && !all(candidate %in% candidate_choices) ||
            length(candidate) == 1L && !candidate %in% candidate_choices)) {
        term_ids <- intersect(candidate, all_res$ID)
    } else {
        candidate <- match.arg(candidate, candidate_choices)
        res <- switch(candidate,
            significant = as.data.frame(x),
            all = all_res,
            top = all_res
        )
        term_ids <- res$ID
        if (candidate == "top") {
            if (!by %in% names(res)) {
                stop("Column specified by 'by' not found in enrichment result.", call. = FALSE)
            }
            vals <- res[[by]]
            ord <- if (is.numeric(vals)) order(vals, na.last = TRUE) else order(as.character(vals), na.last = TRUE)
            term_ids <- res$ID[ord]
        }
    }

    if (is.finite(n_terms) && length(term_ids) > n_terms) {
        term_ids <- term_ids[seq_len(n_terms)]
    }
    list(term_ids = unique(as.character(term_ids)))
}

.build_bayes_enrich_input <- function(gene_sets, observed, universe) {
    universe <- unique(universe)
    observed <- intersect(unique(observed), universe)
    if (length(observed) == 0L) {
        stop("No observed genes overlap the analysis universe.", call. = FALSE)
    }

    gene_sets <- lapply(gene_sets, function(gs) intersect(unique(gs), universe))
    keep <- vapply(gene_sets, length, integer(1)) > 0L
    gene_sets <- gene_sets[keep]

    terms <- names(gene_sets)
    genes <- unique(c(universe, observed))
    gene_index <- stats::setNames(seq_along(genes), genes)

    term_genes <- lapply(gene_sets, function(gs) unname(gene_index[gs]))
    observed_index <- unname(gene_index[observed])

    list(
        terms = terms,
        genes = genes,
        term_genes = term_genes,
        observed = observed_index
    )
}

.sample_active_terms <- function(term_genes,
                                 observed,
                                 prior,
                                 false_positive,
                                 false_negative,
                                 n_iter,
                                 burnin,
                                 thin,
                                 seed,
                                 verbose) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n_terms <- length(term_genes)
    n_genes <- max(unlist(term_genes, use.names = FALSE), observed)
    y <- logical(n_genes)
    y[observed] <- TRUE

    active <- stats::rbinom(n_terms, size = 1L, prob = prior) == 1L
    cover_count <- integer(n_genes)
    if (any(active)) {
        for (idx in which(active)) {
            cover_count[term_genes[[idx]]] <- cover_count[term_genes[[idx]]] + 1L
        }
    }

    log_prior_on <- log(prior)
    log_prior_off <- log1p(-prior)
    log_fp <- log(false_positive)
    log_tn <- log1p(-false_positive)
    log_tp <- log1p(-false_negative)
    log_fn <- log(false_negative)

    gene_loglik <- function(cc, idx = NULL) {
        yy <- if (is.null(idx)) y else y[idx]
        covered <- cc > 0L
        ifelse(yy & covered, log_tp,
            ifelse(yy & !covered, log_fp,
                ifelse(!yy & covered, log_fn, log_tn)))
    }

    current_gene_score <- gene_loglik(cover_count)
    kept <- integer(n_terms)
    n_kept <- 0L

    for (iter in seq_len(n_iter)) {
        j <- sample.int(n_terms, 1L)
        affected <- term_genes[[j]]
        proposal_cover <- cover_count
        if (!active[j]) {
            proposal_cover[affected] <- proposal_cover[affected] + 1L
            prior_delta <- log_prior_on - log_prior_off
        } else {
            proposal_cover[affected] <- proposal_cover[affected] - 1L
            prior_delta <- log_prior_off - log_prior_on
        }

        proposal_gene_score <- gene_loglik(proposal_cover[affected], affected)
        score_delta <- sum(proposal_gene_score - current_gene_score[affected]) +
            prior_delta

        if (log(stats::runif(1L)) < score_delta) {
            active[j] <- !active[j]
            cover_count <- proposal_cover
            current_gene_score[affected] <- proposal_gene_score
        }

        if (iter > burnin && ((iter - burnin) %% thin == 0L)) {
            kept <- kept + as.integer(active)
            n_kept <- n_kept + 1L
        }

        if (verbose && iter %% 1000L == 0L) {
            message("bayes_enrich iteration ", iter, "/", n_iter)
        }
    }

    kept / n_kept
}

.posterior_covered_genes <- function(term_genes, observed, genes) {
    gene <- character(length(term_genes))
    count <- integer(length(term_genes))
    for (i in seq_along(term_genes)) {
        covered <- intersect(term_genes[[i]], observed)
        count[[i]] <- length(covered)
        gene[[i]] <- paste(genes[covered], collapse = "/")
    }
    list(gene = gene, count = count)
}
