## Pathway-score dsea (PATHWAY_DSEA_SPEC.md route A).
##
## counts + design -> VST -> per-sample pathway scores -> test on the design.
##
## The score matrix is a first-class output, not an intermediate: it can be
## reused for clustering, visualisation and heterogeneity analysis without
## re-running the enrichment.  Everything that can change a number is written
## into `meta`; a result whose provenance cannot be reconstructed from `meta` is
## a bug, not a formatting issue.

## ---- input handling --------------------------------------------------------

.dsea_counts <- function(x) {
    if (is(x, "DESeqDataSet")) {
        cnt <- DESeq2::counts(x)
    } else if (is.matrix(x) || is.data.frame(x)) {
        cnt <- as.matrix(x)
    } else {
        stop("'x' must be a DESeqDataSet or a counts matrix; got ", class(x)[1], call. = FALSE)
    }
    storage.mode(cnt) <- "double"
    if (anyNA(cnt)) stop("'x' contains missing counts.", call. = FALSE)
    if (any(cnt < 0)) stop("'x' contains negative counts.", call. = FALSE)
    ## Tolerate one-off floating point noise on counts that are meant to be
    ## integers, but refuse a genuine fractional matrix: silently rounding it
    ## would hide a mis-read input file.
    if (any(abs(cnt - round(cnt)) > 1e-8)) {
        stop("'x' looks like log-scale expression, not raw counts. ",
             "Pass the counts matrix (or a DESeqDataSet), not transformed values.",
             call. = FALSE)
    }
    cnt <- round(cnt)
    if (is.null(rownames(cnt)) || is.null(colnames(cnt))) {
        stop("'x' must have both row names (genes) and column names (samples).", call. = FALSE)
    }
    if (anyDuplicated(rownames(cnt))) stop("'x' has duplicated row names.", call. = FALSE)
    if (anyDuplicated(colnames(cnt))) stop("'x' has duplicated column names.", call. = FALSE)
    cnt
}

## ssGSEA scores are rank-based inside a sample, so duplicated members would be
## counted more than once and would quietly inflate the set's weight.  Reject
## them rather than deduplicating: the caller almost certainly did not mean it.
.dsea_expand_sets <- function(gene_sets) {
    if (inherits(gene_sets, "GeneSet")) gene_sets <- gene_sets@geneSets
    if (!is.list(gene_sets) || is.null(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
        stop("'gene_sets' must be a named list of gene identifier vectors.", call. = FALSE)
    }
    out <- lapply(gene_sets, function(gs) {
        raw <- unlist(gs, use.names = FALSE)
        if (anyDuplicated(as.character(raw))) {
            stop("gene sets must not contain duplicate identifiers: ssGSEA scores a duplicated member twice.", call. = FALSE)
        }
        gs <- unique(as.character(raw))
        gs[!is.na(gs) & nzchar(gs)]
    })
    names(out) <- names(gene_sets)
    if (any(lengths(out) == 0L)) {
        stop("every gene set must contain at least one identifier; ",
             "empty: ", paste0("'", names(out)[lengths(out) == 0L], "'", collapse = ", "),
             call. = FALSE)
    }
    out
}

## ---- expression transform --------------------------------------------------

## `blind = FALSE` is not a default we are free to choose.  DESeq2's
## varianceStabilizingTransformation documentation says blind=FALSE "should be
## used for transforming data for downstream analysis, where the full use of the
## design information should be made".  This route exists to use the design;
## blinding would discard it.  The choice is recorded so a reader can see it.
.dsea_transform <- function(x, transform = c("vst", "rlog", "logcpm"), counts = NULL) {
    transform <- match.arg(transform)
    if (transform == "logcpm") {
        sf <- if (is(x, "DESeqDataSet")) DESeq2::sizeFactors(x) else NULL
        if (is.null(sf)) sf <- .dsea_sizefactors(counts)
        lib <- colSums(counts)
        E <- log2(t(t(counts) / (lib * sf)) + 0.5)
        return(list(E = E, args = list(transform = "logcpm", size_factors = sf)))
    }
    rlang::check_installed("DESeq2", "for variance-stabilising the input.")
    dds <- if (is(x, "DESeqDataSet")) x else DESeq2::DESeqDataSetFromMatrix(counts, DataFrame(counts))
    ## vst()/rlog() need a dispersion trend.  If the caller already ran DESeq()
    ## the object carries one; otherwise vst() estimates it internally.  Either
    ## way, record which happened so the result's provenance is reconstructable.
    disp <- if (is(x, "DESeqDataSet")) mcols(dds)$dispersion else NULL
    fitted <- !is.null(disp) && length(disp) == nrow(dds) && any(is.finite(disp))
    ## vst() is a speed-up that subsets to nsub genes; below that subset size it
    ## refuses outright and hands back advice rather than a matrix.  Fall back to
    ## varianceStabilizingTransformation, which computes every gene.
    use_fast <- nrow(dds) >= 1000L
    args <- list(transform = transform, blind = FALSE, fitted_dds = fitted,
                 engine = if (use_fast) "vst" else "varianceStabilizingTransformation")
    E <- if (transform == "vst") {
        assay(if (use_fast) DESeq2::vst(dds, blind = FALSE)
              else DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))
    } else {
        assay(DESeq2::rlog(dds, blind = FALSE))
    }
    list(E = E, args = args)
}

.dsea_sizefactors <- function(counts) {
    ## median-of-ratios, as DESeq2::estimateSizeFactors does, so that the
    ## fallback path matches the DESeqDataSet path to within rounding.
    lc <- log(counts)
    ok <- is.finite(lc) & lc > 0
    sf <- rep(1, ncol(counts))
    if (all(colSums(counts) > 0)) {
        gm <- exp(rowMeans(lc[ok, , drop = FALSE]))
        ratio <- lc[ok, , drop = FALSE] - log(gm)
        sf <- apply(ratio, 2L, stats::median)
        if (any(!is.finite(sf)) || any(sf <= 0)) sf <- rep(1, ncol(counts))
    }
    sf
}

## ---- gene set filtering ----------------------------------------------------

## Filtering happens on expressed members, BEFORE scoring, so that a set with no
## usable gene produces a dropped record rather than an all-NA score row that
## would later look like a failed test.
.dsea_filter_sets <- function(gene_sets, E, min_size, max_size) {
    expressed <- E[apply(!is.na(E), 1L, any), , drop = FALSE]
    keep <- character(0L)
    sizes <- integer(length(gene_sets))
    for (nm in names(gene_sets)) {
        gs <- intersect(gene_sets[[nm]], rownames(expressed))
        sizes[[nm]] <- length(gs)
        if (length(gs) >= min_size && length(gs) <= max_size) keep <- c(keep, nm)
    }
    list(keep = keep, sizes = sizes)
}

## ---- scoring ---------------------------------------------------------------

.dsea_score <- function(E, gene_sets, method, ...) {
    rlang::check_installed("GSVA", "for pathway scoring.")
    param <- switch(
        method,
        ssgsea = GSVA::ssgseaParam(exprData = E, geneSets = gene_sets,
                                   alpha = 0.25, normalize = TRUE, verbose = FALSE, ...),
        gsva = GSVA::gsvaParam(exprData = E, geneSets = gene_sets,
                               verbose = FALSE, ...)
    )
    S <- GSVA::gsva(param)
    if (is.null(rownames(S))) rownames(S) <- names(gene_sets)
    S
}

#' Pathway scores and pathway-level testing from counts and a design
#'
#' `dsea_scores()` computes a per-sample score for every gene set and returns
#' the score matrix. `dsea()` goes further and tests those scores on an
#' experimental design, returning an effect size, a p-value and an FDR for each
#' pathway. Both start from raw counts: the variance-stabilising transform is
#' fitted inside the call, so a caller does not have to remember to run it, and
#' does not have to get `blind` right.
#'
#' @details
#' \strong{The pipeline.} Counts are variance-stabilised with
#' `DESeq2::vst(blind = FALSE)`, scored per sample by ssGSEA (or GSVA), scaled
#' to z-scores within each pathway, then fitted to the design with ordinary
#' least squares. Residual variances are shrunk towards a common prior by an
#' empirical Bayes step, because with a handful of samples the per-pathway
#' residual variance is otherwise far too noisy to test with.
#'
#' \strong{Why z-score within pathway.} ssGSEA scores are non-linear in the
#' expression values and their scale is set by each pathway's own membership, so
#' scores are not comparable between pathways. Scaling within pathway makes a
#' single linear model able to test all pathways at once, and it is what makes
#' the reported direction meaningful.
#'
#' \strong{What the reported direction means.} `logFC` is the shift of a
#' pathway's score relative to its own mean across samples. It says the score
#' was higher or lower; it is not an absolute statement about pathway activity
#' and it cannot be compared across data sets.
#'
#' \strong{exchangeability.} The test assumes that pathway scores are
#' exchangeable across samples given the design. Sample-label permutation is
#' therefore only valid when treatment is not collinear with batch; `block`
#' should be supplied in that case so that the correlation structure is
#' estimated rather than assumed away.
#'
#' @param x a \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSet}} or a raw
#'   integer counts matrix with gene row names and sample column names.
#' @param gene_sets a named list of gene identifier vectors, identifiers
#'   matching `rownames(x)`. A `GeneSet` object from \pkg{GOSemSim} is also
#'   accepted.
#' @param score_method `"ssgsea"` (default) or `"gsva"`.
#' @param transform `"vst"` (default), `"rlog"` or `"logcpm"`.
#' @param min_size minimum number of expressed members for a set to be tested.
#'   Default 10.
#' @param max_size maximum number of expressed members. Default `Inf`.
#' @param ... passed to the \pkg{GSVA} parameter constructor.
#'
#' @return
#' `dsea_scores()` returns a list with `scores` (pathways by samples matrix,
#' pathway names in rows, sample names in columns) and `meta`.
#'
#' `dsea()` returns `contrasts`, a named list with one entry per tested contrast
#' (`contrast`, `logFC`, `se`, `t`, `df.total`, `p.value`, `adj.p.value` and a
#' `direction` factor taking `"up"`, `"down"` or `"ns"`), plus the shared `scores`
#' matrix, the estimated `correlation`, the fitted `moderation` prior and
#' `meta`. `direction` follows the sign of the effect, not significance;
#' consult `adj.p.value` for that.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' dds <- DESeqDataSetFromMatrix(counts, colData, ~ condition)
#' dds <- DESeq(dds)
#' sets <- list(inflammatory = c("IL6", "TNF", "CXCL8"), stress = c("HSPA1A", "DDIT3"))
#' s <- dsea_scores(dds, sets)
#' design <- model.matrix(~ condition, colData)
#' res <- dsea(dds, sets, design = design, contrast = "conditionB")
#' }
#' @seealso \code{\link{dsea_scores}}
#' @export
dsea_scores <- function(x, gene_sets,
                        score_method = c("ssgsea", "gsva"),
                        transform = c("vst", "rlog", "logcpm"),
                        min_size = 10L, max_size = Inf, ...) {
    score_method <- match.arg(score_method)
    transform <- match.arg(transform)
    counts <- .dsea_counts(x)
    gene_sets <- .dsea_expand_sets(gene_sets)
    tr <- .dsea_transform(x, transform, counts)
    filt <- .dsea_filter_sets(gene_sets, tr$E, min_size, max_size)
    if (length(filt$keep) == 0L) {
        stop("no gene set has between ", min_size, " and ", max_size,
             " expressed members; nothing to score.", call. = FALSE)
    }
    dropped <- setdiff(names(gene_sets), filt$keep)
    if (length(dropped) > 0L) {
        message(sprintf("dsea: %d of %d gene sets dropped by min_size/max_size (%s).",
                        length(dropped), length(gene_sets),
                        paste0("min_size=", min_size, ", max_size=", max_size)))
    }
    S <- .dsea_score(tr$E[rownames(tr$E) %in% unique(unlist(gene_sets[filt$keep], use.names = FALSE)), , drop = FALSE],
                     gene_sets[filt$keep], score_method, ...)
    S <- S[filt$keep, , drop = FALSE]
    na_paths <- rownames(S)[rowSums(!is.na(S)) == 0L]
    if (length(na_paths) > 0L) {
        message(sprintf("dsea: %d pathway(s) scored NA in every sample and are dropped: %s",
                        length(na_paths), paste0("'", head(na_paths, 5L), "'", collapse = ", ")))
        S <- S[setdiff(rownames(S), na_paths), , drop = FALSE]
    }
    if (nrow(S) == 0L) stop("every pathway scored NA in every sample.", call. = FALSE)
    meta <- list(
        transform = transform, transform_args = tr$args,
        score_method = score_method,
        n_input_sets = length(gene_sets), n_tested_sets = nrow(S),
        dropped_sets = length(dropped) + length(na_paths),
        dropped_fraction = (length(dropped) + length(na_paths)) / length(gene_sets),
        set_sizes = filt$sizes[names(S)], min_size = min_size, max_size = max_size,
        genes_used = length(unique(unlist(gene_sets[filt$keep], use.names = FALSE))),
        vst_version = as.character(utils::packageVersion("DESeq2")),
        gsva_version = as.character(utils::packageVersion("GSVA"))
    )
    structure(list(scores = S, meta = meta), class = "DseaScores")
}

#' @rdname dsea_scores
#' @param design an \code{n x p} design matrix, already built with
#'   `model.matrix()`. Must be full rank and must have the same number of rows
#'   as `x` has samples.
#' @param contrast a column name or index of `design`. `NULL` tests every
#'   non-intercept column and returns one result per contrast.
#' @param standardize scale each pathway's scores to mean 0, sd 1 across samples.
#'   Default `TRUE`; see Details.
#' @param moderate `"eb"` (default) shrinks residual variances towards a common
#'   prior; `"none"` uses the raw per-pathway residual variance.
#' @param correlation `"duplicate"` (default) estimates a consensus correlation
#'   between pathway scores and uses it when adjusting; `"none"` assumes
#'   independence.
#' @param block optional grouping factor whose within-block correlation is
#'   preserved when estimating the consensus correlation.
#' @param adjust multiple-testing method passed to `stats::p.adjust`.
#' @export
dsea <- function(x, gene_sets, design, contrast = NULL,
                 score_method = c("ssgsea", "gsva"),
                 transform = c("vst", "rlog"),
                 standardize = TRUE,
                 moderate = c("eb", "none"),
                 correlation = c("duplicate", "none"),
                 block = NULL, adjust = "BH", ...) {
    score_method <- match.arg(score_method)
    transform <- match.arg(transform)
    moderate <- match.arg(moderate)
    correlation <- match.arg(correlation)
    design <- .dsea_check_design(design, colnames(x))
    scores <- dsea_scores(x, gene_sets, score_method = score_method,
                          transform = transform, ...)
    Z <- if (standardize) .dsea_standardize(scores$scores) else scores$scores
    dropped <- attr(Z, "dropped_constant")
    if (length(dropped) > 0L) {
        message(sprintf("dsea: %d pathway(s) have zero variance across samples and are dropped: %s",
                        length(dropped), paste0("'", head(dropped, 5L), "'", collapse = ", ")))
        Z <- Z[setdiff(rownames(Z), dropped), , drop = FALSE]
    }
    if (nrow(Z) == 0L) stop("no pathway has non-zero variance across samples.", call. = FALSE)
    idx <- .dsea_resolve_contrasts(contrast, design)
    rho <- if (correlation == "duplicate") {
        .dsea_consensus_correlation(Z, design, block)
    } else 0
    if (correlation == "duplicate" && rho < 0.1) {
        message(sprintf("dsea: consensus correlation is %.3f; correlation=\"none\" is nearly equivalent here.",
                        rho))
    }
    prior <- NULL
    out <- lapply(idx, function(ci) {
        fit <- .dsea_fit(Z, design, ci, moderate)
        prior <<- fit$moderation
        p <- 2 * stats::pt(-abs(fit$t), fit$df.total)
        adj <- stats::p.adjust(p, method = adjust)
        direction <- rep("ns", length(p))
        direction[is.finite(fit$logFC) & fit$logFC > 0] <- "up"
        direction[is.finite(fit$logFC) & fit$logFC < 0] <- "down"
        direction <- stats::setNames(factor(direction, levels = c("up", "down", "ns")),
                                     rownames(Z))
        list(contrast = colnames(design)[ci],
             logFC = stats::setNames(fit$logFC, rownames(Z)),
             se = stats::setNames(fit$se, rownames(Z)),
             t = stats::setNames(fit$t, rownames(Z)),
             df.total = fit$df.total,
             p.value = stats::setNames(p, rownames(Z)),
             adj.p.value = stats::setNames(adj, rownames(Z)),
             direction = factor(direction, levels = c("up", "down", "ns")))
    })
    names(out) <- vapply(out, `[[`, character(1), "contrast")
    .dsea_guard_scale(out)
    structure(list(contrasts = out, scores = scores$scores,
                   correlation = list(consensus = rho, method = correlation, block = block),
                   moderation = prior,
                   meta = c(scores$meta, list(design = colnames(design), contrast = NULL,
                                              moderate = moderate, adjust = adjust,
                                              standardize = standardize))),
              class = "DseaResult")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.DseaScores <- function(x, ...) {
    cat("<DseaScores>", nrow(x$scores), "pathways x", ncol(x$scores), "samples\n")
    cat("  transform:", x$meta$transform, "| scoring:", x$meta$score_method, "\n")
    cat("  tested:", x$meta$n_tested_sets, "of", x$meta$n_input_sets,
        sprintf("sets (dropped %.1f%%)", 100 * x$meta$dropped_fraction), "\n")
    invisible(x)
}

#' @export
print.DseaResult <- function(x, ...) {
    cat("<DseaResult>", nrow(x$scores), "pathways x", ncol(x$scores), "samples,",
        length(x$contrasts), "contrast(s)\n")
    for (r in x$contrasts) {
        cat(sprintf("  [%s] %d at FDR < 0.05 (%d up, %d down)\n", r$contrast,
                    sum(r$adj.p.value < 0.05, na.rm = TRUE),
                    sum(r$direction == "up" & r$adj.p.value < 0.05, na.rm = TRUE),
                    sum(r$direction == "down" & r$adj.p.value < 0.05, na.rm = TRUE)))
    }
    cat(sprintf("  consensus correlation: %.4f | moderation: %s\n",
                x$correlation$consensus, x$moderation$method))
    invisible(x)
}