## dsea() -- design-aware GSEA, and a self-contained ssGSEA scorer.
##
## WHAT THIS IS.  Counts plus a design, in one call: variance-stabilise, reduce
## the design to a per-gene contrast statistic, rank, and hand the ranking to
## the GSEA engine already in this package.  dsea_scores() is the other half of
## the product: a per-sample score matrix, for analyses that want samples rather
## than one contrast.
##
## WHY GSEA AND NOT A SET TEST OF OUR OWN.
## We built a negative-binomial set test, then an empirical-Bayes moderated one,
## and benchmarked them against limma's fry and against plain GSEA on simulated
## KEGG data with known truth.  Across eight signal structures the ranking-based
## route reached recall 0.93-1.00 and AUC 0.978-0.993, while the moderated set
## test reached 0.00-0.71 and 0.864-0.972.  Our test did beat fry decisively
## (AUC 0.86-0.97 against 0.28-0.37), but fry is not what a user reaches for,
## and a route that loses to the obvious baseline has no reason to exist.  The
## set test was therefore removed rather than tuned.
##
## The ranking statistic was tested too -- plain mean difference, moderated t,
## DESeq2 Wald -- in the project's simulation benchmark (64 replicates per
## scenario, 213 KEGG sets, planted effects with known truth).  Under uniform or
## heterogeneous effects the four are indistinguishable; under an eightfold
## dispersion multiplier the mean difference recalls 0.979 of the planted sets
## and the Wald statistic 0.826.  Variance weighting on top of a rank-based
## statistic double-counts, because GSEA has already equalised the genes by rank,
## and it is worst exactly where the affected genes are the noisiest ones.

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
    ## Tolerate one-off floating point noise on counts meant to be integers, but
    ## refuse a genuinely fractional matrix: rounding it would hide a mis-read
    ## input file rather than report it.
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

## ssGSEA weights a duplicated member twice and would quietly inflate that
## set's contribution, so duplicates are an error rather than a silent fix.
.dsea_expand_sets <- function(gene_sets) {
    if (inherits(gene_sets, "GeneSet")) gene_sets <- gene_sets@geneSets
    if (!is.list(gene_sets) || is.null(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
        stop("'gene_sets' must be a named list of gene identifier vectors.", call. = FALSE)
    }
    out <- lapply(gene_sets, function(gs) {
        raw <- unlist(gs, use.names = FALSE)
        if (anyDuplicated(as.character(raw))) {
            stop("gene sets must not contain duplicate identifiers: ssGSEA scores ",
                 "a duplicated member twice.", call. = FALSE)
        }
        gs <- as.character(raw)
        gs[!is.na(gs) & nzchar(gs)]
    })
    names(out) <- names(gene_sets)
    if (any(lengths(out) == 0L)) {
        stop("every gene set must contain at least one identifier; empty: ",
             paste0("'", names(out)[lengths(out) == 0L], "'", collapse = ", "), call. = FALSE)
    }
    out
}

## ---- expression transform --------------------------------------------------

## blind = FALSE is not a free choice.  DESeq2's own documentation says it
## "should be used for transforming data for downstream analysis, where the full
## use of the design information should be made"; this function exists to use the
## design, so blinding it would discard the thing being offered.
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
    if (is(x, "DESeqDataSet")) {
        dds <- x
    } else {
        ## A bare counts matrix has no colData and vst() needs one row per
        ## sample.  There is no design on this path, which is exactly why
        ## blind = FALSE is still right: nothing is blinded because nothing was
        ## known.
        cds <- S4Vectors::DataFrame(sample = seq_along(colnames(counts)),
                                    row.names = colnames(counts))
        dds <- DESeq2::DESeqDataSetFromMatrix(counts, cds, ~ 1)
    }
    ## DESeq2::dispersions() is the documented accessor; mcols() would need
    ## SummarizedExperiment in Imports and was only resolving in tests because
    ## the test file had attached it.
    disp <- if (is(x, "DESeqDataSet")) DESeq2::dispersions(dds) else NULL
    fitted <- !is.null(disp) && length(disp) == nrow(dds) && any(is.finite(disp))
    ## vst() subsets to nsub genes as a speed-up and refuses outright below it.
    use_fast <- nrow(dds) >= 1000L
    args <- list(transform = transform, blind = FALSE, fitted_dds = fitted,
                 engine = if (use_fast) "vst" else "varianceStabilizingTransformation")
    ## SummarizedExperiment::assay, qualified.  Unqualified it resolves only when
    ## the caller happens to have attached SummarizedExperiment, which the test
    ## file does and the package namespace does not -- so it passes tests and
    ## fails for a user.
    E <- if (transform == "vst") {
        SummarizedExperiment::assay(
            if (use_fast) DESeq2::vst(dds, blind = FALSE)
            else DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))
    } else {
        SummarizedExperiment::assay(DESeq2::rlog(dds, blind = FALSE))
    }
    list(E = E, args = args)
}

## Median-of-ratios, matching DESeq2::estimateSizeFactors on a bare matrix.
.dsea_sizefactors <- function(counts) {
    lc <- log(counts)
    sf <- rep(1, ncol(counts))
    if (all(colSums(counts) > 0)) {
        ## ok flags CELLS that may enter the geometric mean; it has to be
        ## applied to rows so lc[ok, ] and gm share a row set.  Filtering
        ## columns here let the subscript outlive the matrix and the median
        ## fell off the end.
        ok <- apply(lc, 1L, function(v) all(is.finite(v)) && all(v > 0))
        if (any(ok)) {
            lcs <- lc[ok, , drop = FALSE]
            sf <- apply(lcs - log(exp(rowMeans(lcs))), 2L, stats::median)
        }
        if (any(!is.finite(sf)) || any(sf <= 0)) sf <- rep(1, ncol(counts))
    }
    sf
}

## ---- design checks ---------------------------------------------------------

.dsea_check_design <- function(design, samples) {
    if (is.null(design)) stop("'design' is required.", call. = FALSE)
    design <- as.matrix(design)
    if (!is.numeric(design)) {
        stop("'design' must be a numeric matrix. Build it with model.matrix(), ",
             "e.g. model.matrix(~ condition, data = colData); got ", class(design)[1], ".",
             call. = FALSE)
    }
    if (any(!is.finite(design))) stop("'design' contains non-finite values.", call. = FALSE)
    if (nrow(design) != length(samples)) {
        stop("'design' has ", nrow(design), " rows but 'x' has ", length(samples),
             " samples.", call. = FALSE)
    }
    if (is.null(colnames(design))) {
        colnames(design) <- paste0("x", seq_len(ncol(design)))
        warning("'design' had no column names; refer to columns by index only.",
                call. = FALSE)
    }
    if (!is.null(rownames(design)) && !identical(rownames(design), samples)) {
        stop("rownames(design) do not match colnames(x); reorder the design to match the samples.",
             call. = FALSE)
    }
    rank <- qr(design)$rank
    if (rank < ncol(design)) {
        ## Checked before any fitting and named explicitly: a rank-deficient
        ## design is the airway dex/cell trap, and it must not surface later as
        ## an unreadable numerical error.
        stop("'design' is not full rank (rank ", rank, " of ", ncol(design),
             " columns). One or more design terms are linear combinations of ",
             "the others and cannot be estimated; drop the redundant terms first.",
             call. = FALSE)
    }
    design
}

.dsea_resolve_one <- function(contrast, design) {
    if (is.character(contrast)) {
        i <- match(contrast, colnames(design))
        if (is.na(i)) {
            stop("contrast '", contrast, "' not found in design; available: ",
                 paste(colnames(design), collapse = ", "), call. = FALSE)
        }
        return(i)
    }
    i <- as.integer(contrast)
    if (length(i) != 1L || is.na(i) || i < 1L || i > ncol(design)) {
        stop("'contrast' must be a single column name or index in 1..", ncol(design),
             call. = FALSE)
    }
    i
}

## A design column is only a CONTRAST if it encodes a comparison.  Two kinds of
## column do not, and both produce a plausible-looking table while answering a
## different question:
##
##   - a constant column (the intercept).  Ranking by it orders genes by overall
##     mean expression, so GSEA reports highly expressed pathways rather than
##     differentially enriched ones.  This was the DEFAULT (contrast = 1L) until
##     the guard below existed, and nothing errored: on a two-group design the
##     recorded contrast was simply "(Intercept)".
##   - an indicator column in a design with no constant column anywhere, i.e.
##     cell-means coding such as `~ 0 + condition`.  Its coefficient is that
##     group's mean, so ranking by it is again an expression-level ranking.
##
## A continuous column is a real coefficient (a dose slope, say) and stays
## allowed, so the guard is deliberately about constants and group indicators
## rather than about "the design must have an intercept".
.dsea_check_contrast <- function(design, contrast) {
    idx <- .dsea_resolve_one(contrast, design)
    col <- design[, idx]
    if (all(col == col[1L])) {
        others <- setdiff(colnames(design), colnames(design)[idx])
        stop("design column '", colnames(design)[idx], "' is constant (an intercept); ",
             "ranking by it orders genes by mean expression, not by the contrast. ",
             if (length(others)) paste0("Choose the column that encodes the comparison, e.g. '",
                                        others[1L], "'. ") else "",
             call. = FALSE)
    }
    has_constant <- any(vapply(seq_len(ncol(design)),
                               function(j) all(design[, j] == design[1L, j]), logical(1)))
    if (!has_constant && all(col %in% c(0, 1))) {
        stop("'design' has no intercept, so column '", colnames(design)[idx],
             "' is a group indicator and its coefficient is that group's mean, not a ",
             "contrast. Build the design with an intercept (e.g. ~ condition) and select ",
             "the condition column, or supply an explicit contrast vector.", call. = FALSE)
    }
    idx
}

## The transform step reads the design stored in a DESeqDataSet, while the
## ranking step reads the `design` argument.  When those disagree the analysis
## silently mixes two models: measured on a dds fitted with ~ cond and a
## ~ patient + cond matrix passed in, the transformed values differ by up to
## 0.648 log2 units.  There is no way to honour a design matrix without a
## formula, so the honest fix is to say so instead of blending them quietly.
##
## The comparison is on COLUMN SPACE, not on entries.  `vst()` only depends on
## the span of the design, and the same model is routinely written two ways:
## model.matrix(~ cond) codes a dummy for the non-reference level while a
## hand-built cbind(1, as.numeric(cond)) codes the other level, and the two span
## the same space.  Comparing entries would warn on a perfectly consistent call.
.dsea_check_design_agreement <- function(x, design) {
    if (!is(x, "DESeqDataSet")) return(invisible(TRUE))
    X <- tryCatch({
        cd <- as.data.frame(SummarizedExperiment::colData(x))
        stats::model.matrix(DESeq2::design(x), cd)
    }, error = function(e) NULL)
    if (is.null(X) || nrow(X) != nrow(design)) return(invisible(TRUE))
    ## is `design` inside the span of X?  normal equations on a full-rank X
    fitted <- tryCatch({
        X %*% solve(crossprod(X), crossprod(X, design))
    }, error = function(e) NULL)
    if (is.null(fitted)) return(invisible(TRUE))
    same <- max(abs(fitted - design)) < 1e-8
    if (!same) {
        warning("'design' differs from design(x): the variance-stabilising transform uses ",
                "design(x) while the ranking uses 'design'. Refit the DESeqDataSet with the ",
                "intended design (design(dds) <- ...; DESeq(dds)) or pass a counts matrix ",
                "so the transform can use 'design'.", call. = FALSE)
    }
    invisible(same)
}

## ---- ssGSEA ----------------------------------------------------------------

#' Single-sample GSEA scores
#'
#' Score each sample independently against every gene set, giving a gene set by
#' sample matrix that can be clustered, plotted or tested.
#'
#' @details
#' For sample \eqn{j} the genes are ranked by expression and the running
#' enrichment score of a set is accumulated over that ranking, with the same
#' weighting as [gseaScores()]. A set's score is the integral of its running
#' score over the ranking, normalised by the number of ranked genes. That is the
#' ssGSEA definition, and it is computed here rather than imported so that the
#' package needs no additional dependency for it.
#'
#' Scores are comparable ACROSS SAMPLES but not across gene sets, because the
#' scale is set by each set's own membership. Scale within gene set before any
#' cross-set comparison.
#'
#' The caller supplies already-transformed values; `dsea_scores()` does that step
#' and records what it did.
#'
#' @param E a genes by samples matrix of expression values, normalised and on a
#'   log scale. Row names are gene identifiers, column names samples.
#' @param gene_sets a named list of gene identifier vectors.
#' @param exponent weight applied to the rank statistic; 1 is the GSEA default,
#'   0 gives the classic unweighted Kolmogorov-Smirnov form.
#' @param min_size,max_size gene sets outside this range, after intersecting
#'   with `rownames(E)`, are dropped.
#' @return A gene set by sample matrix of scores.
#' @seealso \code{\link{dsea_scores}}, \code{\link{gseaScores}}
#' @export
ssgsea_scores <- function(E, gene_sets, exponent = 1, min_size = 10L, max_size = 500L) {
    E <- as.matrix(E)
    if (is.null(rownames(E)) || is.null(colnames(E))) {
        stop("'E' must have gene row names and sample column names.", call. = FALSE)
    }
    gene_sets <- filter_gene_sets(.dsea_expand_sets(gene_sets), rownames(E))
    keep <- lengths(gene_sets) >= min_size & lengths(gene_sets) <= max_size
    gene_sets <- gene_sets[keep]
    if (length(gene_sets) == 0L) {
        stop("no gene set has between ", min_size, " and ", max_size,
             " genes in common with rownames(E).", call. = FALSE)
    }
    n <- ncol(E)
    out <- matrix(NA_real_, length(gene_sets), n,
                  dimnames = list(names(gene_sets), colnames(E)))
    for (j in seq_len(n)) {
        v <- E[, j]
        ok <- is.finite(v)
        if (!any(ok)) next
        v <- sort(v[ok], decreasing = TRUE)
        stat <- stats::setNames(abs(v), names(v))
        for (i in seq_along(gene_sets)) {
            in_set <- names(stat) %in% gene_sets[[i]]
            ## A set with no member ranked in this sample has a running score of
            ## zero everywhere and an integral of zero.  That is a measurement
            ## ("nothing here"), not a missing value, so it stays 0 rather than
            ## NA and does not silently drop out of later averaging.
            if (!any(in_set)) { out[i, j] <- 0; next }
            r <- gsea_scores_cpp(stat, in_set, exponent)$runningScore
            out[i, j] <- sum(r) / length(r)
        }
    }
    out
}

## ---- design-driven GSEA ----------------------------------------------------

#' Differential set enrichment from a design
#'
#' `dsea()` takes counts and a design and returns GSEA results for a contrast:
#' it transforms the counts, reduces the design to a per-gene statistic, ranks
#' the genes and runs the GSEA engine. `dsea_scores()` stops one step earlier and
#' returns the gene set by sample score matrix instead.
#'
#' @details
#' \strong{Why the ranking is what it is.} Genes are ranked by the difference in
#' mean transformed expression between the two sides of the contrast. That looks
#' like the weakest available choice and it is the one the benchmark supported:
#' ranking by a moderated t or by DESeq2's Wald statistic was no better under
#' uniform or heterogeneous effects, and materially worse when the affected genes
#' were the noisier ones (recall 0.826 for Wald against 0.979 for the mean
#' difference under an eightfold dispersion multiplier; across a 1x/2x/4x series
#' the Wald statistic falls 0.995 to 0.950 while the mean difference stays flat).
#' GSEA's statistic already equalises genes by rank, so layering a variance
#' weight on top discards real signal at the extremes rather than adding power.
#'
#' \strong{Which column may be the contrast.} A design column encodes a
#' comparison only when it is not constant and not a group indicator under
#' cell-means coding. The intercept is rejected because ranking by it orders
#' genes by mean expression, which turns the analysis into "highly expressed
#' pathways"; a column of `~ 0 + condition` is rejected for the same reason,
#' since its coefficient is one group's mean rather than a difference. A
#' continuous column (a dose slope, say) is allowed.
#'
#' \strong{The transform uses design(x), the ranking uses \code{design}.} When
#' `x` is a `DESeqDataSet`, `vst(blind = FALSE)` reads the design stored in the
#' object, while the ranking reads the `design` argument. If the two differ the
#' call warns rather than blending two models: refit the object with the intended
#' design, or pass a counts matrix so the transform can use `design`.
#'
#' \strong{What the design is used for.} For a two-group contrast the ranking is
#' the group mean difference. For a design with nuisance terms it is the contrast
#' coefficient from the fitted linear model, so batch and other covariates are
#' removed before ranking rather than ignored.
#'
#' \strong{What this is not.} It is a design-aware front end to GSEA, not a new
#' set test. Directions are relative to the ranking: a set scoring high sits
#' toward the top, which is the direction the contrast was oriented in.
#'
#' @param x a counts matrix, or a `DESeqDataSet`.
#' @param gene_sets a named list of gene identifier vectors, identifiers matching
#'   `rownames(x)`.
#' @param design for `dsea()`: an `n x p` design matrix from `model.matrix()`,
#'   row-aligned with the samples. Must be full rank.
#' @param contrast a column name or index of `design`. Defaults to the only
#'   non-intercept column of a two-column design; a wider design must name the
#'   contrast explicitly, because the useful column cannot be guessed.
#' @param transform `"vst"` (default) or `"rlog"`; `dsea_scores()` also accepts
#'   `"logcpm"`.
#' @param min_size,max_size gene set size limits, applied to the overlap with the
#'   genes actually present.
#' @param exponent for `dsea_scores()`: rank weighting passed to
#'   [ssgsea_scores()].
#' @param ... for `dsea()`, passed to [gsea()].
#'
#' @return
#' `dsea()` returns the data.frame from [gsea()], sorted by p-value, with two
#' columns added: `padj`, the Benjamini-Hochberg adjusted p-value, and
#' `contrast`, the design column that was tested. Everything else is GSEA's own
#' output, so the table can be filtered and plotted like any other GSEA result.
#'
#' `dsea_scores()` returns a list with `scores` (gene sets by samples) and `meta`
#' recording the transform and scorer settings. `meta$set_sizes` reports the
#' number of members actually scored, after intersecting with the data, and
#' `meta$raw_set_sizes` the input size, so the two can be compared.
#'
#' @examples
#' \dontrun{
#' dds <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(counts, coldata, ~ condition))
#' design <- model.matrix(~ condition, coldata)
#' res <- dsea(dds, gene_sets, design = design, contrast = "conditionB")
#' head(res[order(res$pvalue), ])
#'
#' sc <- dsea_scores(dds, gene_sets)
#' dim(sc$scores)
#' }
#' @seealso \code{\link{gsea}}, \code{\link{ssgsea_scores}}
#' @export
dsea <- function(x, gene_sets, design, contrast = NULL,
                 transform = c("vst", "rlog"),
                 min_size = 10L, max_size = 500L, ...) {
    transform <- match.arg(transform)
    counts <- .dsea_counts(x)
    design <- .dsea_check_design(design, colnames(counts))
    ## The default is only safe in the two-column case, and even there it has to
    ## be the non-intercept column.  Anything wider requires the caller to say
    ## which coefficient is the contrast; guessing produced expression-level
    ## rankings (see .dsea_check_contrast).
    if (is.null(contrast)) {
        if (ncol(design) != 2L) {
            stop("'contrast' must be supplied when 'design' has ", ncol(design),
                 " columns. Name or index the column that encodes the comparison; available: ",
                 paste(colnames(design), collapse = ", "), call. = FALSE)
        }
        contrast <- 2L
    }
    idx <- .dsea_check_contrast(design, contrast)
    .dsea_check_design_agreement(x, design)
    ## The same gene-set validation dsea_scores() applies, so both entry points
    ## fail on the same inputs rather than one of them failing later or not at all.
    gene_sets <- .dsea_expand_sets(gene_sets)
    E <- .dsea_transform(x, transform, counts)$E

    ## Rank by the contrast rather than by a marginal group difference: with
    ## nuisance terms the coefficient is what removes them, and for a plain
    ## two-group design the two agree exactly.
    rank <- .dsea_rank_by_contrast(E, design, idx)
    rank <- sort(rank, decreasing = TRUE)

    res <- gsea(geneList = rank, gene_sets = gene_sets,
                minGSSize = min_size, maxGSSize = max_size, ...)
    if (is.null(res) || nrow(res) == 0L) return(res)
    ## gsea() reports a raw p per set and no multiplicity column.  A caller who
    ## has to remember p.adjust() after every call will eventually forget, so the
    ## adjusted value is computed here.  The raw column is kept untouched.
    res$padj <- stats::p.adjust(res$pvalue, method = "BH")
    res$contrast <- colnames(design)[idx]
    res[order(res$pvalue), , drop = FALSE]
}

## Per-gene contrast statistic. A two-group contrast reduces to the group mean
## difference; anything else needs the fitted coefficient, so the general path
## runs OLS on the transformed matrix and reads off the contrast row.
.dsea_rank_by_contrast <- function(E, design, idx) {
    if (ncol(design) == 2L && idx == 2L && all(design[, 1L] == 1)) {
        g <- design[, 2L]
        u <- sort(unique(g))
        if (length(u) == 2L) {
            return(rowMeans(E[, g == u[2L], drop = FALSE]) -
                   rowMeans(E[, g == u[1L], drop = FALSE]))
        }
    }
    ## Normal equations, not a hand-assembled QR.  qr.Q()/qr.R() do not return
    ## matching factors across R's two QR backends -- with LAPACK = TRUE the
    ## reported Q and R multiply back to something that is not X once the design
    ## has three or more columns -- and the rebuilt coefficient silently differs
    ## from the OLS one.  The design is at most a handful of columns wide here,
    ## so forming X'X explicitly costs nothing and cannot be mis-assembled.
    Y <- t(E)
    XtX <- crossprod(design)
    coef <- solve(XtX, crossprod(design, Y))
    stats::setNames(coef[idx, ], rownames(E))
}

## ---- pathway score matrix --------------------------------------------------

#' @rdname dsea
#' @export
dsea_scores <- function(x, gene_sets, transform = c("vst", "rlog", "logcpm"),
                        min_size = 10L, max_size = 500L, exponent = 1, ...) {
    transform <- match.arg(transform)
    counts <- .dsea_counts(x)
    gene_sets <- .dsea_expand_sets(gene_sets)
    tr <- .dsea_transform(x, transform, counts)
    S <- ssgsea_scores(tr$E, gene_sets, exponent = exponent,
                       min_size = min_size, max_size = max_size)
    dropped <- setdiff(names(gene_sets), rownames(S))
    if (length(dropped) > 0L) {
        message(sprintf("dsea_scores: %d of %d gene sets dropped (size outside [%d, %d] after intersecting with the data).",
                        length(dropped), length(gene_sets), min_size, max_size))
    }
    ## set_sizes must report what the size filter and the scorer ACTUALLY used,
    ## which is the overlap with the data, not the raw input size.  A set of 60
    ## identifiers of which 20 are present scores 20 genes; recording 60 makes
    ## the metadata unable to explain the result it accompanies.  The raw size is
    ## kept alongside so the gap stays visible.
    kept <- rownames(S)
    meta <- list(
        transform = transform, transform_args = tr$args,
        scorer = "ssgsea", exponent = exponent,
        n_input_sets = length(gene_sets), n_scored_sets = nrow(S),
        dropped_sets = length(dropped),
        dropped_fraction = length(dropped) / length(gene_sets),
        set_sizes = lengths(filter_gene_sets(gene_sets[kept], rownames(tr$E))),
        raw_set_sizes = lengths(gene_sets[kept]),
        min_size = min_size, max_size = max_size,
        vst_version = as.character(utils::packageVersion("DESeq2"))
    )
    structure(list(scores = S, meta = meta), class = "DseaScores")
}

#' @export
print.DseaScores <- function(x, ...) {
    cat("<DseaScores>", nrow(x$scores), "gene sets x", ncol(x$scores), "samples\n")
    cat("  transform:", x$meta$transform, "| scorer:", x$meta$scorer, "\n")
    cat(sprintf("  scored %d of %d sets (dropped %.1f%%)\n", x$meta$n_scored_sets,
                x$meta$n_input_sets, 100 * x$meta$dropped_fraction))
    invisible(x)
}
