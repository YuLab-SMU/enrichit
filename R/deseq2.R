## Adapters: extract enrichment-analysis inputs from DESeq2 objects.
##
## `as_ranked()` produces a named, descending-sorted statistic vector for
## gsea(); `as_genes()` produces a significant-gene character vector for
## ora(). Both work directly on DESeqResults, or on a fitted DESeqDataSet
## by forwarding `...` to DESeq2::results(). DESeq2 itself is only needed
## for the DESeqDataSet methods and is checked lazily.

#' Convert DESeq2 results to a ranked gene list for GSEA
#'
#' Extract a per-gene statistic from a \pkg{DESeq2} analysis as a named
#' numeric vector sorted in descending order, ready to be passed to
#' \code{\link{gsea}}. This uses the full ranked list rather than a
#' significance cut-off, so no information is discarded by dichotomising
#' the genes.
#'
#' @details Non-finite statistics (genes with all-zero counts, genes
#'   flagged by outlier replacement, or failed fits) are dropped with a
#'   message. Duplicated gene identifiers are resolved by keeping the
#'   entry with the largest absolute statistic, with a warning. Note
#'   that \code{lfcShrink()} with \code{type = "apeglm"} does not return
#'   a \code{stat} column; use \code{statistic = "log2FoldChange"} for
#'   shrunk results, and keep in mind that shrunk fold changes are
#'   intended for ranking and visualization, not for statistical
#'   testing.
#'
#' @param x a \code{\link[DESeq2:DESeqDataSet-class]{DESeqResults}} or
#'   \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSet}} object. For a
#'   \code{DESeqDataSet}, \code{...} is forwarded to
#'   \code{\link[DESeq2:results]{DESeq2::results()}} (e.g. \code{contrast},
#'   \code{name}, \code{alpha}, \code{lfcThreshold}).
#' @param statistic character, the column to rank by: \code{"stat"} (the
#'   Wald statistic, default) or \code{"log2FoldChange"}. The Wald
#'   statistic is generally preferred for ranking as it incorporates the
#'   per-gene standard error.
#' @param ... additional arguments; for a \code{DESeqDataSet} these are
#'   passed to \code{DESeq2::results()}.
#' @return A named numeric vector sorted in descending order, suitable
#'   for the \code{geneList} argument of \code{\link{gsea}}.
#' @examples
#' \donttest{
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(n = 500)
#' ## strengthen a subset of genes in condition B (keep the rest near-null
#' ## so the size factors have a stable reference)
#' isB <- colData(dds)$condition == "B"
#' ct <- counts(dds)
#' ct[1:60, isB] <- ct[1:60, isB] * 5 + 500
#' storage.mode(ct) <- "integer"
#' counts(dds) <- ct
#' dds <- DESeq(dds, quiet = TRUE)
#' ranked <- as_ranked(dds, contrast = c("condition", "B", "A"))
#' gene_sets <- list(
#'     SET1 = paste0("gene", 1:50),
#'     SET2 = paste0("gene", 201:260)
#' )
#' res <- gsea(ranked, gene_sets, nPerm = 1000)
#' }
#' @export
as_ranked <- function(x, statistic = c("stat", "log2FoldChange"), ...) {
    UseMethod("as_ranked")
}

#' @rdname as_ranked
#' @export
as_ranked.DESeqResults <- function(x, statistic = c("stat", "log2FoldChange"), ...) {
    statistic <- match.arg(statistic)
    .as_ranked_from_results(x, statistic)
}

#' @rdname as_ranked
#' @export
as_ranked.DESeqDataSet <- function(x, statistic = c("stat", "log2FoldChange"), ...) {
    rlang::check_installed("DESeq2", "for converting DESeq2 objects.")
    statistic <- match.arg(statistic)
    res <- DESeq2::results(x, ...)
    as_ranked(res, statistic = statistic)
}

.get_results_column <- function(x, col) {
    ## `[[` works on data.frame / DataFrame (real DESeqResults) but fails
    ## on classed matrices, so fall back to `[` when needed
    out <- tryCatch(x[[col]], error = function(e) NULL)
    if (is.null(out) && !is.null(colnames(x)) && col %in% colnames(x)) {
        out <- x[, col]
    }
    if (!is.null(out)) {
        out <- as.numeric(out)
    }
    out
}

.as_ranked_from_results <- function(x, statistic) {
    stat_col <- switch(statistic,
        stat = "stat",
        log2FoldChange = "log2FoldChange"
    )
    stats <- .get_results_column(x, stat_col)
    if (is.null(stats)) {
        if (statistic == "stat") {
            stop(
                "no 'stat' column found in the DESeqResults object; ",
                "results from lfcShrink(type = 'apeglm') lack it -- ",
                "use statistic = 'log2FoldChange' instead.",
                call. = FALSE
            )
        }
        stop("no '", stat_col, "' column found in the DESeqResults object.", call. = FALSE)
    }

    ids <- rownames(x)
    if (is.null(ids)) {
        stop("the DESeqResults object has no row names; gene identifiers are required.", call. = FALSE)
    }
    names(stats) <- ids

    bad <- !is.finite(stats)
    if (any(bad)) {
        message(
            "dropping ", sum(bad), " gene(s) with non-finite ", statistic,
            " (all-zero counts, Cook's-distance outlier replacement, or failed fits)."
        )
        stats <- stats[!bad]
    }
    if (length(stats) == 0) {
        stop("no finite ", statistic, " values available.", call. = FALSE)
    }

    if (anyDuplicated(names(stats)) > 0) {
        n_dup <- sum(duplicated(names(stats)))
        warning(
            n_dup, " duplicated gene identifier(s) found; for each id ",
            "keeping the entry with the largest absolute ", statistic, ".",
            call. = FALSE
        )
        stats <- stats[order(abs(stats), decreasing = TRUE)]
        stats <- stats[!duplicated(names(stats))]
    }

    sort(stats, decreasing = TRUE)
}

#' Extract significant genes from DESeq2 results for ORA
#'
#' Select the genes passing an adjusted-p-value cut-off (optionally
#' split by direction of change) from a \pkg{DESeq2} analysis, as a
#' character vector of gene identifiers for \code{\link{ora}}.
#'
#' @details Selection uses \code{padj < cutoff} on the \code{padj}
#'   column as produced by \code{\link[DESeq2:results]{DESeq2::results()}}
#'   (which includes independent filtering). \code{minAbsLFC} is a
#'   post-hoc filter applied to the reported \code{log2FoldChange}; it
#'   does not change the underlying test. To test against a minimum
#'   effect size, set \code{lfcThreshold} in \code{results()} /
#'   \code{DESeq2::results()} instead.
#'
#' @param x a \code{\link[DESeq2:DESeqDataSet-class]{DESeqResults}} or
#'   \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSet}} object. For a
#'   \code{DESeqDataSet}, \code{...} is forwarded to
#'   \code{DESeq2::results()}.
#' @param cutoff numeric, adjusted p-value cut-off (default 0.05); genes
#'   with \code{padj < cutoff} are kept.
#' @param direction character, \code{"both"} (default), \code{"up"} or
#'   \code{"down"}; selects genes with positive or negative
#'   \code{log2FoldChange} among the significant ones.
#' @param minAbsLFC numeric, keep only significant genes with
#'   \code{abs(log2FoldChange) > minAbsLFC} (default 0, no filtering).
#' @param ... additional arguments; for a \code{DESeqDataSet} these are
#'   passed to \code{DESeq2::results()}.
#' @return A character vector of gene identifiers (empty if no gene
#'   passes, with a warning).
#' @examples
#' \donttest{
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(n = 500)
#' isB <- colData(dds)$condition == "B"
#' ct <- counts(dds)
#' ct[1:60, isB] <- ct[1:60, isB] * 5 + 500
#' storage.mode(ct) <- "integer"
#' counts(dds) <- ct
#' dds <- DESeq(dds, quiet = TRUE)
#' up_genes <- as_genes(dds, contrast = c("condition", "B", "A"), direction = "up")
#' gene_sets <- list(
#'     SET1 = paste0("gene", 1:50),
#'     SET2 = paste0("gene", 201:260)
#' )
#' res <- ora(gene = up_genes, gene_sets = gene_sets, universe = rownames(dds))
#' }
#' @export
as_genes <- function(x,
                     cutoff = 0.05,
                     direction = c("both", "up", "down"),
                     minAbsLFC = 0,
                     ...) {
    UseMethod("as_genes")
}

#' @rdname as_genes
#' @export
as_genes.DESeqResults <- function(x,
                                  cutoff = 0.05,
                                  direction = c("both", "up", "down"),
                                  minAbsLFC = 0,
                                  ...) {
    direction <- match.arg(direction)
    .extract_sig_genes(x, cutoff = cutoff, direction = direction, minAbsLFC = minAbsLFC)
}

#' @rdname as_genes
#' @export
as_genes.DESeqDataSet <- function(x,
                                  cutoff = 0.05,
                                  direction = c("both", "up", "down"),
                                  minAbsLFC = 0,
                                  ...) {
    rlang::check_installed("DESeq2", "for converting DESeq2 objects.")
    direction <- match.arg(direction)
    res <- DESeq2::results(x, ...)
    as_genes(res, cutoff = cutoff, direction = direction, minAbsLFC = minAbsLFC)
}

.extract_sig_genes <- function(x, cutoff, direction, minAbsLFC) {
    padj <- .get_results_column(x, "padj")
    if (is.null(padj)) {
        stop("no 'padj' column found in the DESeqResults object.", call. = FALSE)
    }
    ids <- rownames(x)
    if (is.null(ids)) {
        stop("the DESeqResults object has no row names; gene identifiers are required.", call. = FALSE)
    }

    if (any(is.na(padj))) {
        message(
            sum(is.na(padj)), " gene(s) with NA adjusted p-value ",
            "(not passing independent filtering or flagged as outliers) are excluded."
        )
    }

    keep <- which(!is.na(padj) & padj < cutoff)

    if (length(keep) > 0 && (direction != "both" || minAbsLFC > 0)) {
        lfc <- .get_results_column(x, "log2FoldChange")
        if (is.null(lfc)) {
            stop("no 'log2FoldChange' column found in the DESeqResults object.", call. = FALSE)
        }
        lfc <- lfc[keep]
        ok <- rep(TRUE, length(keep))
        if (direction == "up") {
            ok <- ok & lfc > 0
        } else if (direction == "down") {
            ok <- ok & lfc < 0
        }
        if (minAbsLFC > 0) {
            ok <- ok & abs(lfc) > minAbsLFC
        }
        keep <- keep[ok & !is.na(ok)]
    }

    out <- ids[keep]
    if (length(out) == 0) {
        warning(
            "no gene has padj < ", cutoff,
            if (direction != "both") paste0(" with direction = '", direction, "'"),
            if (minAbsLFC > 0) paste0(" and abs(log2FoldChange) > ", minAbsLFC),
            "; returning an empty vector.",
            call. = FALSE
        )
    }
    out
}


## ---- DseaData container ---------------------------------------------------
## as_elist() / limma EList were removed (see NEWS 0.2.6.9000): dsea() is
## zero-dependency, so the container is the package's own DseaData list. as_dsea()
## replaces as_elist(); the bridge keeps DESeq2's fitted mu/dispersion as inputs.

#' Convert a DESeq2 dataset into a DseaData container
#'
#' Build a log2 expression matrix and per-gene precision weights from a
#' fitted \code{DESeqDataSet}, for use with \code{\link{dsea}}.
#'
#' @param x A \code{\link[DESeq2]{DESeqDataSet}} that has been fitted with
#'   \code{\link[DESeq2]{DESeq}} (so that \code{assay(x, "mu")} and
#'   \code{mcols(x)$dispersion} are populated).
#' @param norm.factors Sample-wise normalisation factors. Defaults to
#'   \code{DESeq2::sizeFactors(x)}. These are \emph{depth-like} factors (see
#'   Details), so they are rescaled to the mean sequencing depth before being
#'   used as a denominator. Pass \code{NULL} to use raw sequencing depth only.
#' @param mu.floor Lower bound applied to the fitted means before forming the
#'   delta-method weights. A single positive finite number; default \code{1}.
#' @param ... Ignored, for S3 generic consistency.
#'
#' @return An object of S3 class \code{"DseaData"}: a list with elements
#' \describe{
#'   \item{\code{E}}{G x n matrix of log2 expression values.}
#'   \item{\code{weights}}{G x n matrix of precision weights, or \code{NULL}.}
#'   \item{\code{design}}{The n x p design matrix, with \code{colnames} used
#'     for contrast resolution.}
#' }
#'
#' @details
#' \strong{Why the library size is \code{mean(colSums(ct)) * norm.factors}.}
#'
#' \code{DESeq2::sizeFactors()} are estimated by the median-of-ratios method and
#' are \emph{themselves a measure of sequencing depth}, merely rescaled so that
#' their geometric mean is 1:
#'
#' \preformatted{
#'   airway: exp(mean(log(sizeFactors))) = 1.0048
#'           cor(log sizeFactors, log colSums) = 0.983
#'           sizeFactors range 0.670 - 1.399 (2.09x)
#'           colSums     range 15.2M - 30.8M (2.03x)   <- the same quantity
#'   simulated data with depth decoupled from group:
#'           cor(log sizeFactors, log depth) = 0.9999
#' }
#'
#' They are \strong{not} edgeR-style relative correction factors. DESeq2's
#' model is \code{count_ij ~ sf_j * q_ij}, \code{sf_j} being the depth. Writing
#' that into the expression matrix, the old denominator
#' \code{colSums(ct) * norm.factors} cancels \code{sf} against itself and
#' leaves a per-sample offset that the linear design cannot absorb:
#'
#' \preformatted{
#'   E_ij ~ log2(q_ij) + log2(sf_j) - log2(colSums_j) - log2(sf_j)
#'        = log2(q_ij) - log2(colSums_j)
#' }
#'
#' Dividing by \code{mean(colSums(ct)) * norm.factors} instead keeps the CPM
#' scale while removing exactly one copy of the depth. Measured on \code{airway}
#' (\code{~ cell + dex}), same delta-method weights, only the denominator
#' changed:
#'
#' \preformatted{
#'   metric                                colSums*sf   mean(colSums)*sf
#'   unmodelled per-sample sd (log2)           0.3440            0.1267
#'   median per-gene residual sigma            0.1690            0.1378
#'   cor(colMeans(E), -log2(colSums))          0.9975            0.9829
#' }
#'
#' The residual per-sample offset drops by a factor of 2.7. Note also that
#' \code{colMeans(E)} remains strongly correlated with \code{-log2(colSums)}
#' (\eqn{r = 0.983}) -- the depth term is genuinely present in the data; it is
#' the \emph{unmodelled} part that the design cannot explain, and only 52\% of
#' \code{log2(colSums)} is explainable by the design, so that part has nowhere
#' to go. Because the offset has the same sign for every gene, it acts as a
#' rank-1 sample effect: competitive statistics (camera) centre it away, while
#' fry/roast, which do not centre across the genome, inherit it.
#'
#' \strong{Delta-method weights.}
#'
#' With \eqn{\tilde\mu = \max(\mu, \code{mu.floor})}, the weights are the
#' reciprocal of the delta-method variance of
#' \eqn{\log_2((c + 0.5) / D)}:
#'
#' \preformatted{
#'   w = (mu_tilde + 0.5)^2 * log(2)^2 / (mu_tilde + alpha * mu_tilde^2)
#' }
#'
#' \code{mu} is taken from \code{assay(x, "mu")} and \code{alpha} from
#' \code{mcols(x)$dispersion}; both require a prior \code{DESeq2::DESeq()} run.
#'
#' Because \code{mu} is \code{NA} for genes dropped by DESeq2's independent
#' filtering and the raw formula yields non-finite weights there, non-finite
#' weights are replaced by the minimum finite weight. This replacement is
#' required: downstream weighted least squares rejects non-finite weights.
#'
#' @examples
#' \dontrun{
#' library(DESeq2)
#' dds <- DESeq(dds)
#' dd <- as_dsea(dds)
#' dsea(dd, gene_sets, contrast = "dex_untrt_vs_trt")
#' }
#'
#' @seealso \code{\link{dsea}}
#' @export
as_dsea <- function(x, ...) {
    UseMethod("as_dsea")
}

#' @rdname as_dsea
#' @export
as_dsea.DESeqDataSet <- function(x, norm.factors = DESeq2::sizeFactors(x),
                                 mu.floor = 1, ...) {
    ## ---- validate mu.floor -------------------------------------------------
    if (!is.numeric(mu.floor) || length(mu.floor) != 1L ||
        !is.finite(mu.floor) || mu.floor <= 0) {
        stop("'mu.floor' must be a single positive finite number.", call. = FALSE)
    }

    ## ---- required DESeq2 quantities ---------------------------------------
    ## Both are produced by DESeq2::DESeq(); without them the delta-method
    ## weights are undefined, so fail loudly rather than silently degrading.
    mu <- tryCatch(SummarizedExperiment::assay(x, "mu"), error = function(e) NULL)
    if (is.null(mu)) {
        stop("assay(x, 'mu') is not available; run DESeq2::DESeq(x) first ",
             "so that the fitted means are stored.", call. = FALSE)
    }

    alpha <- SummarizedExperiment::mcols(x)$dispersion
    if (is.null(alpha)) {
        stop("mcols(x)$dispersion is not available; run DESeq2::DESeq(x) first ",
             "so that per-gene dispersions are estimated.", call. = FALSE)
    }

    ct <- DESeq2::counts(x)

    ## ---- library size ------------------------------------------------------
    ## norm.factors are DEPTH-LIKE (median-of-ratios, geometric mean ~1), so
    ## they are rescaled to the mean depth rather than multiplied by the raw
    ## per-sample depth (which would count depth twice). See @details.
    if (is.null(norm.factors)) {
        lib.size <- rep(mean(colSums(ct)), ncol(ct))
    } else {
        if (!is.numeric(norm.factors) || length(norm.factors) != ncol(ct)) {
            stop("'norm.factors' must be NULL or a numeric vector with one ",
                 "value per sample (", ncol(ct), ").", call. = FALSE)
        }
        lib.size <- mean(colSums(ct)) * norm.factors
    }

    ## Explicit byrow = TRUE: `matrix * length-n vector` recycles COLUMN-major,
    ## which silently broadcasts onto genes whenever nrow(ct) is a multiple of
    ## ncol(ct). Same class of bug as the batch-effect one in the benchmark.
    denom <- matrix(lib.size, nrow(ct), ncol(ct), byrow = TRUE) + 1

    E <- log2((ct + 0.5) / denom * 1e6)
    dimnames(E) <- dimnames(ct)

    ## ---- delta-method precision weights ------------------------------------
    mu_tilde <- pmax(mu, mu.floor)
    w <- (mu_tilde + 0.5)^2 * log(2)^2 /
        (mu_tilde + alpha * mu_tilde^2)

    bad <- !is.finite(w)
    if (any(bad)) {
        ## mu is NA for genes dropped by independent filtering; the formula
        ## then yields NaN. Fall back to the smallest finite weight.
        fin <- w[!bad]
        w[bad] <- if (length(fin)) min(fin) else 0
    }
    dimnames(w) <- dimnames(E)

    ## ---- design ------------------------------------------------------------
    design <- attr(x, "modelMatrix")
    if (is.null(design)) {
        design <- stats::model.matrix(DESeq2::design(x),
                                     SummarizedExperiment::colData(x))
    }
    design <- as.matrix(design)

    structure(
        list(E = E, weights = w, design = design),
        class = "DseaData"
    )
}
