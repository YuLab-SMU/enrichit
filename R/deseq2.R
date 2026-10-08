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

## ---- EList bridge ---------------------------------------------------------

#' Convert a DESeqDataSet to a limma EList with delta-method precision weights
#'
#' Build a limma \code{EList} (log-expression matrix plus per-observation
#' precision weights and the design matrix) from a fitted
#' \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSet}}, ready for limma's
#' weighted linear-model machinery and gene-set tests such as
#' \code{\link[limma:camera]{camera}}, \code{\link[limma:roast]{roast}} and
#' \code{\link[limma:fry]{fry}} (see \code{\link{set_test}}).
#'
#' @details The expression matrix follows the voom convention,
#'   \code{E[i,j] = log2((count[i,j] + 0.5) / (lib.size[j] * norm.factor[j] + 1) * 1e6)}
#'   with \code{lib.size = colSums(counts)}.
#'
#'   The precision weights, however, are not re-learned from an empirical
#'   mean-variance trend: they reuse the negative-binomial mean-variance law
#'   that DESeq2 has already estimated genome-wide (fitted means in
#'   \code{assay(x, "mu")}, empirical-Bayes-shrunk dispersions in
#'   \code{mcols(x)$dispersion}). The delta method on
#'   \code{log2(c + 0.5)} gives
#'   \code{Var(y[i,j]) ~ (mu[i,j] + alpha[i] * mu[i,j]^2) / ((mu[i,j] + 0.5)^2 * log(2)^2)},
#'   and the weight is the reciprocal of this variance. Global scaling of
#'   the weights is absorbed by the residual-variance estimate downstream,
#'   so the \code{log(2)^2} factor only fixes the absolute scale.
#'
#'   As \code{mu[i,j]} approaches zero the delta-method variance vanishes
#'   while the discreteness scatter of near-zero counts does not, so fitted
#'   means are floored at \code{mu.floor} (default 1 count) before the
#'   weights are computed. Genes without an estimated dispersion get
#'   \code{alpha = 1}; all-zero genes additionally have NA fitted means,
#'   and any non-finite weights are replaced by the minimum finite
#'   weight, with a message.
#'
#'   DESeq2 must have been run on the object: both the \code{"mu"} assay
#'   and the dispersions are required.
#'
#' @param x a fitted \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSet}}.
#' @param norm.factors numeric vector of length \code{ncol(x)} of
#'   normalization factors multiplied into the library sizes; defaults to
#'   \code{DESeq2::sizeFactors(x)}. \code{NULL} means no normalization
#'   (raw library sizes only).
#' @param mu.floor positive numeric, fitted means are floored at this value
#'   before computing weights (default 1).
#' @param ... currently unused.
#' @return An object of class \code{\link[limma:EList-class]{EList}} with
#'   components \code{E} (genes x samples log-expression matrix),
#'   \code{weights} (genes x samples precision weights) and \code{design}
#'   (the model matrix used by \code{DESeq2::DESeq()}, with columns named
#'   as in \code{DESeq2::resultsNames()}, e.g.
#'   \code{condition_trt_vs_ctrl}).
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
#' elist <- as_elist(dds)
#' }
#' @export
as_elist <- function(x, ...) {
    UseMethod("as_elist")
}

#' @rdname as_elist
#' @export
as_elist.default <- function(x, ...) {
    stop(
        "no as_elist() method for objects of class '",
        paste(class(x), collapse = ", "),
        "'; currently DESeqDataSet is supported.",
        call. = FALSE
    )
}

#' @rdname as_elist
#' @export
as_elist.DESeqDataSet <- function(x,
                                  norm.factors = DESeq2::sizeFactors(x),
                                  mu.floor = 1,
                                  ...) {
    rlang::check_installed("DESeq2", "for converting DESeq2 objects.")
    rlang::check_installed("limma", "for constructing an EList.")

    if (is.null(mu.floor) || !is.numeric(mu.floor) || length(mu.floor) != 1 ||
        !is.finite(mu.floor) || mu.floor <= 0) {
        stop("mu.floor must be a single positive finite number.", call. = FALSE)
    }

    if (!"mu" %in% SummarizedExperiment::assayNames(x)) {
        stop(
            "no 'mu' assay found; run DESeq2::DESeq() on the object first -- ",
            "as_elist() reuses its fitted means and dispersions.",
            call. = FALSE
        )
    }

    ct <- DESeq2::counts(x)
    ids <- rownames(ct)
    if (is.null(ids)) {
        stop("the DESeqDataSet has no row names; gene identifiers are required.", call. = FALSE)
    }
    n <- ncol(ct)
    if (is.null(norm.factors)) {
        norm.factors <- rep(1, n)
    }
    if (!is.numeric(norm.factors) || length(norm.factors) != n ||
        anyNA(norm.factors) || any(norm.factors <= 0)) {
        stop("norm.factors must be NULL or a positive numeric vector of length ncol(x).", call. = FALSE)
    }
    lib.size <- colSums(ct) * norm.factors
    if (any(lib.size <= 0)) {
        stop("all library sizes must be positive.", call. = FALSE)
    }

    E <- log2((ct + 0.5) / (matrix(lib.size, nrow(ct), n, byrow = TRUE) + 1) * 1e6)

    ## DESeq2 stores the model matrix it fitted with, with columns renamed
    ## to the resultsNames() convention (e.g. condition_trt_vs_ctrl);
    ## a fresh model.matrix() would give unrenamed names (conditiontrt)
    design <- attr(x, "modelMatrix")
    if (is.null(design)) {
        design <- stats::model.matrix(
            DESeq2::design(x),
            data = as.data.frame(SummarizedExperiment::colData(x))
        )
    }

    alpha <- SummarizedExperiment::mcols(x)$dispersion
    if (is.null(alpha)) {
        stop("no dispersions found; run DESeq2::DESeq() on the object first.", call. = FALSE)
    }
    mu <- SummarizedExperiment::assay(x, "mu")

    w <- .deseq2_precision_weights(mu, alpha, mu.floor = mu.floor)
    dimnames(w) <- dimnames(E)

    methods::new("EList", list(E = E, weights = w, design = design))
}

.deseq2_precision_weights <- function(mu, alpha, mu.floor) {
    alpha[is.na(alpha)] <- 1
    mu_f <- pmax(mu, mu.floor)
    w <- (mu_f + 0.5)^2 * log(2)^2 / (mu_f + alpha * mu_f^2)
    bad <- !is.finite(w)
    if (any(bad)) {
        message(
            "replacing ", sum(bad),
            " non-finite precision weight(s) with the minimum finite weight."
        )
        w[bad] <- min(w[!bad])
    }
    w
}
