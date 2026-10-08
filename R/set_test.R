## Set-level competitive gene-set tests on DESeq2 / EList input.
##
## `set_test()` dispatches to limma's camera / roast / fry. A fitted
## DESeqDataSet is converted by `as_elist()`, which reuses DESeq2's
## genome-wide fitted means and dispersions as per-observation precision
## weights, so no separate voom step is needed. limma is a lazy Suggests
## dependency and provides the reference backend; a native engine may
## replace it later.

#' Competitive gene-set tests on a DESeqDataSet or EList
#'
#' Run limma's competitive gene-set tests -- \code{camera},
#' \code{roast} or \code{fry} -- directly on a fitted
#' \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSet}}, or on any
#' \code{\link[limma:EList-class]{EList}} such as the one produced by
#' \code{\link{as_elist}}.
#'
#' @details Competitive tests ask whether the genes of a set shift
#'   consistently, relative to all other genes, under a contrast of the
#'   experimental design. This is the design-aware, no-cut-off alternative
#'   to \code{\link{ora}} (which requires a significance cut-off) and to
#'   \code{\link{gsea}} (which reduces the design to a single ranked list):
#'   any contrast expressible in the DESeq2 design formula can be tested,
#'   and the inter-gene correlation is accounted for as in limma.
#'
#'   For a \code{DESeqDataSet}, the object is converted with
#'   \code{\link{as_elist}}, whose per-observation precision weights are
#'   derived from DESeq2's own fitted means and dispersions rather than
#'   from a re-learned mean-variance trend. Arguments in \code{...} are
#'   forwarded to the underlying limma function (e.g.
#'   \code{inter.gene.cor}, \code{use.ranks}, \code{nrot}, \code{robust});
#'   to tune the weight construction instead, call
#'   \code{\link{as_elist}} explicitly and pass the resulting
#'   \code{EList}.
#'
#'   Gene sets with fewer than 2 genes present in the data are dropped,
#'   with a warning.
#'
#' @param x a fitted \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSet}}
#'   (must have been processed by \code{DESeq2::DESeq()}) or an
#'   \code{\link[limma:EList-class]{EList}} with a \code{design} component.
#' @param gene_sets named list of character vectors of gene identifiers
#'   matching the row names of the data; a single unnamed character vector
#'   is treated as one set.
#' @param contrast the contrast to test. Either a design-column name
#'   (character length 1), a numeric contrast vector of length
#'   \code{ncol(design)}, a DESeq2-style three-token contrast
#'   \code{c(factor, numerator, denominator)} which is mapped to the
#'   corresponding design column, or \code{NULL} (default) to test the
#'   last column of the design, matching limma's and DESeq2's defaults.
#' @param method character, the limma test to run: \code{"camera"}
#'   (default), \code{"roast"} or \code{"fry"}.
#' @param ... additional arguments forwarded to the selected limma
#'   function.
#' @return The result object of the selected limma test, unchanged:
#'   a data frame for \code{camera} and \code{fry}, a matrix for
#'   \code{roast}, with one row per gene set.
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
#' gene_sets <- list(
#'     SET1 = paste0("gene", 1:50),
#'     SET2 = paste0("gene", 201:260)
#' )
#' res <- set_test(dds, gene_sets, contrast = c("condition", "B", "A"))
#' }
#' @export
set_test <- function(x, gene_sets, contrast = NULL,
                     method = c("camera", "roast", "fry"), ...) {
    method <- match.arg(method)
    rlang::check_installed("limma", "for set_test().")

    if (inherits(x, "EList")) {
        elist <- x
    } else {
        elist <- as_elist(x)
    }

    design <- elist$design
    if (is.null(design)) {
        stop(
            "the EList does not carry a design matrix; build it with ",
            "as_elist() or assign elist$design directly.",
            call. = FALSE
        )
    }

    idx <- .validate_gene_sets(gene_sets, rownames(elist$E))
    ct <- .resolve_contrast(contrast, design)

    switch(method,
        camera = limma::camera(elist, idx, design = design, contrast = ct, ...),
        roast = limma::roast(elist, idx, design = design, contrast = ct, ...),
        fry = limma::fry(elist, idx, design = design, contrast = ct, ...)
    )
}

.validate_gene_sets <- function(gene_sets, ids) {
    if (!is.list(gene_sets)) {
        gene_sets <- list(set1 = gene_sets)
    }
    if (length(gene_sets) == 0) {
        stop("gene_sets is empty.", call. = FALSE)
    }
    if (is.null(names(gene_sets))) {
        names(gene_sets) <- paste0("set", seq_along(gene_sets))
    }
    keep <- vapply(
        gene_sets,
        function(s) sum(!is.na(s) & s %in% ids) >= 2L,
        logical(1)
    )
    if (!any(keep)) {
        stop(
            "none of the gene sets has at least 2 genes present in the data; ",
            "check that the identifiers match the row names.",
            call. = FALSE
        )
    }
    if (!all(keep)) {
        warning(
            sum(!keep), " gene set(s) dropped: fewer than 2 genes present ",
            "in the data.",
            call. = FALSE
        )
    }
    gene_sets[keep]
}

.resolve_contrast <- function(contrast, design) {
    if (is.null(contrast)) {
        return(ncol(design))
    }
    if (is.character(contrast) && length(contrast) == 3L) {
        col <- paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3])
        if (col %in% colnames(design)) {
            return(col)
        }
        stop(
            "contrast c(", paste(contrast, collapse = ", "),
            ") maps to design column '", col, "', which does not exist; ",
            "available columns: ", paste(colnames(design), collapse = ", "),
            ".",
            call. = FALSE
        )
    }
    contrast
}
