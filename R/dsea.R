# dsea.R -- user-facing entry point for the native set-level test engine
#
# This file is the R half of the dsea() pipeline.  Everything statistical
# lives in src/set_test.cpp; what happens here is input normalisation,
# contrast resolution and argument validation.
#
# This package has exactly one engine, its own C++ one (commit ca09bf0 removed
# the third-party alternative).  The `engine` argument is kept so that old
# scripts fail with an actionable message instead of an "unused argument"
# error, but "native" is the only accepted value.
#
# Dispatch uses getFromNamespace() so the C++ entry points are resolved at
# call time rather than build time.

# ---------------------------------------------------------------------------
# C++ entry point lookup
# ---------------------------------------------------------------------------

## The three set-test entry points are not exported from the package
## namespace: they take raw index lists and a pre-resolved contrast column, and
## are not part of the user-facing API.  Each one is bound in the namespace by
## RcppExports.R, so getFromNamespace() finds it after load_all().

.get_cpp <- function(name, hint) {
    fn <- tryCatch(getFromNamespace(name, "enrichit"),
                   error = function(e) NULL)
    if (is.null(fn)) {
        stop(sprintf("internal error: '%s' is not available in this build%s",
                     name, hint), call. = FALSE)
    }
    fn
}

.cpp_hint <- paste0(
    ".\n  The native set-test engine is built from src/set_test.cpp; ",
    "reinstall or run devtools::load_all('.', recompile = TRUE)."
)

# ---------------------------------------------------------------------------
# Contrast resolution
# ---------------------------------------------------------------------------

#' Resolve a contrast specification to a design-matrix column
#'
#' Internal helper used by \code{\link{dsea}} to turn the three accepted
#' \code{contrast} spellings into the 1-based column index that the C++
#' engine expects.
#'
#' @param contrast A design column name (character), a 1-based column index
#'   (numeric of length 1), a length-\code{p} contrast vector, or \code{NULL}.
#' @param design The \code{n x p} design matrix.
#' @return A single integer: the 1-based index of \code{design} to use as the
#'   contrast column.
#'
#' @details
#' The returned index encodes the \emph{signed} comparison, not just which
#' column of \code{design} is of interest. The C++ engine puts the named
#' column last and takes the sign of its (signed) code from the QR pivot, so
#' the effects coordinate is the log2 fold change of that column's level
#' relative to the reference level. For a numeric contrast vector with a
#' single non-zero entry \eqn{c_i}:
#'
#' \itemize{
#'   \item \eqn{c_i > 0}: the contrast asks for \eqn{+\beta_i}, which is
#'     exactly what the engine computes for column \eqn{i}; return \eqn{i}.
#'   \item \eqn{c_i < 0}: the contrast asks for \eqn{-\beta_i}. Rather than
#'     negating the (unnamed) effect columns and guessing the statistic's
#'     sign convention downstream, the design column is \emph{reflected}
#'     before the call: returning the index of the \emph{indicator} column of
#'     the level that the contrast names as the numerator makes the engine
#'     compute the requested difference directly. This requires the caller to
#'     hand the engine the reflected contrast column; see
#'     \code{\link{.dsea_design}}.
#' }
#'
#' A length-1 negative integer is treated the same way (a column index with a
#' minus sign is how some callers spell "the other direction").
#'
#' @seealso \code{\link{.dsea_design}}
#' @noRd
.resolve_contrast_idx <- function(contrast, design) {
    design <- as.matrix(design)
    p <- ncol(design)
    if (p < 1L) {
        stop("'design' must have at least one column.", call. = FALSE)
    }
    cn <- colnames(design)
    if (is.null(cn)) cn <- rep(NA_character_, p)

    ## ---- NULL: only unambiguous when there is a single non-intercept column
    if (is.null(contrast)) {
        cand <- which(!(cn %in% c("(Intercept)", "Intercept")))
        if (length(cand) == 1L) {
            return(list(idx = as.integer(cand), negate = FALSE))
        }
        if (length(cand) == 0L && p == 1L) {
            return(list(idx = 1L, negate = FALSE))
        }
        stop("'contrast' is NULL but the design has ", length(cand),
             " non-intercept columns (",
             if (any(!is.na(cn))) paste(cn, collapse = ", ") else "unnamed",
             "); specify one of them by name.", call. = FALSE)
    }

    ## ---- character: design column name
    if (is.character(contrast)) {
        if (length(contrast) != 1L || is.na(contrast) || !nzchar(contrast)) {
            stop("'contrast' must be a single non-empty column name.",
                 call. = FALSE)
        }
        hit <- match(contrast, cn)
        if (is.na(hit)) {
            stop("'contrast' = \"", contrast, "\" is not a column of the ",
                 "design matrix. Available columns: ",
                 if (any(!is.na(cn))) paste(cn, collapse = ", ")
                 else "<unnamed>.", call. = FALSE)
        }
        return(list(idx = as.integer(hit), negate = FALSE))
    }

    ## ---- numeric
    if (is.numeric(contrast)) {
        negate <- FALSE
        if (length(contrast) == 1L) {
            ci <- contrast
            if (!is.finite(ci) || ci != as.integer(ci)) {
                stop("a length-1 numeric 'contrast' must be an integer column ",
                     "index.", call. = FALSE)
            }
            ci <- as.integer(ci)
            if (ci < 0L) {
                ci <- -ci
                negate <- TRUE
            }
            if (ci < 1L || ci > p) {
                stop("'contrast' = ", contrast, " is out of range; the design ",
                     "has ", p, " column", if (p == 1L) "" else "s", ".",
                     call. = FALSE)
            }
            return(list(idx = ci, negate = negate))
        }

        if (length(contrast) != p) {
            stop("a contrast vector must have one entry per design column (",
                 p, "), got ", length(contrast), ".", call. = FALSE)
        }
        if (any(!is.finite(contrast))) {
            stop("'contrast' contains non-finite values.", call. = FALSE)
        }

        nz <- which(contrast != 0)
        if (length(nz) == 1L) {
            return(list(idx = as.integer(nz), negate = contrast[nz] < 0))
        }
        if (length(nz) == 0L) {
            stop("'contrast' is all zeros; it defines an empty comparison.",
                 call. = FALSE)
        }
        stop("'contrast' must be a single design column name, a column index, ",
             "or a vector with exactly one non-zero entry (got ", length(nz),
             " non-zero entries). Contrasts that span several columns are not ",
             "supported by the native set-test engine: refit with\n",
             "  DESeq2::results(dds, contrast = ...)\n",
             "or form the comparison as a single design column first.",
             call. = FALSE)
    }

    stop("'contrast' must be NULL, a column name (character), a column index ",
         "(integer) or a numeric contrast vector; got class '",
         paste(class(contrast), collapse = "/"), "'.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Design reflection for a negated contrast
# ---------------------------------------------------------------------------

#' Reflect a design matrix for a negated contrast
#'
#' Applies the \eqn{\pm 1} diagonal transformation that turns a "test
#' \eqn{-\beta_i}" request into a "test \eqn{+\beta_{i'}}" request on the same
#' column space, which is what the C++ engine can compute directly.
#'
#' @param design The \eqn{n \times p} design matrix.
#' @param idx 1-based column index, as returned by
#'   \code{\link{.resolve_contrast_idx}}.
#' @param negate Logical: whether to reflect; when \code{TRUE} all columns
#'   other than \code{idx} are negated.
#' @return The (possibly reflected) design matrix.
#' @noRd
.dsea_design <- function(design, idx, negate) {
    p <- ncol(design)
    if (!isTRUE(negate)) {
        return(design)
    }
    ## -1 on every other column: X -> X diag(1,..,-1@i,..,1). The column space
    ## - and therefore the residual coordinates and their degrees of freedom -
    ## is unchanged; only the sign of the contrast coefficient flips.
    s <- rep(1, p)
    s[-idx] <- -1
    reflected <- sweep(design, 2L, s, `*`)
    colnames(reflected) <- colnames(design)
    reflected
}

# ---------------------------------------------------------------------------
# Container normalisation
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.as_dsea_data <- function(x, ...) {
    ## A fitted DESeqDataSet (or anything as_dsea() understands) is funnelled
    ## through the S3 generic; the DseaData object and the bare list are used
    ## as they are.
    if (inherits(x, "DseaData")) {
        dd <- x
    } else if (is.list(x) && all(c("E", "design") %in% names(x))) {
        ## Backward compatibility: bare list(E=, weights=, design=).
        dd <- structure(
            list(E = x$E,
                 weights = if (is.null(x$weights)) NULL else x$weights,
                 design = x$design),
            class = "DseaData"
        )
    } else {
        dd <- as_dsea(x, ...)
    }

    if (!is.list(dd) || is.null(dd$E) || is.null(dd$design)) {
        stop("'x' must provide both an expression matrix 'E' and a 'design' ",
             "matrix.", call. = FALSE)
    }

    E <- as.matrix(dd$E)
    storage.mode(E) <- "double"
    if (!is.numeric(E)) {
        stop("'E' must be a numeric matrix of log2 expression values.",
             call. = FALSE)
    }

    design <- as.matrix(dd$design)
    storage.mode(design) <- "double"

    if (nrow(design) != ncol(E)) {
        stop("nrow(design) (", nrow(design), ") must equal ncol(E) (",
             ncol(E), ").", call. = FALSE)
    }
    if (is.null(rownames(E))) {
        stop("'E' must have rownames (gene identifiers).", call. = FALSE)
    }

    w <- dd$weights
    if (!is.null(w)) {
        w <- as.matrix(w)
        storage.mode(w) <- "double"
        if (!identical(dim(w), dim(E))) {
            stop("'weights' must have the same dimensions as 'E' (",
                 paste(dim(E), collapse = " x "), ").", call. = FALSE)
        }
    }

    list(E = E, weights = w, design = design)
}

## Map a user-supplied gene set (character vector, or list of them) onto
## 1-based row indices of E.
##
## NOTE ON INDEXING. `contrast_idx` and `index` are both 1-based at this
## boundary: the C++ helper st_resolve_set() subtracts 1 from every entry it is
## handed (`iset[j] - 1`), i.e. it takes R-style 1-based row numbers, exactly
## as the reference set loops do. Feeding it 0-based indices shifts every set by one
## gene and pushes the last row of the matrix out of range.
.build_index <- function(gene_sets, gene_ids) {
    if (length(gene_sets) == 0L) {
        stop("'gene_sets' is empty.", call. = FALSE)
    }
    if (is.character(gene_sets)) {
        gene_sets <- list(gene_sets)
    }
    if (!is.list(gene_sets)) {
        stop("'gene_sets' must be a list of character vectors, or a single ",
             "character vector.", call. = FALSE)
    }
    nms <- names(gene_sets)
    if (is.null(nms) || any(!nzchar(nms))) {
        nms <- paste0("set", seq_along(gene_sets))
    }
    nms <- make.unique(nms, sep = "_")

    sets <- lapply(gene_sets, function(s) {
        s <- as.character(s)
        s <- s[!is.na(s)]
        if (is.null(s)) character(0) else unique(s)
    })
    names(sets) <- nms

    index <- lapply(sets, function(s) {
        i <- match(s, gene_ids)
        as.integer(i[!is.na(i)])
    })
    names(index) <- nms

    n <- vapply(index, length, integer(1))
    index <- index[n > 0L]
    if (length(index) == 0L) {
        stop("none of the gene set members are present in 'x'.", call. = FALSE)
    }
    index
}

## ---------------------------------------------------------------------------
## Per-method argument resolution
## ---------------------------------------------------------------------------

.nrot_value <- function(nrot) {
    if (!is.numeric(nrot) || length(nrot) != 1L || !is.finite(nrot) ||
        nrot < 1 || nrot != as.integer(nrot)) {
        stop("'nrot' must be a single positive integer (number of rotations).",
             call. = FALSE)
    }
    as.integer(nrot)
}

.gene_weights_value <- function(gene.weights, ids, method) {
    if (is.null(gene.weights)) {
        return(NULL)
    }
    if (!is.numeric(gene.weights) || is.null(names(gene.weights))) {
        stop("'gene.weights' must be a named numeric vector for method = \"",
             method, "\".", call. = FALSE)
    }
    gw <- gene.weights[ids]
    if (anyNA(gw)) {
        stop("'gene.weights' does not cover all genes in 'x' (",
             sum(is.na(gw)), " missing).", call. = FALSE)
    }
    gw <- as.numeric(gw)
    if (any(!is.finite(gw))) {
        stop("'gene.weights' contains non-finite values.", call. = FALSE)
    }
    if (all(gw == gw[1])) {
        ## Constant weights carry no information, so they are dropped.
        return(NULL)
    }
    gw
}

# ---------------------------------------------------------------------------
# dsea()
# ---------------------------------------------------------------------------

#' Set-level differential expression testing for RNA-seq data
#'
#' Test whether gene sets are shifted between the two groups of a
#' \code{DESeq2} contrast, using one of three self-contained competitive
#' set tests: \code{camera}, \code{roast} or \code{fry}. All three are
#' implemented in C++ in this package; no third-party statistics package is
#' called.
#'
#' @param x The data to test. One of
#'   \describe{
#'     \item{a fitted \code{\link[DESeq2]{DESeqDataSet}}}{(after
#'       \code{DESeq2::DESeq()}); converted on the fly with
#'       \code{\link{as_dsea}}.}
#'     \item{a \code{"DseaData"} object}{as returned by
#'       \code{\link{as_dsea}}.}
#'     \item{a bare list}{with elements \code{E} (genes x samples log2
#'       expression matrix), \code{weights} (optional) and \code{design}.
#'       Kept for backward compatibility with pre-0.2.7 scripts.}
#'   }
#' @param gene_sets A named list of character vectors of gene identifiers, or
#'   a single character vector (treated as one set). Identifiers are matched
#'   against \code{rownames(E)}; members that are absent are dropped.
#' @param contrast Which design-matrix column encodes the comparison. Accepts
#'   the DESeq2-style contrast name (e.g. \code{"dex_trt_vs_untrt"}), a column
#'   index, or a numeric vector with exactly one non-zero entry. \code{NULL}
#'   is allowed only when the design has a single non-intercept column.
#'   Note that the design column names have the sign implicit in the column
#'   itself: the requested direction is read from the column, not from the
#'   order of the name.
#' @param method Which set test to run: \code{"camera"}, \code{"roast"} or
#'   \code{"fry"}. See Details.
#' @param engine The statistical engine. Only \code{"native"} is supported.
#'   The argument exists so that scripts written against the older two-engine
#'   API fail with an actionable message instead of an "unused argument"
#'   error; any other value raises an error.
#' @param sort Logical: sort the result by \code{PValue} (default \code{TRUE}).
#' @param ... Method-specific arguments. These are passed through \code{...}
#'   rather than named in the signature, so they do not appear in the usage
#'   section above; they are validated against the selected \code{method} and
#'   an unknown name is an error. For \code{method = "camera"}:
#'   \code{inter.gene.cor} (default \code{0.01}; \code{NA} estimates it) and
#'   \code{allow.neg.cor} (default \code{FALSE}). For
#'   \code{method = "roast"}: \code{nrot} (default \code{999}),
#'   \code{set.statistic} (default \code{"mean"}), \code{midp} (default
#'   \code{FALSE}) and \code{gene.weights} (default \code{NULL}). For
#'   \code{method = "fry"}: \code{standardize} (default
#'   \code{"posterior.sd"}) and \code{gene.weights} (default \code{NULL}).
#' \code{inter.gene.cor}: (\code{camera}) Inter-gene correlation used to set
#'   the variance inflation factor. Default \code{0.01}. Pass \code{NA} to
#'   estimate it from the data.
#' \code{allow.neg.cor}: (\code{camera}) Allow a negative inter-gene
#'   correlation estimate. Default \code{FALSE}; a negative estimate is then
#'   treated as zero.
#' \code{nrot}: (\code{roast}) Number of rotations. Default \code{999}.
#' \code{set.statistic}: (\code{roast}) Set-level statistic: one of
#'   \code{"mean"}, \code{"floormean"}, \code{"mean50"} or \code{"msq"}.
#'   Default \code{"mean"}.
#' \code{midp}: (\code{roast}) Accepted for drop-in compatibility, but it
#'   \strong{cannot change any value this function returns}. In
#'   \code{limma::mroast} the mid-p correction is applied to the multiple-testing
#'   columns only (\code{tab$FDR} and \code{tab$FDR.Mixed} are computed from
#'   \code{PValue - 1/2/(nrot + 1)} and then floored back at the unadjusted
#'   \code{PValue}); the raw \code{PValue} column is untouched. \code{dsea()}
#'   returns unadjusted p-values and no FDR column at all -- the caller chooses the
#'   testing family (see the \emph{Value} section) -- so there is nothing for
#'   \code{midp} to act
#'   on. Verified: with an identical seed, \code{midp = TRUE} and
#'   \code{midp = FALSE} return \code{identical()} p-values, which matches limma's
#'   own \code{PValue} under both settings. If you need limma's mid-p FDR, apply
#'   \code{p.adjust(PValue - 1/2/(nrot + 1), "BH")} yourself and floor it at
#'   \code{PValue}.
#'   Default \code{FALSE} (limma's own default is \code{TRUE}).
#' \code{standardize}: (\code{fry}) Standardisation of the residual variances:
#'   one of \code{"posterior.sd"}, \code{"none"}, \code{"residual.sd"} or
#'   \code{"p2"}. Default \code{"posterior.sd"}, which is the only choice that
#'   borrows strength across genes.
#' \code{gene.weights}: (\code{roast}, \code{fry}) Optional named numeric
#'   vector of per-gene weights, covering every gene in \code{x}. Default
#'   \code{NULL} (unweighted). A constant vector is ignored.
#' @param verbose Logical: emit progress messages. Default \code{FALSE}.
#'
#' @return A \code{data.frame} with one row per tested gene set and columns
#'   \describe{
#'     \item{\code{set}}{Gene set name.}
#'     \item{\code{NGenes}}{Number of set members found in \code{x}.}
#'     \item{\code{Direction}}{\code{"Up"} if the set's genes are on average
#'       higher in the numerator group of the contrast, \code{"Down"} if
#'       lower, and \code{NA} when the method cannot report a direction
#'       (\code{fry} on a single-gene set). The two-sided p-value in
#'       \code{PValue} does not depend on the direction label.}
#'     \item{\code{PValue}}{Two-sided p-value for the set. See the note below.}
#'     \item{\code{Statistic}}{The set-level statistic carrying the direction:
#'       the mean moderated t for \code{camera} (recovered from the two tail
#'       probabilities), the moderated set t for \code{fry}, and the
#'       difference in active proportions for \code{roast}.}
#'     \item{\code{Method}}{The method that produced the row.}
#'   }
#'
#' @section What \code{dsea()} returns, and what it does not:
#'
#' \code{dsea()} returns \strong{p-values}, not per-sample or per-gene
#' measurement values. The set-level tests summarise each gene set into a
#' single statistic and a single p-value; no expression values, fold changes
#' or per-gene statistics are returned, and none should be reconstructed
#' from the output. This boundary is deliberate: the output is a
#' \code{data.frame} that plugs into the same downstream reporting as
#' \code{\link{gsea}} and \code{\link{ora}}, and downstream code must not
#' assume it carries measurements. Use \code{\link{as_dsea}} directly if you
#' want the \code{E} matrix itself.
#'
#' P-values are \emph{not} adjusted for multiple testing. Apply
#' \code{\link[stats]{p.adjust}} on \code{PValue} as appropriate; this is
#' left to the caller because the right family depends on whether the sets
#' form one collection or several.
#'
#' @section Resampling methods and reproducibility:
#'
#' \code{roast} (and \code{fry} passes no rotations) use R's random number
#' generator. The C++ engine draws rotations without reseeding, so results
#' are reproducible exactly when the caller controls the seed, e.g.
#' \code{set.seed(42)} before \code{dsea(...)}. The engine deliberately does
#' \emph{not} call \code{set.seed()} itself, because that would silently
#' ignore the caller's stream.
#'
#' @section Which method to use:
#'
#' \describe{
#'   \item{camera}{Competitive: tests whether the set's genes are more
#'     differentially expressed than the rest of the genome. Fast and robust,
#'     but centres out any genome-wide shift, so it cannot see a global
#'     effect.}
#'   \item{fry}{Self-contained: tests whether the set differs between the two
#'     groups using only the set's own genes. Does not centre across the
#'     genome, so a global shift shows up in every set.}
#'   \item{roast}{Self-contained like \code{fry}, and reports a direction and
#'     a set-level statistic via rotation; useful for small sets where
#'     \code{fry}'s residual degrees of freedom are limiting.}
#' }
#'
#' @examples
#' \donttest{
#' library(DESeq2)
#' library(airway)
#'
#' data(airway)
#' dds <- DESeqDataSet(airway, design = ~ cell + dex)
#' dds <- DESeq(dds)
#'
#' ## Always ask DESeq2 which contrasts exist rather than assuming a name.
#' resultsNames(dds)
#' ct <- "dex_trt_vs_untrt"
#'
#' gene_sets <- list(
#'     TNFa_Signaling = c("TNFAIP3", "NFKB1", "NFKB2", "REL", "RELA"),
#'     OXPHOS         = c("NDUFA1", "NDUFB1", "COX5A", "ATP5F1A", "SDHA")
#' )
#'
#' ## The three methods. Each returns a data.frame of p-values.
#' cam <- dsea(dds, gene_sets, contrast = ct, method = "camera")
#' fry <- dsea(dds, gene_sets, contrast = ct, method = "fry")
#' rst <- dsea(dds, gene_sets, contrast = ct, method = "roast", nrot = 99)
#'
#' ## Multiple-testing correction is left to the caller.
#' cam$FDR <- p.adjust(cam$PValue, method = "BH")
#'
#' ## A bare list works too (backward compatibility).
#' dd   <- as_dsea(dds)
#' bare <- dsea(list(E = dd$E, weights = dd$weights, design = dd$design),
#'              gene_sets, contrast = colnames(dd$design)[ncol(dd$design)],
#'              method = "fry")
#' }
#'
#' ## Only one engine exists; asking for another one is an error, not a
#' ## silent fallback.
#' tryCatch(dsea(list(E = matrix(0, 2, 2), design = matrix(1, 2, 1)),
#'               list(a = "x"), engine = "other"),
#'          error = function(e) conditionMessage(e))
#'
#' @seealso \code{\link{as_dsea}}, \code{\link{gsea}}
#' @importFrom stats p.adjust
#' @importFrom utils getFromNamespace modifyList
#' @export
dsea <- function(x, gene_sets, contrast = NULL,
                 method = c("camera", "roast", "fry"),
                 engine = "native", sort = TRUE, verbose = FALSE, ...) {

    ## ---- engine: "native" only -------------------------------------------
    ## Do NOT remove this check. The alternative engine was deleted on purpose
    ## (commit ca09bf0); silently falling through would change every p-value.
    if (missing(engine) || is.null(engine)) {
        engine <- "native"
    }
    if (!is.character(engine) || length(engine) != 1L || is.na(engine)) {
        stop("'engine' must be a single character string.", call. = FALSE)
    }
    if (!identical(engine, "native")) {
        stop("engine = \"", engine, "\" is not available. ",
             "This package has exactly one statistical engine, the native ",
             "C++ one, so 'engine' accepts only \"native\". The alternative ",
             "engine was removed on purpose: it was a third-party ",
             "dependency the package no longer carries. If you need that ",
             "comparison, install the other package yourself and call its ",
             "set tests directly -- nothing in enrichit will route there.",
             call. = FALSE)
    }

    method <- match.arg(method)
    verbose <- isTRUE(verbose)

    ## ---- normalise the container -----------------------------------------
    dd <- .as_dsea_data(x)
    E <- dd$E
    W <- dd$weights
    design <- dd$design

    if (anyNA(E)) {
        stop("'E' contains NA values; the set tests require a complete ",
             "expression matrix. Drop or impute the affected genes first.",
             call. = FALSE)
    }
    if (any(!is.finite(E))) {
        stop("'E' contains non-finite values; the set tests require finite ",
             "expression values.", call. = FALSE)
    }
    if (!is.null(W) && any(!is.finite(W))) {
        stop("'weights' contains non-finite values.", call. = FALSE)
    }

    n <- ncol(E)
    p <- ncol(design)
    if (n <= p) {
        stop("the design has ", p, " column", if (p == 1L) "" else "s",
             " for only ", n, " sample", if (n == 1L) "" else "s",
             "; there are no residual degrees of freedom.", call. = FALSE)
    }

    ci <- .resolve_contrast_idx(contrast, design)
    contrast_idx <- ci$idx
    if (isTRUE(ci$negate)) {
        design <- .dsea_design(design, contrast_idx, TRUE)
    }

    index <- .build_index(gene_sets, rownames(E))

    if (verbose) {
        message(sprintf("dsea: method=%s, %d genes, %d sets (%d-%d genes), ",
                        method, nrow(E), length(index),
                        min(lengths(index)), max(lengths(index))),
                sprintf("contrast column %d (%s)",
                        contrast_idx,
                        if (is.null(colnames(design))) "<unnamed>"
                        else colnames(design)[contrast_idx]))
    }

    ## ---- dispatch ---------------------------------------------------------
    dots <- list(...)
    known <- switch(method,
                    camera = c("inter.gene.cor", "allow.neg.cor"),
                    roast  = c("nrot", "set.statistic", "midp", "gene.weights"),
                    fry    = c("standardize", "gene.weights"))
    unknown <- setdiff(names(dots), known)
    if (length(unknown)) {
        stop("unused argument", if (length(unknown) == 1L) "" else "s",
             " for method = \"", method, "\": ",
             paste(sprintf("'%s'", unknown), collapse = ", "), call. = FALSE)
    }
    dots <- utils::modifyList(
        list(inter.gene.cor = 0.01, allow.neg.cor = FALSE,
             nrot = 999, set.statistic = "mean", midp = FALSE,
             standardize = "posterior.sd", gene.weights = NULL,
             verbose = verbose),
        dots
    )

    res <- switch(
        method,
        camera = .dsea_camera(E, W, design, contrast_idx, index,
                              dots$inter.gene.cor, dots$allow.neg.cor),
        roast  = {
            ## midp is validated but deliberately NOT forwarded to the engine.
            ## In limma:::mroast.default it only touches the multiple-testing
            ## columns (TwoSidedP2 <- TwoSidedP2 - 1/2/(nrot + 1), then
            ## tab$FDR <- pmax(tab$FDR, pv[,"UpOrDown"])); dsea() returns
            ## unadjusted p-values and no FDR, so it cannot change anything here
            ## -- see the @param midp note. Rejecting a nonsense value is still
            ## the right behaviour rather than silently swallowing it.
            if (!is.logical(dots$midp) || length(dots$midp) != 1L ||
                is.na(dots$midp)) {
                stop("'midp' must be TRUE or FALSE.", call. = FALSE)
            }
            .dsea_roast(E, W, design, contrast_idx, index,
                        dots$set.statistic, dots$gene.weights, dots$nrot)
        },
        fry    = .dsea_fry(E, W, design, contrast_idx, index,
                           dots$gene.weights, dots$standardize)
    )

    res <- .dsea_frame(res, names(index), method)

    if (isTRUE(sort) && nrow(res) > 0L) {
        res <- res[order(res$PValue, -res$NGenes, res$set), , drop = FALSE]
    }
    # Keep the set names as row names, exactly as limma's camera/fry/mroast do.
    # Clearing them here discards the only key that ties a row back to its gene
    # set, which breaks any caller that matches arms by row name (the benchmark
    # parity scripts do exactly that).
    if (nrow(res) > 0L) rownames(res) <- res$set
    res
}

## ---------------------------------------------------------------------------
## Method front-ends: validate, call C++, coerce the result to a frame
## ---------------------------------------------------------------------------

.dsea_camera <- function(E, W, design, contrast_idx, index,
                         inter.gene.cor, allow.neg.cor) {
    ## NA can reach us as logical NA (e.g. straight from a user's
    ## inter.gene.cor = NA) or numeric NA; treat both as "estimate it".
    if (is.na(inter.gene.cor)) {
        est <- TRUE
    } else if (!is.numeric(inter.gene.cor) || length(inter.gene.cor) != 1L) {
        stop("'inter.gene.cor' must be a single number or NA.", call. = FALSE)
    } else if (!is.finite(inter.gene.cor)) {
        stop("'inter.gene.cor' must be NA or a finite number.", call. = FALSE)
    } else if (inter.gene.cor < 0 || inter.gene.cor > 1) {
        ## A correlation outside [0, 1] makes the variance inflation factor
        ## undefined; refuse rather than return a meaningless p-value.
        stop("'inter.gene.cor' must lie in [0, 1] (or be NA to estimate it).",
             call. = FALSE)
    } else est <- FALSE
    if (!is.logical(allow.neg.cor) || length(allow.neg.cor) != 1L) {
        stop("'allow.neg.cor' must be TRUE or FALSE.", call. = FALSE)
    }
    f <- .get_cpp("st_camera_cpp", .cpp_hint)
    ## NA means "estimate the inter-gene correlation per set". The C++ engine
    ## flags that with Rcpp::Nullable, which is null ONLY for R NULL -- an
    ## NA_real_ (or even a logical NA) is a *present* value and would take the
    ## fixed-correlation branch with correlation = NA, making the variance
    ## inflation factor NaN and every p-value NA. So hand over NULL.
    ## (REBUILD_ACCEPTANCE.md, D2.)
    cor_arg <- if (est) NULL else as.numeric(inter.gene.cor)
    f(E, W, design, as.integer(contrast_idx), index, cor_arg,
      isTRUE(allow.neg.cor))
}

.dsea_roast <- function(E, W, design, contrast_idx, index,
                        set.statistic, gene.weights, nrot) {
    set.statistic <- match.arg(set.statistic,
                               c("mean", "floormean", "mean50", "msq"))
    nrot <- .nrot_value(nrot)
    gw <- .gene_weights_value(gene.weights, rownames(E), "roast")
    f <- .get_cpp("st_roast_cpp", .cpp_hint)
    f(E, W, design, as.integer(contrast_idx), index, set.statistic, gw, nrot)
}

.dsea_fry <- function(E, W, design, contrast_idx, index,
                      gene.weights, standardize) {
    standardize <- match.arg(standardize,
                             c("posterior.sd", "none", "residual.sd", "p2"))
    gw <- .gene_weights_value(gene.weights, rownames(E), "fry")
    f <- .get_cpp("st_fry_cpp", .cpp_hint)
    f(E, W, design, as.integer(contrast_idx), index, gw, standardize)
}

## ---------------------------------------------------------------------------
## Result frame
##
## The three engines return the conventional set-test shapes, so the
## mapping is explicit per method rather than sniffed:
##
##   camera: list(NGenes, Correlation, Down, Up, TwoSided, df.camera, df.prior)
##           Down / Up are one-sided tail probabilities; the statistic with the
##           sign of the set's mean log-fold-change is recovered from which
##           tail is the smaller one. TwoSided is the reported p-value.
##   fry:    list(NGenes, t, PValue)          -- t keeps the direction sign.
##   roast:  list(NGenes, Active, PV, var.prior, df.prior)
##           PV is an nsets x 4 matrix: [,1] Down, [,2] Up, [,3] UpOrDown,
##           [,4] Mixed. Sign of t is not returned,
##           so the direction is read from Active.Prop, which roasts the same
##           statistic: Active[,1] is the "Down" proportion and [,2] the "Up".
## ---------------------------------------------------------------------------

.dsea_frame <- function(res, set_names, method) {
    if (!is.list(res) || is.null(res$NGenes)) {
        stop("internal error: the ", method, " engine returned an unexpected ",
             "object (class ", paste(class(res), collapse = "/"), ").",
             call. = FALSE)
    }

    nsets <- length(set_names)
    ngenes <- as.integer(res$NGenes)
    if (length(ngenes) != nsets) {
        stop("internal error: the ", method, " engine returned ", length(ngenes),
             " set sizes for ", nsets, " gene sets.", call. = FALSE)
    }

    if (method == "camera") {
        down <- as.numeric(res$Down)
        up <- as.numeric(res$Up)
        pv <- as.numeric(res$TwoSided)
        ## camera has no signed statistic; the direction is the tail
        ## that carries the smaller probability (equivalently, the sign of the
        ## mean moderated t among the set's genes).
        stat <- down - up

    } else if (method == "fry") {
        stat <- as.numeric(res$t)
        pv <- as.numeric(res$PValue)

    } else {  # roast
        pvmat <- as.matrix(res$PV)
        act <- as.matrix(res$Active)
        if (nrow(pvmat) != nsets || ncol(pvmat) < 3L) {
            stop("internal error: the roast engine returned a ",
                 nrow(pvmat), " x ", ncol(pvmat), " p-value matrix for ",
                 nsets, " gene sets.", call. = FALSE)
        }
        ## Column 3 is "UpOrDown"; with the default set.statistic = "mean" the
        ## engine sets it to min(count[down], count[up]), i.e. exactly the
        ## two-sided p-value.
        pv <- pvmat[, 3L]
        stat <- if (ncol(act) >= 2L) act[, 2L] - act[, 1L] else NA_real_
        ## Direction follows limma::mroast's rule, NOT the sign of `stat`:
        ## limma sets "Up" iff the Up rotation p-value is smaller than the
        ## Down one (limma:::mroast.default, `Up <- pv[,"Up"] < pv[,"Down"]`).
        ## The active-proportion difference and the two one-sided rotation
        ## tails can disagree in sign for a near-null set -- that is where the
        ## two rules part company -- and a drop-in engine must label sets the
        ## same way the reference implementation does. `Statistic` keeps
        ## reporting the proportion difference, so nothing is lost:
        ## Direction = which tail won, Statistic = how much of the set moved.
        up_tail <- if (ncol(pvmat) >= 2L) pvmat[, "Up"] else NA_real_
        down_tail <- if (ncol(pvmat) >= 1L) pvmat[, "Down"] else NA_real_
        ## NA left as NA: a set with no usable tails has no direction.
        dir_rule <- ifelse(is.na(up_tail) | is.na(down_tail), NA_character_,
                           ifelse(up_tail < down_tail, "Up", "Down"))
    }

    if (length(pv) != nsets || length(stat) != nsets) {
        stop("internal error: the ", method, " engine returned ",
             length(pv), " p-values and ", length(stat),
             " statistics for ", nsets, " gene sets.", call. = FALSE)
    }

    bad <- !is.na(pv) & (pv < 0 | pv > 1)
    if (any(bad)) {
        stop("internal error: the ", method, " engine returned p-values ",
             "outside [0, 1] for set", if (sum(bad) == 1L) "" else "s", " ",
             paste(set_names[bad], collapse = ", "), ".", call. = FALSE)
    }

    ## A set whose statistic is unavailable has no direction to report; label
    ## it NA rather than guessing "Down". roast overrides the sign-of-statistic
    ## rule with limma's which-tail-won rule (see dir_rule above).
    if (method == "roast") {
        ## Which one-sided tail won -- limma's rule, already computed as
        ## dir_rule. Sets with no usable tails keep NA.
        direction <- dir_rule
    } else {
        direction <- ifelse(is.na(stat), NA_character_,
                            ifelse(stat > 0, "Up", "Down"))
    }

    data.frame(
        set = set_names,
        NGenes = ngenes,
        Direction = direction,
        PValue = as.numeric(pv),
        Statistic = as.numeric(stat),
        Method = method,
        stringsAsFactors = FALSE
    )
}
