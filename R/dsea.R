## Set-level competitive gene-set tests on DESeq2 / EList input.
##
## `dsea()` runs camera / roast / fry. The default engine is a
## self-contained C++ implementation that mirrors limma's algorithms
## (contrast-column QR, eBayes-style variance moderation, rotation or
## SVD-based set statistics); `engine = "limma"` dispatches to limma's
## own functions and supports their full argument sets. limma remains a
## lazy Suggests dependency and serves as the reference backend.

#' Competitive gene-set tests on a DESeqDataSet or EList
#'
#' Run competitive gene-set tests -- \code{camera}, \code{roast} or
#' \code{fry} -- directly on a fitted
#' \code{\link[DESeq2:DESeqDataSet-class]{DESeqDataSet}}, or on any
#' \code{\link[limma:EList-class]{EList}} such as the one produced by
#' \code{\link{as_elist}}.
#'
#' @details Competitive tests ask whether the genes of a set shift
#'   consistently, relative to all other genes, under a contrast of the
#'   experimental design. This is the design-aware, no-cut-off alternative
#'   to \code{\link{ora}} (which requires a significance cut-off) and to
#'   \code{\link{gsea}} (which reduces the design to a single ranked list):
#'   any contrast expressible in the DESeq2 design formula can be tested.
#'
#'   By default (\code{engine = "native"}) the tests run in a
#'   self-contained C++ engine that mirrors limma's algorithms, including
#'   the eBayes-style variance moderation and (for \code{roast}) the
#'   rotation scheme, which consumes R's random-number stream in the same
#'   order as limma -- so the same \code{set.seed()} reproduces limma's
#'   p-values exactly. \code{engine = "limma"} dispatches to
#'   \code{\link[limma:camera]{limma::camera}},
#'   \code{\link[limma:roast]{limma::roast}} /
#'   \code{\link[limma:roast]{limma::mroast}} or
#'   \code{\link[limma:fry]{limma::fry}} and supports their full argument
#'   sets through \code{...}; it requires the limma package.
#'
#'   For a \code{DESeqDataSet}, the object is converted with
#'   \code{\link{as_elist}}, whose per-observation precision weights are
#'   derived from DESeq2's own fitted means and dispersions rather than
#'   from a re-learned mean-variance trend.
#'
#'   Arguments supported by the native engine, with limma-identical
#'   defaults:
#'   \itemize{
#'     \item \code{camera}: \code{inter.gene.cor = 0.01} (pass \code{NA}
#'       to estimate it from the residuals), \code{allow.neg.cor = FALSE},
#'       \code{sort = TRUE}, \code{directional = TRUE}.
#'     \item \code{roast}: \code{set.statistic = "mean"},
#'       \code{nrot = 1999}, \code{gene.weights}, \code{midp = TRUE},
#'       \code{sort = "directional"}. The result is \emph{mroast}-style
#'       (limma's \code{roast} also delegates to \code{mroast} when given
#'       multiple sets); the \code{PValue} column always holds the
#'       unadjusted rotation p-value and \code{FDR}/\code{FDR.Mixed} use
#'       Benjamini-Hochberg (limma hardcodes the same method).
#'     \item \code{fry}: \code{gene.weights},
#'       \code{standardize = "posterior.sd"}, \code{sort = "directional"}.
#'   }
#'   Any other limma argument (e.g. \code{use.ranks}, \code{trend.var},
#'   \code{robust}, \code{array.weights}, \code{block}) is rejected by the
#'   native engine with a pointer to \code{engine = "limma"}.
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
#'   \code{ncol(design)} (\code{engine = "limma"} only), a DESeq2-style
#'   three-token contrast \code{c(factor, numerator, denominator)} which
#'   is mapped to the corresponding design column, or \code{NULL}
#'   (default) to test the last column of the design, matching limma's
#'   and DESeq2's defaults.
#' @param method character, the test to run: \code{"camera"} (default),
#'   \code{"roast"} or \code{"fry"}.
#' @param engine character, \code{"native"} (default) for the built-in
#'   C++ engine, or \code{"limma"} to delegate to the limma package.
#' @param ... additional arguments forwarded to the underlying engine;
#'   see Details for what each method supports.
#' @return \code{engine = "limma"} returns the result object of the
#'   selected limma test, unchanged. \code{engine = "native"} returns a
#'   data frame with one row per gene set, sorted as limma does:
#'   \code{camera} gives \code{NGenes}, \code{Direction}, \code{PValue}
#'   (plus \code{Correlation} when the inter-gene correlation is
#'   estimated, and \code{FDR} for multiple sets); \code{roast} gives the
#'   \emph{mroast} columns \code{NGenes}, \code{PropDown}, \code{PropUp},
#'   \code{Direction}, \code{PValue}, \code{FDR}, \code{PValue.Mixed},
#'   \code{FDR.Mixed}; \code{fry} gives \code{NGenes}, \code{Direction},
#'   \code{PValue}, \code{PValue.Mixed} (plus \code{FDR} and
#'   \code{FDR.Mixed} for multiple sets).
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
#' res <- dsea(dds, gene_sets, contrast = c("condition", "B", "A"))
#' ## the limma backend with its full argument set:
#' res2 <- dsea(dds, gene_sets, engine = "limma", trend.var = TRUE)
#' }
#' @export
dsea <- function(x, gene_sets, contrast = NULL,
                     method = c("camera", "roast", "fry"),
                     engine = c("native", "limma"), ...) {
    method <- match.arg(method)
    engine <- match.arg(engine)

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
    design <- as.matrix(design)

    ids <- rownames(elist$E)
    idx <- .validate_gene_sets(gene_sets, ids)

    if (engine == "limma") {
        rlang::check_installed("limma", "for dsea(engine = \"limma\").")
        ct <- .resolve_contrast(contrast, design)
        return(switch(method,
            camera = limma::camera(elist, idx, design = design, contrast = ct, ...),
            roast = limma::roast(elist, idx, design = design, contrast = ct, ...),
            fry = limma::fry(elist, idx, design = design, contrast = ct, ...)
        ))
    }

    ci <- .resolve_contrast_idx(contrast, design)
    W <- .expand_weights(elist$weights, nrow(elist$E), ncol(elist$E))
    idx0 <- lapply(idx, function(s) {
        p <- match(s, ids)
        sort(unique(p[!is.na(p)])) - 1L
    })
    dots <- list(...)
    if (length(dots) && (is.null(names(dots)) || "" %in% names(dots))) {
        stop("all ... arguments for engine = \"native\" must be named.",
             call. = FALSE)
    }

    switch(method,
        camera = .camera_native(elist$E, W, design, ci, idx0, dots,
            "inter.gene.cor" %in% names(match.call(expand.dots = TRUE))[-1L]),
        roast = .roast_native(elist$E, W, design, ci, idx0, dots),
        fry = .fry_native(elist$E, W, design, ci, idx0, dots)
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

.resolve_contrast_idx <- function(contrast, design) {
    if (is.null(contrast)) {
        return(ncol(design))
    }
    if (is.character(contrast) && length(contrast) == 3L) {
        col <- paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3])
        j <- match(col, colnames(design))
        if (is.na(j)) {
            stop(
                "contrast c(", paste(contrast, collapse = ", "),
                ") maps to design column '", col, "', which does not exist; ",
                "available columns: ", paste(colnames(design), collapse = ", "),
                ".",
                call. = FALSE
            )
        }
        return(j)
    }
    if (length(contrast) == 1L && is.numeric(contrast)) {
        return(as.integer(contrast))
    }
    if (length(contrast) == 1L && is.character(contrast)) {
        j <- match(contrast, colnames(design))
        if (is.na(j)) {
            stop(
                "contrast '", contrast, "' not found in the design ",
                "column names: ", paste(colnames(design), collapse = ", "),
                ".",
                call. = FALSE
            )
        }
        return(j)
    }
    stop(
        "engine = \"native\" only supports contrasts that name a single ",
        "design column; pass a numeric contrast vector with ",
        "engine = \"limma\".",
        call. = FALSE
    )
}

.expand_weights <- function(w, G, n) {
    if (is.null(w)) {
        return(NULL)
    }
    if (is.matrix(w)) {
        if (nrow(w) != G || ncol(w) != n) {
            stop("weights must have the same dimensions as the data.",
                 call. = FALSE)
        }
        return(as.matrix(w))
    }
    if (is.numeric(w) && length(w) == G) {
        return(matrix(w, G, n))
    }
    if (is.numeric(w) && length(w) == n) {
        return(matrix(w, G, n, byrow = TRUE))
    }
    stop(
        "weights must be a genes x samples matrix or a vector of length ",
        "nrow or ncol of the data.",
        call. = FALSE
    )
}

.check_dots_native <- function(dots, allowed, method) {
    bad <- setdiff(names(dots), allowed)
    if (length(bad)) {
        stop(
            "argument(s) ", paste0("'", bad, "'", collapse = ", "),
            " not supported by dsea(method = \"", method,
            "\", engine = \"native\"); use engine = \"limma\" instead.",
            call. = FALSE
        )
    }
}

.coerce_sort <- function(sort, default) {
    if (is.logical(sort)) {
        sort <- if (sort) "directional" else "none"
    }
    if (is.null(sort)) {
        sort <- default
    }
    match.arg(sort, c("directional", "mixed", "none"))
}

.check_gene_weights <- function(gw, G) {
    if (is.null(gw)) {
        return(NULL)
    }
    if (!is.numeric(gw)) {
        stop("gene.weights must be numeric.", call. = FALSE)
    }
    if (length(gw) != G) {
        stop("gene.weights vector should be of length nrow(E).",
             call. = FALSE)
    }
    gw
}

.camera_native <- function(E, W, design, ci, idx0, dots, igc_supplied) {
    .check_dots_native(
        dots, c("inter.gene.cor", "allow.neg.cor", "sort", "directional"),
        "camera"
    )
    if (!is.null(dots$directional) && !isTRUE(dots$directional)) {
        stop(
            "directional = FALSE requires the rank-based test, which ",
            "engine = \"native\" does not implement; use ",
            "engine = \"limma\".",
            call. = FALSE
        )
    }
    igc <- dots$inter.gene.cor
    fixed <- TRUE
    if (!igc_supplied) {
        igc <- 0.01
    } else if (is.null(igc) || is.na(igc)) {
        fixed <- FALSE
        igc <- NULL
    }
    sortv <- dots$sort
    if (is.null(sortv)) sortv <- TRUE
    if (!is.logical(sortv) || length(sortv) != 1L || is.na(sortv)) {
        stop("sort must be TRUE or FALSE.", call. = FALSE)
    }

    res <- st_camera_cpp(
        E, W, design, ci, idx0, igc, isTRUE(dots$allow.neg.cor)
    )
    nsets <- length(idx0)
    tab <- data.frame(NGenes = as.integer(res$ngenes))
    if (!fixed) {
        tab$Correlation <- res$correlation
    }
    tab$Direction <- ifelse(res$down < res$up, "Down", "Up")
    tab$PValue <- res$twosided
    if (nsets > 1) {
        tab$FDR <- p.adjust(tab$PValue, method = "BH")
    }
    rownames(tab) <- names(idx0)
    if (sortv && nsets > 1) {
        tab <- tab[order(tab$PValue), , drop = FALSE]
    }
    tab
}

.roast_native <- function(E, W, design, ci, idx0, dots) {
    .check_dots_native(
        dots,
        c("set.statistic", "nrot", "gene.weights", "sort", "midp"),
        "roast"
    )
    ss <- dots$set.statistic
    ss <- if (is.null(ss)) {
        "mean"
    } else {
        match.arg(ss, c("mean", "floormean", "mean50", "msq"))
    }
    nrot <- dots$nrot
    if (is.null(nrot)) {
        nrot <- 1999L
    } else {
        nrot <- as.integer(nrot)
        if (is.na(nrot) || nrot < 1L) {
            stop("nrot must be a positive integer.", call. = FALSE)
        }
    }
    gw <- .check_gene_weights(dots$gene.weights, nrow(E))
    midp <- if (is.null(dots$midp)) TRUE else isTRUE(dots$midp)
    sortv <- .coerce_sort(dots$sort, "directional")

    res <- st_roast_cpp(E, W, design, ci, idx0, ss, gw, nrot)
    pv <- res$pv
    act <- res$active
    nsets <- length(idx0)
    direction <- rep_len("Down", nsets)
    direction[pv[, 2] < pv[, 1]] <- "Up"
    pvalue <- pv[, 3]
    pmixed <- pv[, 4]
    if (midp) {
        tab_fdr <- pmax(
            p.adjust(pvalue - 0.5 / (nrot + 1), method = "BH"), pvalue
        )
        tab_fdr_m <- pmax(
            p.adjust(pmixed - 0.5 / (nrot + 1), method = "BH"), pmixed
        )
    } else {
        tab_fdr <- p.adjust(pvalue, method = "BH")
        tab_fdr_m <- p.adjust(pmixed, method = "BH")
    }
    tab <- data.frame(
        NGenes = as.integer(res$ngenes),
        PropDown = act[, 1],
        PropUp = act[, 2],
        Direction = direction,
        PValue = pvalue,
        FDR = tab_fdr,
        PValue.Mixed = pmixed,
        FDR.Mixed = tab_fdr_m,
        row.names = names(idx0)
    )
    if (midp) {
        tab$FDR <- pmax(tab$FDR, pvalue)
        tab$FDR.Mixed <- pmax(tab$FDR.Mixed, pmixed)
    }
    if (sortv == "directional") {
        prop <- pmax(tab$PropUp, tab$PropDown)
        o <- order(tab$PValue, -prop, -tab$NGenes, tab$PValue.Mixed)
    } else if (sortv == "mixed") {
        prop <- tab$PropUp + tab$PropDown
        o <- order(tab$PValue.Mixed, -prop, -tab$NGenes, tab$PValue)
    } else {
        o <- seq_len(nsets)
    }
    tab[o, , drop = FALSE]
}

.fry_native <- function(E, W, design, ci, idx0, dots) {
    .check_dots_native(dots, c("gene.weights", "standardize", "sort"), "fry")
    gw <- .check_gene_weights(dots$gene.weights, nrow(E))
    std <- dots$standardize
    std <- if (is.null(std)) {
        "posterior.sd"
    } else {
        match.arg(std, c("none", "residual.sd", "posterior.sd", "p2"))
    }
    sortv <- .coerce_sort(dots$sort, "directional")

    res <- st_fry_cpp(E, W, design, ci, idx0, gw, std)
    nsets <- length(idx0)
    ng <- as.integer(res$ngenes)
    tstat <- res$t.stat
    pvalue <- 2 * stats::pt(-abs(tstat), df = res$df.residual)
    pmixed <- res$pvalue.mixed
    pmixed[ng == 1] <- pvalue[ng == 1]
    tab <- data.frame(
        NGenes = ng,
        Direction = ifelse(tstat < 0, "Down", "Up"),
        PValue = pvalue
    )
    if (nsets > 1) {
        tab$FDR <- p.adjust(pvalue, method = "BH")
        tab$PValue.Mixed <- pmixed
        tab$FDR.Mixed <- p.adjust(pmixed, method = "BH")
    } else {
        tab$PValue.Mixed <- pmixed
    }
    rownames(tab) <- names(idx0)
    if (sortv == "directional") {
        o <- order(tab$PValue, -tab$NGenes, tab$PValue.Mixed)
    } else if (sortv == "mixed") {
        o <- order(tab$PValue.Mixed, -tab$NGenes, tab$PValue)
    } else {
        o <- seq_len(nsets)
    }
    tab[o, , drop = FALSE]
}
