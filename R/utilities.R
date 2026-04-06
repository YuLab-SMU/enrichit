#' @importFrom yulab.utils yulab_msg
check_gene_id <- function(gene, gsid2gene) {
    if (!any(gene %in% gsid2gene$gene)) {
        yulab_msg("--> No gene can be mapped....")
        sg <- unique(gsid2gene$gene[1:100])
        sg <- sample(sg, min(length(sg), 6))
        yulab_msg("--> Expected input gene ID: ", paste0(sg, collapse=','))
        yulab_msg("--> return NULL...")
        return(FALSE)
    }
    return(TRUE)
}

calculate_qvalue <- function(pvals) {
    if (length(pvals) == 0)
        return(numeric(0))

    # Check if qvalue package is available
    if (!requireNamespace("qvalue", quietly = TRUE)) {
        warning("qvalue package not installed. qvalue column will be NA.")
        return(rep(NA_real_, length(pvals)))
    }

    qv <- rep(NA_real_, length(pvals))
    ok <- is.finite(pvals) & !is.na(pvals) & pvals >= 0 & pvals <= 1
    if (any(!ok)) {
        warning("Invalid p-values detected (NA, non-finite, <0, or >1). qvalue will be computed on valid p-values only.")
    }

    p <- pvals[ok]
    if (length(p) == 0) {
        return(qv)
    }

    if (all(p == 1)) {
        qv[ok] <- 1
        return(qv)
    }

    if (all(p == 0)) {
        qv[ok] <- 0
        return(qv)
    }

    qobj <- NULL
    last_err <- NULL
    try_qvalue <- function(expr) {
        tryCatch(expr, error = function(e) {
            last_err <<- e
            NULL
        })
    }

    qobj <- try_qvalue(qvalue::qvalue(p))
    if (!inherits(qobj, "qvalue")) {
        qobj <- try_qvalue(qvalue::qvalue(p, lambda = 0.05, pi0.method = "bootstrap"))
    }
    if (!inherits(qobj, "qvalue")) {
        qobj <- try_qvalue(qvalue::qvalue(p, lambda = seq(0, 0.9, 0.05), pi0.method = "bootstrap"))
    }
    if (!inherits(qobj, "qvalue")) {
        qobj <- try_qvalue(qvalue::qvalue(p, lambda = seq(0, 0.95, 0.05), pi0.method = "bootstrap"))
    }

    if (inherits(qobj, "qvalue")) {
        qv[ok] <- qobj$qvalues
        return(qv)
    }

    if (inherits(last_err, "error")) {
        warning("qvalue::qvalue() failed, returning NA for qvalue. Error: ", conditionMessage(last_err))
    } else {
        warning("qvalue::qvalue() failed, returning NA for qvalue. This may be due to p-value distribution.")
    }
    return(qv)
}

get_geneSet_index <- function(geneSets, minGSSize, maxGSSize) {
    if (is.na(minGSSize) || is.null(minGSSize))
        minGSSize <- 1
    if (is.na(maxGSSize) || is.null(maxGSSize))
        maxGSSize <- Inf 

    ## index of geneSets in used.
    ## logical
    geneSet_size <- sapply(geneSets, length)
    idx <-  minGSSize <= geneSet_size & geneSet_size <= maxGSSize
    return(idx)
}

validate_gene_sets <- function(gene_sets) {
    if (!is.list(gene_sets) || is.null(names(gene_sets))) {
        stop("gene_sets must be a named list")
    }
    gene_sets <- lapply(gene_sets, function(x) {
        if (!is.character(x)) {
            stop("Each element in gene_sets must be a character vector")
        }
        unique(x)
    })
    return(gene_sets)
}

gsea_leading_edge_details <- function(geneList, geneSet, exponent, scoreType) {
    genes <- names(geneList)
    N <- length(geneList)
    geneSet <- unique(intersect(as.character(geneSet), genes))
    if (length(geneSet) == 0 || N == 0) {
        return(list(rank = 0L, leading_edge = "tags=0%, list=0%, signal=0%", core_enrichment = ""))
    }

    in_set <- !is.na(match(genes, geneSet))
    N_H <- sum(in_set)
    if (N_H == 0) {
        return(list(rank = 0L, leading_edge = "tags=0%, list=0%, signal=0%", core_enrichment = ""))
    }

    weights <- abs(geneList)
    if (!isTRUE(all.equal(exponent, 1.0))) {
        weights <- weights^exponent
    }
    N_R <- sum(weights[in_set])
    N_miss <- N - N_H

    hit_inc <- if (N_R == 0) rep(0, N) else (weights * in_set) / N_R
    miss_inc <- if (N_miss == 0) rep(0, N) else (!in_set) / N_miss
    running <- cumsum(hit_inc) - cumsum(miss_inc)

    if (scoreType == "pos") {
        peak_idx <- which.max(running)
        es <- running[peak_idx]
    } else if (scoreType == "neg") {
        peak_idx <- which.min(running)
        es <- running[peak_idx]
    } else {
        max_es <- max(running)
        min_es <- min(running)
        if (abs(max_es) >= abs(min_es)) {
            peak_idx <- which.max(running)
            es <- max_es
        } else {
            peak_idx <- which.min(running)
            es <- min_es
        }
    }

    if (is.na(peak_idx) || length(peak_idx) == 0) {
        return(list(rank = 0L, leading_edge = "tags=0%, list=0%, signal=0%", core_enrichment = ""))
    }
    peak_idx <- peak_idx[[1]]

    if (es >= 0) {
        rank <- as.integer(peak_idx)
        ledge_range <- seq_len(peak_idx)
        list_frac <- peak_idx / N
    } else {
        rank <- as.integer(N - peak_idx + 1)
        ledge_range <- peak_idx:N
        list_frac <- (N - peak_idx + 1) / N
    }

    core_genes <- genes[ledge_range][in_set[ledge_range]]
    hits_in_ledge <- length(core_genes)

    tags_frac <- hits_in_ledge / N_H
    signal_frac <- tags_frac * (1 - list_frac) * (N / max(N_miss, 1))
    leading_edge <- paste0(
        "tags=", round(tags_frac * 100), "%",
        ", list=", round(list_frac * 100), "%",
        ", signal=", round(signal_frac * 100), "%"
    )

    list(
        rank = rank,
        leading_edge = leading_edge,
        core_enrichment = paste0(core_genes, collapse = "/")
    )
}

#' Common parameters for enrichit functions
#'
#' @param geneList A named numeric vector of gene statistics (e.g., log fold change), ranked in descending order.
#' @param gene_sets A named list of gene sets. Each element is a character vector of genes.
#' @param nPerm Number of permutations for p-value calculation (default: 1000).
#' @param exponent Weighting exponent for enrichment score (default: 1.0).
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param pvalueCutoff P-value cutoff.
#' @param pAdjustMethod P-value adjustment method (e.g., "BH").
#' @param verbose Logical. Print progress messages.
#' @param gson A GSON object containing gene set information.
#' @param method Permutation method.
#' @param adaptive Logical. Use adaptive permutation.
#' @param minPerm Minimum permutations for adaptive mode.
#' @param maxPerm Maximum permutations for adaptive mode.
#' @param pvalThreshold P-value threshold for early stopping.
#' @name enrichit_params
NULL

get_enriched <- function(object) {

    Over <- object@result

    pvalueCutoff <- object@pvalueCutoff
    if (length(pvalueCutoff) != 0) {
        Over <- Over[ Over$pvalue <= pvalueCutoff, ]
        Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
    }

    qvalueCutoff <- object@qvalueCutoff
    if (length(qvalueCutoff) != 0) {
        if (! any(is.na(Over$qvalue))) {
            if (length(qvalueCutoff) > 0)
                Over <- Over[ Over$qvalue <= qvalueCutoff, ]
        }
    }

    object@result <- Over
    return(object)
}

TERM2NAME <- function(term, gson) {
    if (inherits(gson, "environment")) { 
        PATHID2NAME <- get("PATHID2NAME", envir = gson)
        #if (is.null(PATHID2NAME) || is.na(PATHID2NAME)) {
        if (is.null(PATHID2NAME) || all(is.na(PATHID2NAME))) {
            return(as.character(term))
        }
        res <- PATHID2NAME[term]
        i <-  is.na(res)
        res[i] <- term[i]
    } else if (inherits(gson, "GSON")) {
        gsid2name <- gson@gsid2name
        res <- as.character(term)
        if (!is.null(gsid2name) && nrow(gsid2name) > 0L) {
            i <- match(as.character(term), as.character(gsid2name$gsid))
            j <- !is.na(i)
            res[j] <- as.character(gsid2name$name[i[j]])
        }
    } else {
        res <- as.character(term)
    }

    names(res) <- term
    return(res) 
}

TERMID2EXTID <- function(term, gson) {
    if (inherits(gson, "GSON")) {
        gsid2gene <- gson@gsid2gene
        res <- split(
            as.character(gsid2gene$gene),
            as.character(gsid2gene$gsid)
        )
        res <- res[as.character(term)]
        names(res) <- as.character(term)
        i <- vapply(res, is.null, logical(1))
        res[i] <- list(character(0))
        return(res)
    } else if (inherits(gson, "environment")) {
        PATHID2EXTID <- get("PATHID2EXTID", envir = gson)
        res <- PATHID2EXTID[term]
        return(res)
    } else {
        stop("gson not supported")
    }
}
