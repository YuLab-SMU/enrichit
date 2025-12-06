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

    qobj <- tryCatch(qvalue::qvalue(pvals, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)

    if (inherits(qobj, "qvalue")) {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }
    return(qvalues)
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
        i <- match(term, gsid2name)
        j <- !is.na(i)
        res <- term
        res[j] <- gsid2name[i[j]]
    } else {
        res <- as.character(term)
    }

    names(res) <- term
    return(res) 
}

TERMID2EXTID <- function(term, gson) {
    if (inherits(gson, "GSON")) {
        gsid2gene <- gson@gsid2gene
        gsid2gene <- gsid2gene[gsid2gene %in% term, ]
        res <- split(gsid2gene, gsid2gene)
        return(res)
    } else if (inherits(gson, "environment")) {
        PATHID2EXTID <- get("PATHID2EXTID", envir = gson)
        res <- PATHID2EXTID[term]
        return(res)
    } else {
        stop("gson not supported")
    }
}
