#' interal method for enrichment analysis
#'
#' using the hypergeometric model
#' @title ora-gson
#' @param gene a vector of entrez gene id.
#' @param pvalueCutoff Cutoff value of pvalue.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes, default is the intersection of the 'universe' with genes that have annotations. 
#' Users can set `options(enrichment_force_universe = TRUE)` to force the 'universe' untouched.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param qvalueCutoff cutoff of qvalue
#' @param gson ontology information
#' @return  A `enrichResult` instance.
#' @importClassesFrom methods data.frame
#' @importFrom methods new
#' @importFrom stats p.adjust
#' @keywords manip
#' @author Guangchuang Yu <https://yulab-smu.top>
#' @export
ora_gson <- function(gene,
                              pvalueCutoff,
                              pAdjustMethod="BH",
                              universe = NULL,
                              minGSSize=10,
                              maxGSSize=500,
                              qvalueCutoff=0.2,
                              gson){

    if (!inherits(gson, "GSON")) {
        stop("gson should be a GSON object")
    }

    ## query external ID to Term ID
    gene <- as.character(unique(gene))
    
    # Extract gene sets from GSON
    gsid2gene <- gson@gsid2gene
    
    # ID Match Check
    if (!check_gene_id(gene, gsid2gene)) {
        return(NULL)
    }

    # Handle universe
    extID <- unique(gsid2gene$gene)
    if (missing(universe))
        universe <- NULL
    if(!is.null(universe)) {
        if (is.character(universe)) {
            force_universe <- getOption("enrichment_force_universe", FALSE)
            if (force_universe) {
                extID <- universe
            } else {
                extID <- intersect(extID, universe)
            }
        } else {
            ## https://github.com/YuLab-SMU/clusterProfiler/issues/217
            message("`universe` is not in character and will be ignored...")
        }
    }

    # Prepare Gene Sets
    geneSets <- split(gsid2gene$gene, gsid2gene$gsid)
    # Intersect with universe
    geneSets <- lapply(geneSets, intersect, extID)
    
    # Filter by size
    idx <- get_geneSet_index(geneSets, minGSSize, maxGSSize)
    if (sum(idx) == 0) {
        msg <- paste("No gene sets have size between", minGSSize, "and", maxGSSize, "...")
        message(msg)
        message("--> return NULL...")
        return (NULL)
    }
    geneSets <- geneSets[idx]

    ora_res <- ora(gene, geneSets, universe = extID)

    if (is.null(ora_res) || nrow(ora_res) == 0) {
        return(NULL)
    }

    # Calculate p.adjust
    ora_res$p.adjust <- p.adjust(ora_res$pvalue, method=pAdjustMethod)

    # Calculate qvalue
    ora_res$qvalue <- calculate_qvalue(ora_res$pvalue)

    # Calculate zScore
    # Need k, M, n, N    
    N <- length(extID)
    n <- length(intersect(gene, extID))
    k <- ora_res$Count
    
    # M = size of gene set (in universe)
    # Map ID to M
    # ora_res$ID should match names(geneSets)
    M <- sapply(geneSets[ora_res$ID], length)
    
    mu <- M * n / N
    sigma <- mu * (N - n) * (N - M) / N / (N - 1)
    zScore <- (k - mu) / sqrt(sigma)
    ora_res$zScore <- zScore

    # Add Description
    gsid2name <- gson@gsid2name
    if (!is.null(gsid2name) && "ID" %in% names(ora_res)) {
        description <- gsid2name$name[match(ora_res$ID, gsid2name$gsid)]
        na_idx <- is.na(description)
        description[na_idx] <- ora_res$ID[na_idx]
        ora_res$Description <- description
    } else {
        if (!"Description" %in% names(ora_res)) {
            ora_res$Description <- ora_res$ID
        }
    }

    # Reorder columns
    expected_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    other_cols <- setdiff(names(ora_res), expected_cols)
    ora_res <- ora_res[, c(expected_cols, other_cols)]
    
    # Sort by pvalue
    ora_res <- ora_res[order(ora_res$pvalue), ]
    
    # Set row names
    row.names(ora_res) <- ora_res$ID

    x <- new("enrichResult",
             result         = ora_res,
             pvalueCutoff   = pvalueCutoff,
             pAdjustMethod  = pAdjustMethod,
             qvalueCutoff   = qvalueCutoff,
             gene           = as.character(gene),
             universe       = extID,
             geneSets       = geneSets,
             organism       = if (!is.null(gson@species)) gson@species else "UNKNOWN",
             keytype        = if (!is.null(gson@keytype)) gson@keytype else "UNKNOWN",
             ontology       = if (!is.null(gson@gsname)) gsub(".*;", "", gson@gsname) else "UNKNOWN",
             readable       = FALSE
             )
             
    return(x)
}


get_enriched <- function(object) {

    Over <- object@result

    pvalueCutoff <- object@pvalueCutoff
    if (length(pvalueCutoff) != 0) {
        ## if groupGO result, numeric(0)
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
        i <- match(term, gsid2name$gsid)
        j <- !is.na(i)
        res <- term
        res[j] <- gsid2name$name[i[j]]
    } else {
        res <- as.character(term)
    }

    names(res) <- term
    return(res) 
}

get_geneSet_index <- function(geneSets, minGSSize, maxGSSize) {
    if (is.na(minGSSize) || is.null(minGSSize))
        minGSSize <- 1
    if (is.na(maxGSSize) || is.null(maxGSSize))
        maxGSSize <- Inf #.Machine$integer.max

    ## index of geneSets in used.
    ## logical
    geneSet_size <- sapply(geneSets, length)
    idx <-  minGSSize <= geneSet_size & geneSet_size <= maxGSSize
    return(idx)
}

TERMID2EXTID <- function(term, gson) {
    if (inherits(gson, "GSON")) {
        gsid2gene <- gson@gsid2gene
        gsid2gene <- gsid2gene[gsid2gene$gsid %in% term, ]
        res <- split(gsid2gene$gene, gsid2gene$gsid)
        return(res)
    } else if (inherits(gson, "environment")) {
        PATHID2EXTID <- get("PATHID2EXTID", envir = gson)
        res <- PATHID2EXTID[term]
        return(res)
    } else {
        stop("gson not supported")
    }
}


calculate_qvalue <- function(pvals) {
    if (length(pvals) == 0)
        return(numeric(0))

    qobj <- tryCatch(qvalue::qvalue(pvals, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)

    # if (class(qobj) == "qvalue") {
    if (inherits(qobj, "qvalue")) {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }
    return(qvalues)
}

# https://github.com/YuLab-SMU/ReactomePA/issues/43
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
