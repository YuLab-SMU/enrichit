#' Gene Set Enrichment Analysis (GSEA)
#'
#' Perform Gene Set Enrichment Analysis (GSEA) using a ranked gene list.
#'
#' @param geneList A named numeric vector of gene statistics (e.g., log fold change), ranked in descending order.
#' @param gene_sets A named list of gene sets. Each element is a character vector of genes.
#' @param nPerm Number of permutations for p-value calculation (default: 1000). Only used when adaptive=FALSE.
#' @param exponent Weighting exponent for enrichment score (default: 1.0).
#' @param method Permutation method: "sample" (default) for random gene set sampling (faster, similar to fgsea), 
#' or "permute" for label permutation (slower, standard GSEA).
#' @param adaptive Logical. If TRUE, use adaptive early-stopping permutation for more accurate p-values 
#' on significant gene sets. Default: FALSE for backward compatibility.
#' @param minPerm Minimum number of permutations for adaptive mode (default: 1000).
#' @param maxPerm Maximum number of permutations for adaptive mode (default: 100000).
#' @param pvalThreshold P-value threshold for early stopping in adaptive mode (default: 0.1). 
#' Gene sets with p-value > threshold will stop after minPerm permutations.
#'
#' @return A data.frame with columns:
#' - **ID**: Gene set name
#' - **enrichmentScore**: Enrichment Score
#' - **NES**: Normalized Enrichment Score
#' - **pvalue**: Empirical p-value from permutation test
#' - **setSize**: Size of the gene set (number of genes found in geneList)
#' - **nPerm**: (adaptive mode only) Actual number of permutations used
#' - **rank**: Rank at which the maximum enrichment score is attained
#' - **leading_edge**: Leading edge statistics (tags, list, signal)
#' - **core_enrichment**: Genes in the leading edge, separated by '/'
#'
#' @examples
#' # Example data
#' stats <- rnorm(1000)
#' names(stats) <- paste0("Gene", 1:1000)
#' stats <- sort(stats, decreasing = TRUE)
#' 
#' gs1 <- paste0("Gene", 1:50)
#' gs2 <- paste0("Gene", 500:550)
#' gene_sets <- list(Pathway1 = gs1, Pathway2 = gs2)
#' 
#' # Use default fixed permutation method
#' result <- gsea(geneList=stats, gene_sets=gene_sets, nPerm=100)
#' 
#' # Use adaptive permutation for more accurate p-values
#' \dontrun{
#' result_adaptive <- gsea(geneList=stats, gene_sets=gene_sets, adaptive=TRUE)
#' }
#'
#' @export
gsea <- function(geneList, gene_sets, nPerm = 1000, exponent = 1.0, method = "sample",
                 adaptive = FALSE, minPerm = 1000, maxPerm = 100000, pvalThreshold = 0.1) {
    
    # Validate inputs
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("geneList must be a named numeric vector")
    }
    if (!is.list(gene_sets) || is.null(names(gene_sets))) {
        stop("gene_sets must be a named list")
    }
    
    method <- match.arg(method, c("sample", "permute"))
    
    # Ensure geneList is sorted
    if (is.unsorted(rev(geneList))) {
        warning("geneList is not sorted in descending order. Sorting it now.")
        geneList <- sort(geneList, decreasing = TRUE)
    }
    
    # Ensure gene_sets are character vectors
    gene_sets <- lapply(gene_sets, function(x) {
        if (!is.character(x)) {
            stop("Each element in gene_sets must be a character vector")
        }
        unique(x)
    })
    
    gene_set_names <- names(gene_sets)
    
    # Call appropriate C++ function
    if (adaptive) {
        result <- gsea_adaptive_cpp(geneList, gene_sets, gene_set_names, 
                                    minPerm, maxPerm, pvalThreshold, exponent, method)
    } else {
        result <- gsea_cpp(geneList, gene_sets, gene_set_names, nPerm, exponent, method)
    }
    
    # Rename columns to standard names
    names(result)[names(result) == "GeneSet"] <- "ID"
    names(result)[names(result) == "ES"] <- "enrichmentScore"
    names(result)[names(result) == "PValue"] <- "pvalue"
    names(result)[names(result) == "Size"] <- "setSize"

    # Sort by absolute NES (descending)
    if (nrow(result) > 0) {
        result <- result[order(abs(result$NES), decreasing = TRUE), ]
        rownames(result) <- NULL
    }
    
    return(result)
}


#' generic function for gene set enrichment analysis
#'
#'
#' @title gsea_gson
#' @param geneList order ranked geneList
#' @param gson GSON object
#' @param nPerm Number of permutations for p-value calculation (default: 1000). Only used when adaptive=FALSE.
#' @param exponent weight of each step
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param pvalueCutoff p value Cutoff
#' @param pAdjustMethod p value adjustment method
#' @param method Permutation method: "sample" (default) or "permute"
#' @param adaptive Logical. If TRUE, use adaptive early-stopping permutation. Default: FALSE.
#' @param minPerm Minimum permutations for adaptive mode (default: 1000).
#' @param maxPerm Maximum permutations for adaptive mode (default: 100000).
#' @param pvalThreshold P-value threshold for early stopping (default: 0.1).
#' @param verbose print message or not
#' @return gseaResult object
#' @author Guangchuang Yu
#' @export
gsea_gson <- function(geneList,
                 gson,
                 nPerm = 1000,
                 exponent = 1.0,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 method = "sample",
                 adaptive = FALSE,
                 minPerm = 1000,
                 maxPerm = 100000,
                 pvalThreshold = 0.1,
                 verbose = TRUE) {

    if (!inherits(gson, "GSON")) {
        stop("gson should be a GSON object")
    }

    ## query external ID to Term ID
    gene <- names(geneList)
    
    # Extract gene sets from GSON
    gsid2gene <- gson@gsid2gene
    
    # ID Match Check
    if (!check_gene_id(gene, gsid2gene)) {
        return(NULL)
    }

    # Prepare Gene Sets
    geneSets <- split(gsid2gene$gene, gsid2gene$gsid)
    
    # Filter by size
    idx <- get_geneSet_index(geneSets, minGSSize, maxGSSize)
    if (sum(idx) == 0) {
        if (verbose) {
            msg <- paste("No gene sets have size between", minGSSize, "and", maxGSSize, "...")
            message(msg)
            message("--> return NULL...")
        }
        return (NULL)
    }
    geneSets <- geneSets[idx]
    
    gsea_res <- gsea(geneList = geneList, 
                     gene_sets = geneSets, 
                     nPerm = nPerm, 
                     exponent = exponent, 
                     method = method,
                     adaptive = adaptive,
                     minPerm = minPerm,
                     maxPerm = maxPerm,
                     pvalThreshold = pvalThreshold)
                     
    if (is.null(gsea_res) || nrow(gsea_res) == 0) {
        return(NULL)
    }

    # Add Description
    gsid2name <- gson@gsid2name
    if (!is.null(gsid2name) && "ID" %in% names(gsea_res)) {
        description <- gsid2name$name[match(gsea_res$ID, gsid2name$gsid)]
        na_idx <- is.na(description)
        description[na_idx] <- gsea_res$ID[na_idx]
        gsea_res$Description <- description
    } else {
        if (!"Description" %in% names(gsea_res)) {
            gsea_res$Description <- gsea_res$ID
        }
    }

    # Calculate p.adjust
    gsea_res$p.adjust <- p.adjust(gsea_res$pvalue, method=pAdjustMethod)
    
    # Calculate qvalue
    gsea_res$qvalues <- calculate_qvalue(gsea_res$pvalue)
    
    # Filter by pvalueCutoff
    if (!is.null(pvalueCutoff)) {
        gsea_res <- gsea_res[gsea_res$pvalue <= pvalueCutoff, ]
    }
    
    if (nrow(gsea_res) == 0) {
        return(NULL)
    }
    
    # Reorder columns
    expected_cols <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues", "rank", "leading_edge", "core_enrichment")
    other_cols <- setdiff(names(gsea_res), expected_cols)
    gsea_res <- gsea_res[, c(expected_cols, other_cols)]
    
    # Set row names
    rownames(gsea_res) <- gsea_res$ID
    
    params <- list(pvalueCutoff = pvalueCutoff,
                   nPerm = nPerm,
                   pAdjustMethod = pAdjustMethod,
                   exponent = exponent,
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize)
                   
    res <- new("gseaResult",
               result = gsea_res,
               organism = if (!is.null(gson@species)) gson@species else "UNKNOWN",
               setType = if (!is.null(gson@gsname)) gsub(".*;", "", gson@gsname) else "UNKNOWN",
               geneSets = geneSets,
               geneList = geneList,
               keytype = if (!is.null(gson@keytype)) gson@keytype else "UNKNOWN",
               permScores = matrix(), 
               params = params,
               gene2Symbol = character(), 
               readable = FALSE
              )
              
    return(res)
}


