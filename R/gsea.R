#' Gene Set Enrichment Analysis (GSEA)
#'
#' Perform Gene Set Enrichment Analysis (GSEA) using a ranked gene list.
#'
#' @inheritParams enrichit_params
#' @param eps Epsilon for multilevel methods (default: 1e-10). Sets the smallest p-value that can be estimated.
#' @param sampleSize Sample size for multilevel methods (default: 101).
#' @param seed Random seed for reproducibility (default: FALSE). If FALSE, a random seed is generated.
#' @param nPermSimple Number of permutations for the simple method (default: 1000).
#' @param scoreType Type of enrichment score calculation: "std", "pos", "neg" (default: "std").
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
gsea <- function(geneList, gene_sets, 
                 minGSSize = 10,
                 maxGSSize = 500,
                 nPerm = 1000, 
                 exponent = 1.0, 
                 method = "multilevel",
                 adaptive = FALSE, 
                 minPerm = 101, 
                 maxPerm = 100000, 
                 pvalThreshold = 0.1, 
                 eps = 1e-10,
                 sampleSize = 101,
                 seed = FALSE,
                 nPermSimple = 1000,
                 scoreType = "std",
                 verbose = TRUE) {
    
    # Validate inputs
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("geneList must be a named numeric vector")
    }
    
    gene_sets <- validate_gene_sets(gene_sets)
    
    method <- match.arg(method, c("sample", "permute", "multilevel"))
    scoreType <- match.arg(scoreType, c("std", "pos", "neg"))
    
    # Ensure geneList is sorted
    if (is.unsorted(rev(geneList))) {
        warning("geneList is not sorted in descending order. Sorting it now.")
        geneList <- sort(geneList, decreasing = TRUE)
    }
    
    # Filter by size
    idx <- get_geneSet_index(gene_sets, minGSSize, maxGSSize)
    if (sum(idx) == 0) {
        if (verbose) {
            msg <- paste("No gene sets have size between", minGSSize, "and", maxGSSize, "...")
            message(msg)
            message("--> return NULL...")
        }
        return (NULL)
    }
    gene_sets <- gene_sets[idx]   

    gene_set_names <- names(gene_sets)
    
    # Call appropriate C++ function
    if (method == "multilevel") {
        if (isFALSE(seed)) {
            seed <- sample.int(1e9, 1)
        }
        result <- gsea_multilevel_cpp(geneList = geneList, gene_sets = gene_sets, gene_set_names = gene_set_names,
                                 minPerm = minPerm, maxPerm = maxPerm, pvalThreshold = pvalThreshold,
                                 exponent = exponent, method = method, eps = eps, sampleSize = sampleSize, seed = seed,
                                 nPermSimple = nPermSimple, scoreType = scoreType)
    } else if (adaptive) {
        result <- gsea_adaptive_cpp(geneList, gene_sets, gene_set_names, 
                                    minPerm, maxPerm, pvalThreshold, exponent, method)
    } else {
        result <- gsea_cpp(geneList, gene_sets, gene_set_names, nPerm, exponent, method)
    }
    
    # Rename columns to standard names
    if (!"ID" %in% names(result) && "GeneSet" %in% names(result)) names(result)[names(result) == "GeneSet"] <- "ID"
    if (!"enrichmentScore" %in% names(result) && "ES" %in% names(result)) names(result)[names(result) == "ES"] <- "enrichmentScore"
    if (!"pvalue" %in% names(result) && "PValue" %in% names(result)) names(result)[names(result) == "PValue"] <- "pvalue"
    if (!"setSize" %in% names(result) && "Size" %in% names(result)) names(result)[names(result) == "Size"] <- "setSize"

    # Add setSize if missing
    if (!"setSize" %in% names(result)) {
        # Calculate size based on gene_sets (filtered)
        # Note: gene_sets order must match result ID order
        # Assuming result ID matches gene_set_names order passed to C++
        # If C++ preserves order (it does), we can map
        set_sizes <- sapply(gene_sets, length)
        result$setSize <- set_sizes[result$ID]
    }

    # Add NES if missing (placeholder)
    if (!"NES" %in% names(result)) {
        result$NES <- NA
    }

    # Sort by absolute NES (descending) or pvalue (ascending)
    if (nrow(result) > 0) {
        if (!all(is.na(result$NES))) {
            result <- result[order(abs(result$NES), decreasing = TRUE), ]
        } else {
             result <- result[order(result$pvalue, decreasing = FALSE), ]
        }
    }
    
    return(result)
}


#' generic function for gene set enrichment analysis
#'
#'
#' @title gsea_gson
#' @inheritParams enrichit_params
#' @param ... Additional parameters passed to gsea()
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
                 method = "multilevel",
                 adaptive = FALSE,
                 minPerm = 101,
                 maxPerm = 100000,
                 pvalThreshold = 0.1,
                 verbose = TRUE,
                 ...) {

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
    
    gsea_res <- gsea(geneList = geneList, 
                     gene_sets = geneSets, 
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     nPerm = nPerm, 
                     exponent = exponent, 
                     method = method,
                     adaptive = adaptive,
                     minPerm = minPerm,
                     maxPerm = maxPerm,
                     pvalThreshold = pvalThreshold,
                     verbose = verbose,
                     ...)
                     
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


