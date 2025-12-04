#' Gene Set Enrichment Analysis (GSEA)
#'
#' Perform Gene Set Enrichment Analysis (GSEA) using a ranked gene list.
#'
#' @param genelist A named numeric vector of gene statistics (e.g., log fold change), ranked in descending order.
#' @param gene_sets A named list of gene sets. Each element is a character vector of genes.
#' @param nPerm Number of permutations for p-value calculation (default: 1000).
#' @param exponent Weighting exponent for enrichment score (default: 1.0).
#' @param method Permutation method: "sample" (default) for random gene set sampling (faster, similar to fgsea), 
#' or "permute" for label permutation (slower, standard GSEA).
#'
#' @return A data.frame with columns:
#' - **GeneSet**: Gene set name
#' - **ES**: Enrichment Score
#' - **NES**: Normalized Enrichment Score
#' - **PValue**: Empirical p-value from permutation test
#' - **Size**: Size of the gene set (number of genes found in genelist)
#' - **rank**: Rank at which the maximum enrichment score is attained
#' - **leading_edge**: Leading edge statistics (tags, list, signal)
#' - **core_enrichment**: Genes in the leading edge, separated by '/'
#'
#' @examples
#' ```{r}
#' # Example data
#' stats <- rnorm(1000)
#' names(stats) <- paste0("Gene", 1:1000)
#' stats <- sort(stats, decreasing = TRUE)
#' 
#' gs1 <- paste0("Gene", 1:50)
#' gs2 <- paste0("Gene", 500:550)
#' gene_sets <- list(Pathway1 = gs1, Pathway2 = gs2)
#' 
#' # Use default sampling method
#' result <- gsea(genelist=stats, gene_sets=gene_sets, nPerm=100)
#' 
#' # Use label permutation method
#' result_perm <- gsea(genelist=stats, gene_sets=gene_sets, nPerm=100, method="permute")
#' ```
#'
#' @export
gsea <- function(genelist, gene_sets, nPerm = 1000, exponent = 1.0, method = "sample") {
    
    # Validate inputs
    if (!is.numeric(genelist) || is.null(names(genelist))) {
        stop("genelist must be a named numeric vector")
    }
    if (!is.list(gene_sets) || is.null(names(gene_sets))) {
        stop("gene_sets must be a named list")
    }
    
    method <- match.arg(method, c("sample", "permute"))
    
    # Ensure genelist is sorted
    if (is.unsorted(rev(genelist))) {
        warning("genelist is not sorted in descending order. Sorting it now.")
        genelist <- sort(genelist, decreasing = TRUE)
    }
    
    # Ensure gene_sets are character vectors
    gene_sets <- lapply(gene_sets, function(x) {
        if (!is.character(x)) {
            stop("Each element in gene_sets must be a character vector")
        }
        unique(x)
    })
    
    gene_set_names <- names(gene_sets)
    
    # Call C++ function
    result <- gsea_cpp(genelist, gene_sets, gene_set_names, nPerm, exponent, method)
    
    # Sort by absolute NES (descending)
    if (nrow(result) > 0) {
        result <- result[order(abs(result$NES), decreasing = TRUE), ]
        rownames(result) <- NULL
    }
    
    return(result)
}
