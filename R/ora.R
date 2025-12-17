#' Over-Representation Analysis (ORA)
#'
#' Perform over-representation analysis using hypergeometric test (Fisher's exact test).
#'
#' @param gene Character vector of differentially expressed genes (or gene list of interest).
#' @param universe Character vector of background genes (e.g., all genes in the platform).
#' @inheritParams enrichit_params
#'
#' @return A data.frame with columns:
#' \item{GeneSet}{Gene set name}
#' \item{SetSize}{Number of genes in the gene set (intersected with universe)}
#' \item{DEInSet}{Number of differentially expressed genes in the gene set}
#' \item{DESize}{Total number of differentially expressed genes in universe}
#' \item{PValue}{Raw p-value from hypergeometric test}
#'
#' @examples
#' # Example data
#' de_genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
#' all_genes <- paste0("Gene", 1:1000)
#' 
#' gs1 <- paste0("Gene", 1:50)
#' gs2 <- paste0("Gene", 51:150)
#' gs3 <- paste0("Gene", 151:300)
#' gene_sets <- list(Pathway1 = gs1, Pathway2 = gs2, Pathway3 = gs3)
#' 
#' result <- ora(gene=de_genes, gene_sets=gene_sets, universe=all_genes)
#' head(result)
#'
#' @export
ora <- function(gene, gene_sets, universe) {
    
    # Validate inputs
    if (!is.character(gene)) {
        stop("gene must be a character vector")
    }
    if (!is.character(universe)) {
        stop("universe must be a character vector")
    }
    
    gene_sets <- validate_gene_sets(gene_sets)
    
    # Remove duplicates
    gene <- unique(gene)
    universe <- unique(universe)
    
    gene_set_names <- names(gene_sets)
    
    # Call C++ function through Rcpp (using the Rcpp-generated ora_cpp function)
    result <- ora_cpp(gene, universe, gene_sets, gene_set_names)
    
    # Rename columns to standard names
    # C++ returns: GeneSet, SetSize, DEInSet, DESize, UniverseSize, PValue, geneID
    names(result)[names(result) == "GeneSet"] <- "ID"
    names(result)[names(result) == "PValue"] <- "pvalue"
    names(result)[names(result) == "DEInSet"] <- "Count"
    
    # Calculate derived columns
    # GeneRatio = Count / DESize
    result$GeneRatio <- paste0(result$Count, "/", result$DESize)
    
    # BgRatio = SetSize / UniverseSize
    result$BgRatio <- paste0(result$SetSize, "/", result$UniverseSize)
    
    # RichFactor = Count / SetSize
    result$RichFactor <- result$Count / result$SetSize
    
    # FoldEnrichment = (Count/DESize) / (SetSize/UniverseSize)
    result$FoldEnrichment <- (result$Count / result$DESize) / (result$SetSize / result$UniverseSize)
    
    # Sort by p-value
    result <- result[order(result$pvalue), ]
    rownames(result) <- NULL
    
    return(result)
}


