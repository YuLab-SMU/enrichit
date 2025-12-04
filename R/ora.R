#' Over-Representation Analysis
#'
#' Perform Over-Representation Analysis (ORA) for a set of genes against a universe and gene sets.
#'
#' @param gene A character vector of gene IDs (the query list).
#' @param universe A character vector of gene IDs (the background universe).
#' @param gene_sets A named list of gene sets. Each element is a character vector of genes.
#'
#' @return A data.frame with columns:
#' - **ID**: Gene set ID
#' - **GeneRatio**: Ratio of input genes that are in the gene set
#' - **BgRatio**: Ratio of background genes that are in the gene set
#' - **RichFactor**: Enrichment factor (Count/SetSize)
#' - **FoldEnrichment**: Fold enrichment (GeneRatio/BgRatio)
#' - **pvalue**: P-value from hypergeometric test
#' - **geneID**: Genes in the gene set that overlap with the input list
#' - **Count**: Number of overlapping genes
#'
#' @export
ora <- function(gene, universe, gene_sets) {
    
    # Validate inputs
    if (!is.character(gene)) stop("gene must be a character vector")
    if (!is.character(universe)) stop("universe must be a character vector")
    if (!is.list(gene_sets) || is.null(names(gene_sets))) stop("gene_sets must be a named list")
    
    # Ensure gene_sets are character vectors
    gene_sets <- lapply(gene_sets, function(x) {
        if (!is.character(x)) stop("Each element in gene_sets must be a character vector")
        unique(x)
    })
    
    gene_set_names <- names(gene_sets)
    
    # Call C++ function
    # Returns: GeneSet, SetSize, DEInSet, DESize, UniverseSize, PValue, geneID
    result <- ora_cpp(gene, universe, gene_sets, gene_set_names)
    
    # Rename columns to match clusterProfiler
    names(result)[names(result) == "GeneSet"] <- "ID"
    names(result)[names(result) == "PValue"] <- "pvalue"
    names(result)[names(result) == "DEInSet"] <- "Count"
    
    # GeneRatio: k/n (Count / DESize)
    result$GeneRatio <- paste0(result$Count, "/", result$DESize)
    
    # BgRatio: M/N (SetSize / UniverseSize)
    result$BgRatio <- paste0(result$SetSize, "/", result$UniverseSize)
    
    # RichFactor: Count / SetSize
    result$RichFactor <- result$Count / result$SetSize
    
    # FoldEnrichment: (Count/DESize) / (SetSize/UniverseSize)
    gene_ratio_num <- result$Count / result$DESize
    bg_ratio_num <- result$SetSize / result$UniverseSize
    result$FoldEnrichment <- gene_ratio_num / bg_ratio_num
    
    # Reorder columns
    cols <- c("ID", "GeneRatio", "BgRatio", "RichFactor", 
              "FoldEnrichment", "pvalue", "geneID", "Count")
    
    result <- result[, cols]
    
    # Sort by pvalue
    if (nrow(result) > 0) {
        result <- result[order(result$pvalue), ]
        rownames(result) <- NULL
    }
    
    return(result)
}
