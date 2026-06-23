#' Over-Representation Analysis (ORA)
#'
#' Perform over-representation analysis using hypergeometric test (Fisher's exact test).
#'
#' @param gene Character vector of differentially expressed genes (or gene list of interest).
#' @param universe Character vector of background genes (e.g., all genes in the platform).
#' @param weight A named numeric vector of weights for background genes. If provided, Weighted ORA will be performed using Wallenius' noncentral hypergeometric distribution (requires 'BiasedUrn' package). The names should match the universe genes.
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
ora <- function(gene, gene_sets, universe, weight = NULL) {
    
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
    
    if (!is.null(weight)) {
        if (!is.numeric(weight) || is.null(names(weight))) {
            stop("weight must be a named numeric vector")
        }
        rlang::check_installed('BiasedUrn', 'for weighted ORA.')
        
        weight_map <- weight[universe]
        names(weight_map) <- universe
        if (any(is.na(weight_map))) {
            warning("Some background genes do not have weights. Assigning median weight.")
            weight_map[is.na(weight_map)] <- median(weight, na.rm = TRUE)
        }
        if (any(weight_map <= 0, na.rm = TRUE)) {
            warning("Weights must be strictly positive for BiasedUrn. Adjusting non-positive weights.")
            min_pos <- min(weight_map[weight_map > 0], na.rm = TRUE)
            if (is.na(min_pos) || min_pos <= 0) min_pos <- 1e-6
            weight_map[weight_map <= 0] <- min_pos
        }
        
        pvals <- numeric(nrow(result))
        for (i in seq_len(nrow(result))) {
            gs_name <- result$GeneSet[i]
            gset <- gene_sets[[gs_name]]
            genes_in_set <- intersect(gset, universe)
            
            m1 <- length(genes_in_set)
            m2 <- length(universe) - m1
            x <- result$DEInSet[i]
            n <- result$DESize[i]
            
            if (m1 == 0 || x == 0) {
                pvals[i] <- 1.0
                next
            }
            
            w_in <- mean(weight_map[genes_in_set], na.rm = TRUE)
            w_out <- mean(weight_map[setdiff(universe, genes_in_set)], na.rm = TRUE)
            
            odds <- w_in / w_out
            if (is.na(odds) || !is.finite(odds)) odds <- 1.0
            
            pvals[i] <- 1 - BiasedUrn::pWNCHypergeo(x - 1, m1, m2, n, odds)
        }
        result$PValue <- pvals
    }
    
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


