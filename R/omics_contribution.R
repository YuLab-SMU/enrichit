#' Get gene-level omics contribution for a specific pathway
#'
#' Extract the original multi-omics statistics for genes in a specific enriched pathway.
#'
#' @param res An \code{enrichResult} or \code{gseaResult} object.
#' @param agg An \code{omics_aggregated} object from \code{aggregate_omics()}.
#' @param pathway_id Character, the ID of the pathway to extract. If NULL, the top pathway is used.
#'
#' @return A data.frame containing the genes, their original omics statistics, the aggregated score, and whether they belong to the core enrichment.
#' @export
get_omics_contribution <- function(res, agg, pathway_id = NULL) {
    if (!inherits(agg, "omics_aggregated") || is.null(agg$original_matrix)) {
        stop("agg must be an omics_aggregated object containing original_matrix.")
    }
    
    res_df <- as.data.frame(res)
    if (nrow(res_df) == 0) {
        stop("The enrichment result is empty.")
    }
    
    if (is.null(pathway_id)) {
        pathway_id <- res_df$ID[1]
        message("pathway_id not provided. Using the top pathway: ", pathway_id)
    }
    
    if (!pathway_id %in% res_df$ID) {
        stop("pathway_id not found in the enrichment result.")
    }
    
    pathway_row <- res_df[res_df$ID == pathway_id, , drop = FALSE]
    
    # Extract genes in the pathway
    gene_col <- if ("core_enrichment" %in% colnames(pathway_row)) "core_enrichment" else "geneID"
    if (!gene_col %in% colnames(pathway_row)) {
        stop("Cannot find 'core_enrichment' or 'geneID' column in the result.")
    }
    
    pathway_genes <- unlist(strsplit(as.character(pathway_row[[gene_col]]), "/"))
    
    # Get the matrix and score
    mat <- agg$original_matrix
    all_genes <- rownames(mat)
    
    valid_genes <- intersect(pathway_genes, all_genes)
    if (length(valid_genes) == 0) {
        warning("No matching genes found between the pathway and the aggregated object.")
        return(NULL)
    }
    
    df <- as.data.frame(mat[valid_genes, , drop = FALSE])
    df$Aggregated_Score <- agg$score[valid_genes]
    
    if (!is.null(agg$pvalue)) {
        df$Aggregated_Pvalue <- agg$pvalue[valid_genes]
    }
    
    df$Feature <- valid_genes
    rownames(df) <- NULL
    
    # Reorder columns
    cols <- c("Feature", colnames(mat), "Aggregated_Score")
    if (!is.null(agg$pvalue)) cols <- c(cols, "Aggregated_Pvalue")
    
    return(df[, cols, drop = FALSE])
}

#' Classify pathway-level multi-omics patterns
#'
#' Compare merged enrichment results with single-omics enrichment results to classify the contribution pattern of each pathway.
#'
#' @param merged_res An \code{enrichResult} or \code{gseaResult} object from the merged multi-omics analysis.
#' @param single_res A named list of \code{enrichResult} or \code{gseaResult} objects from single-omics analyses.
#' @param p_cutoff Numeric, the significance cutoff. Default is 0.05.
#' @param by Character, the column to use for significance threshold. Default is "p.adjust".
#'
#' @return The \code{merged_res} object with an additional column \code{Omics_Pattern} in its result data.frame.
#' @export
classify_omics_pattern <- function(merged_res, single_res, p_cutoff = 0.05, by = "p.adjust") {
    if (!is.list(single_res) || is.null(names(single_res))) {
        stop("single_res must be a named list of enrichment results.")
    }
    
    merged_df <- as.data.frame(merged_res)
    if (!by %in% colnames(merged_df)) {
        stop(sprintf("Column '%s' not found in merged_res.", by))
    }
    
    # Prepare a lookup table for single-omics results
    single_lookups <- lapply(single_res, function(res) {
        df <- as.data.frame(res)
        if (nrow(df) > 0 && by %in% colnames(df)) {
            stats::setNames(df[[by]], df$ID)
        } else {
            stats::setNames(numeric(0), character(0))
        }
    })
    
    omics_names <- names(single_res)
    
    patterns <- sapply(seq_len(nrow(merged_df)), function(i) {
        path_id <- merged_df$ID[i]
        merged_p <- merged_df[[by]][i]
        
        if (is.na(merged_p) || merged_p > p_cutoff) {
            return("Not Significant")
        }
        
        sig_count <- 0
        sig_omics <- character(0)
        
        for (om in omics_names) {
            p_val <- single_lookups[[om]][path_id]
            if (!is.na(p_val) && p_val <= p_cutoff) {
                sig_count <- sig_count + 1
                sig_omics <- c(sig_omics, om)
            }
        }
        
        if (sig_count == 0) {
            return("Enhanced (1+1>2)")
        } else if (sig_count == 1) {
            return(paste0(sig_omics[1], "-Specific"))
        } else {
            return("Shared")
        }
    })
    
    if (inherits(merged_res, "enrichResult") || inherits(merged_res, "gseaResult")) {
        merged_res@result$Omics_Pattern <- patterns
    } else if (is.data.frame(merged_res)) {
        merged_res$Omics_Pattern <- patterns
    }
    return(merged_res)
}