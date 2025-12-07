
#' Calculate GSEA Running Enrichment Scores
#'
#' @param geneList a named numeric vector of gene statistics (e.g., t-statistics or log-fold changes), sorted in decreasing order.
#' @param geneSet a character vector of gene IDs belonging to the gene set.
#' @param exponent a numeric value defining the weight of the running enrichment score. Default is 1.
#' @param fortify logical. If TRUE, returns a data frame with columns `x`, `runningScore`, and `position`. 
#' If FALSE (default), returns the enrichment score (ES).
#'
#' @return If `fortify = TRUE`, a data frame containing the running enrichment scores and positions. 
#' If `fortify = FALSE`, a numeric value representing the Enrichment Score (ES).
#'
#' @export
#' @author Guangchuang Yu
gseaScores <- function(geneList, geneSet, exponent=1, fortify=FALSE) {
    geneSet <- intersect(geneSet, names(geneList))
    
    if (length(geneSet) == 0) {
        # If no genes in set match geneList, ES is 0? Or NA?
        # User code doesn't handle this explicitly, but intersection would be empty.
        # If empty, hits are all false.
        # But let's let the C++ code handle it or return 0 here.
        if (fortify) {
            return(data.frame(x = seq_along(geneList), runningScore = 0, position = 0))
        }
        return(0)
    }
    
    in_set <- names(geneList) %in% geneSet
    
    res <- gsea_scores_cpp(geneList, in_set, exponent)
    
    if (fortify) {
        return(res)
    }
    
    # Calculate ES
    # ES is the maximum deviation from zero
    max.ES <- max(res$runningScore)
    min.ES <- min(res$runningScore)
    
    if (abs(max.ES) > abs(min.ES)) {
        return(max.ES)
    } else {
        return(min.ES)
    }
}
