#' Get cached contribution tables from a `mnseaResult`
#'
#' @param res A `mnseaResult` object.
#' @param pathway_id Optional pathway ID. If `NULL`, returns all pathways for
#'   `level = "pathway"` and uses the top pathway for `level = "feature"`.
#' @param level One of `"pathway"` or `"feature"`.
#'
#' @return A data.frame containing cached contribution information.
#' @export
get_mnsea_contribution <- function(res, pathway_id = NULL, level = c("pathway", "feature")) {
    level <- match.arg(level)
    
    if (!inherits(res, "mnseaResult")) {
        stop("res must be a mnseaResult object.")
    }
    
    if (level == "pathway") {
        df <- res@pathway_contribution
        if (is.null(pathway_id)) {
            return(df)
        }
        return(df[df$ID == pathway_id, , drop = FALSE])
    }
    
    if (is.null(pathway_id)) {
        if (nrow(res@result) == 0) {
            stop("The mnseaResult is empty.")
        }
        pathway_id <- as.character(res@result$ID[1])
        message("pathway_id not provided. Using the top pathway: ", pathway_id)
    }
    
    df <- res@feature_contribution
    df[df$ID == pathway_id, , drop = FALSE]
}
