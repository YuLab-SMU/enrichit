#' geneID generic
#'
#' @param x enrichResult object
#' @return 'geneID' return the 'geneID' column of the enriched result which can be converted to data.frame via 'as.data.frame'
#' @export
#' @examples
#' \donttest{
#' data(geneList, package="DOSE")
#' de <- names(geneList)[1:100]
#' x <- DOSE::enrichDO(de)
#' geneID(x)
#' }
geneID <- function(x) {
   UseMethod("geneID", x)
}

#' geneInCategory generic
#'
#' @param x enrichResult
#' @return 'geneInCategory' return a list of genes, by spliting the input gene vector to enriched functional categories
#' @export
#' @examples
#' \donttest{
#' data(geneList, package="DOSE")
#' de <- names(geneList)[1:100]
#' x <- DOSE::enrichDO(de)
#' geneInCategory(x)
#' }
geneInCategory <- function(x) {
   UseMethod("geneInCategory", x)
}

