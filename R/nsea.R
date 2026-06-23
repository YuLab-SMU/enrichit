#' Prepare network for repeated NSEA runs
#'
#' @param network edge list (data.frame with 2 or 3 columns) or sparse matrix.
#' @param directed logical, whether the network is directed. Default is FALSE.
#' @param normalize one of "column", "row", or "none". Default is "column".
#'
#' @return A sparse matrix (dgCMatrix) that has been properly formatted and normalized.
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#' @importFrom Matrix Diagonal
#' @importFrom stats setNames
#' @export
prepare_network <- function(network, directed = FALSE, normalize = "column") {
    if (inherits(network, "sparseMatrix")) {
        A <- network
    } else if (is.data.frame(network) || is.matrix(network)) {
        network <- as.data.frame(network)
        if (ncol(network) < 2) {
            stop("network must have at least 2 columns")
        }
        if (ncol(network) == 2) {
            network$weight <- 1
        }
        
        nodes <- unique(c(as.character(network[[1]]), as.character(network[[2]])))
        node_idx <- setNames(seq_along(nodes), nodes)
        
        i <- node_idx[as.character(network[[1]])]
        j <- node_idx[as.character(network[[2]])]
        x <- as.numeric(network[[3]])
        
        if (!directed) {
            i_all <- c(i, j)
            j_all <- c(j, i)
            x_all <- c(x, x)
        } else {
            i_all <- i
            j_all <- j
            x_all <- x
        }
        
        A <- Matrix::sparseMatrix(i = i_all, j = j_all, x = x_all,
                                  dims = c(length(nodes), length(nodes)),
                                  dimnames = list(nodes, nodes))
    } else {
        stop("network must be a data.frame, matrix or sparseMatrix")
    }
    
    if (normalize == "column") {
        cs <- Matrix::colSums(A)
        cs[cs == 0] <- 1 # avoid division by zero
        A <- A %*% Matrix::Diagonal(x = 1/cs)
    } else if (normalize == "row") {
        rs <- Matrix::rowSums(A)
        rs[rs == 0] <- 1
        A <- Matrix::Diagonal(x = 1/rs) %*% A
    }
    
    return(A)
}

#' Network-based Gene Set Enrichment Analysis
#'
#' @param geneList named numeric vector. Must be non-negative evidence scores.
#' @param network edge list (data.frame) or sparse matrix.
#' @param gene_sets list of gene sets.
#' @param p restart probability for RWR (default is 0.5).
#' @param minGSSize minimal size of each gene set for analyzing. default here is 10.
#' @param maxGSSize maximal size of genes annotated for testing. default here is 500.
#' @param threshold convergence threshold for RWR (default is 1e-9).
#' @param maxIter maximal number of RWR iterations (default is 100).
#' @param verbose logical, print messages.
#' @param ... other arguments passed to `gsea()`.
#'
#' @return A `data.frame` of NSEA results.
#' @export
nsea <- function(geneList,
                 network,
                 gene_sets,
                 p = 0.5,
                 minGSSize = 10,
                 maxGSSize = 500,
                 threshold = 1e-9,
                 maxIter = 100,
                 verbose = TRUE,
                 ...) {
    
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("geneList must be a named numeric vector")
    }
    
    if (any(geneList < 0)) {
        warning("geneList contains negative values. NSEA mode 'evidence' expects non-negative scores. Negative values will be propagated as is, but might violate RWR assumptions.")
    }
    
    if (verbose) message("Preparing network...")
    A <- prepare_network(network)
    
    nodes <- rownames(A)
    v <- rep(0, length(nodes))
    names(v) <- nodes
    
    common_nodes <- intersect(names(geneList), nodes)
    if (length(common_nodes) == 0) {
        stop("No overlapping genes between geneList and network.")
    }
    v[common_nodes] <- geneList[common_nodes]
    
    # Normalize restart vector
    sum_v <- sum(v)
    if (sum_v > 0) {
        v <- v / sum_v
    } else {
        stop("The sum of geneList scores in the network is zero.")
    }
    
    if (verbose) message("Running Random Walk with Restart (RWR)...")
    rwr_scores <- rwr_eigen_cpp(A, v, restart = p, threshold = threshold, max_iter = maxIter)
    names(rwr_scores) <- nodes
    
    rwr_scores <- sort(rwr_scores, decreasing = TRUE)
    
    if (verbose) message("Running GSEA...")
    res <- gsea(geneList = rwr_scores,
                gene_sets = gene_sets,
                minGSSize = minGSSize,
                maxGSSize = maxGSSize,
                scoreType = "pos",
                ...)
    
    return(res)
}

#' Network-based GSEA using a GSON object
#'
#' @param geneList named numeric vector. Must be non-negative evidence scores.
#' @param network edge list (data.frame) or sparse matrix.
#' @param gson a GSON object.
#' @param p restart probability for RWR (default is 0.5).
#' @param minGSSize minimal size of each gene set for analyzing. default here is 10.
#' @param maxGSSize maximal size of genes annotated for testing. default here is 500.
#' @param threshold convergence threshold for RWR (default is 1e-9).
#' @param maxIter maximal number of RWR iterations (default is 100).
#' @param verbose logical, print messages.
#' @param ... other arguments passed to `gsea_gson()`.
#'
#' @return A `gseaResult` object.
#' @export
nsea_gson <- function(geneList,
                      network,
                      gson,
                      p = 0.5,
                      minGSSize = 10,
                      maxGSSize = 500,
                      threshold = 1e-9,
                      maxIter = 100,
                      verbose = TRUE,
                      ...) {
    
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("geneList must be a named numeric vector")
    }
    
    if (verbose) message("Preparing network...")
    A <- prepare_network(network)
    
    nodes <- rownames(A)
    v <- rep(0, length(nodes))
    names(v) <- nodes
    
    common_nodes <- intersect(names(geneList), nodes)
    if (length(common_nodes) == 0) {
        stop("No overlapping genes between geneList and network.")
    }
    v[common_nodes] <- geneList[common_nodes]
    
    sum_v <- sum(v)
    if (sum_v > 0) {
        v <- v / sum_v
    } else {
        stop("The sum of geneList scores in the network is zero.")
    }
    
    if (verbose) message("Running Random Walk with Restart (RWR)...")
    rwr_scores <- rwr_eigen_cpp(A, v, restart = p, threshold = threshold, max_iter = maxIter)
    names(rwr_scores) <- nodes
    
    rwr_scores <- sort(rwr_scores, decreasing = TRUE)
    
    if (verbose) message("Running GSEA...")
    res <- gsea_gson(geneList = rwr_scores,
                     gson = gson,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     scoreType = "pos",
                     ...)
    
    return(res)
}
