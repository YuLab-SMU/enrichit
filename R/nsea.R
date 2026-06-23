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
#' @param geneList named numeric vector. In "evidence" mode, must be non-negative. In "signed" mode, can contain both positive and negative values.
#' @param network edge list (data.frame) or sparse matrix.
#' @param gene_sets list of gene sets.
#' @param mode character, either "evidence" (default) or "signed". If "signed", the network propagation runs separately for positive and negative values.
#' @param p restart probability for RWR (default is 0.5).
#' @param minGSSize minimal size of each gene set for analyzing. default here is 10.
#' @param maxGSSize maximal size of genes annotated for testing. default here is 500.
#' @param threshold convergence threshold for RWR (default is 1e-9).
#' @param maxIter maximal number of RWR iterations (default is 100).
#' @param verbose logical, print messages.
#' @param ... other arguments passed to `gsea()`.
#'
#' @return A `gseaResult` object of NSEA results.
#' @export
nsea <- function(geneList,
                 network,
                 gene_sets,
                 mode = c("evidence", "signed"),
                 p = 0.5,
                 minGSSize = 10,
                 maxGSSize = 500,
                 threshold = 1e-9,
                 maxIter = 100,
                 verbose = TRUE,
                 ...) {
    
    mode <- match.arg(mode)
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("geneList must be a named numeric vector")
    }
    
    if (mode == "evidence" && any(geneList < 0)) {
        warning("geneList contains negative values but mode is 'evidence'. Negative values will be propagated as is, which might violate RWR assumptions. Consider using mode = 'signed'.")
    }
    
    if (verbose) message("Preparing network...")
    A <- prepare_network(network)
    nodes <- rownames(A)
    
    common_nodes <- intersect(names(geneList), nodes)
    if (length(common_nodes) == 0) {
        stop("No overlapping genes between geneList and network.")
    }
    
    if (mode == "evidence") {
        v <- rep(0, length(nodes))
        names(v) <- nodes
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
        
    } else {
        # signed mode
        if (verbose) message("Running Signed RWR (Up and Down separately)...")
        
        v_up <- rep(0, length(nodes))
        names(v_up) <- nodes
        v_down <- rep(0, length(nodes))
        names(v_down) <- nodes
        
        genes_up <- common_nodes[geneList[common_nodes] > 0]
        genes_down <- common_nodes[geneList[common_nodes] < 0]
        
        v_up[genes_up] <- geneList[genes_up]
        v_down[genes_down] <- abs(geneList[genes_down])
        
        if (sum(v_up) > 0) v_up <- v_up / sum(v_up)
        if (sum(v_down) > 0) v_down <- v_down / sum(v_down)
        
        rwr_up <- rep(0, length(nodes))
        rwr_down <- rep(0, length(nodes))
        
        if (sum(v_up) > 0) rwr_up <- rwr_eigen_cpp(A, v_up, restart = p, threshold = threshold, max_iter = maxIter)
        if (sum(v_down) > 0) rwr_down <- rwr_eigen_cpp(A, v_down, restart = p, threshold = threshold, max_iter = maxIter)
        
        rwr_scores <- rwr_up - rwr_down
        names(rwr_scores) <- nodes
    }
    
    rwr_scores <- sort(rwr_scores, decreasing = TRUE)
    
    if (verbose) message("Running GSEA...")
    if (mode == "evidence") {
        res <- gsea(geneList = rwr_scores,
                    gene_sets = gene_sets,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    scoreType = "pos",
                    ...)
    } else {
        res <- gsea(geneList = rwr_scores,
                    gene_sets = gene_sets,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    scoreType = "std",
                    ...)
    }
    
    res_nsea <- new("nseaResult",
                    res,
                    network = network,
                    diffusion_scores = rwr_scores,
                    mode = mode)
    
    return(res_nsea)
}

#' Network-based GSEA using a GSON object
#'
#' @param geneList named numeric vector. In "evidence" mode, must be non-negative. In "signed" mode, can contain both positive and negative values.
#' @param network edge list (data.frame) or sparse matrix.
#' @param gson a GSON object.
#' @param mode character, either "evidence" (default) or "signed".
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
                      mode = c("evidence", "signed"),
                      p = 0.5,
                      minGSSize = 10,
                      maxGSSize = 500,
                      threshold = 1e-9,
                      maxIter = 100,
                      verbose = TRUE,
                      ...) {
    
    mode <- match.arg(mode)
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("geneList must be a named numeric vector")
    }
    
    if (verbose) message("Preparing network...")
    A <- prepare_network(network)
    nodes <- rownames(A)
    
    common_nodes <- intersect(names(geneList), nodes)
    if (length(common_nodes) == 0) {
        stop("No overlapping genes between geneList and network.")
    }
    
    if (mode == "evidence") {
        v <- rep(0, length(nodes))
        names(v) <- nodes
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
        
    } else {
        if (verbose) message("Running Signed RWR (Up and Down separately)...")
        
        v_up <- rep(0, length(nodes))
        names(v_up) <- nodes
        v_down <- rep(0, length(nodes))
        names(v_down) <- nodes
        
        genes_up <- common_nodes[geneList[common_nodes] > 0]
        genes_down <- common_nodes[geneList[common_nodes] < 0]
        
        v_up[genes_up] <- geneList[genes_up]
        v_down[genes_down] <- abs(geneList[genes_down])
        
        if (sum(v_up) > 0) v_up <- v_up / sum(v_up)
        if (sum(v_down) > 0) v_down <- v_down / sum(v_down)
        
        rwr_up <- rep(0, length(nodes))
        rwr_down <- rep(0, length(nodes))
        
        if (sum(v_up) > 0) rwr_up <- rwr_eigen_cpp(A, v_up, restart = p, threshold = threshold, max_iter = maxIter)
        if (sum(v_down) > 0) rwr_down <- rwr_eigen_cpp(A, v_down, restart = p, threshold = threshold, max_iter = maxIter)
        
        rwr_scores <- rwr_up - rwr_down
        names(rwr_scores) <- nodes
    }
    
    rwr_scores <- sort(rwr_scores, decreasing = TRUE)
    
    if (verbose) message("Running GSEA...")
    if (mode == "evidence") {
        res <- gsea_gson(geneList = rwr_scores,
                         gson = gson,
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize,
                         scoreType = "pos",
                         ...)
    } else {
        res <- gsea_gson(geneList = rwr_scores,
                         gson = gson,
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize,
                         scoreType = "std",
                         ...)
    }
    
    res_nsea <- new("nseaResult",
                    res,
                    network = network,
                    diffusion_scores = rwr_scores,
                    mode = mode)
    
    return(res_nsea)
}
