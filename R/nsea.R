#' Prepare network for repeated NSEA runs
#'
#' @param network edge list (data.frame with 2 or 3 columns), a sparse
#'   matrix, or a `mechgraph` object. For `mechgraph`, its `edges` table must
#'   have `from`/`to` columns; a numeric `score` or `weight` column, if
#'   present, is used as the edge weight, otherwise unit weights are used.
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
    } else {
        # mechgraph objects carry the network as an edge table with
        # from/to/type plus optional evidence columns; convert them to the
        # plain edge-list form consumed below (duck-typed on the class name,
        # so enrichit does not hard-depend on the mechgraph package).
        if (inherits(network, "mechgraph")) {
            edges <- network$edges
            if (is.null(edges) || !is.data.frame(edges) ||
                !all(c("from", "to") %in% names(edges))) {
                stop("mechgraph network must have edges with 'from' and 'to' columns")
            }
            wt_col <- if ("score" %in% names(edges) && is.numeric(edges$score)) {
                "score"
            } else if ("weight" %in% names(edges) && is.numeric(edges$weight)) {
                "weight"
            } else {
                NULL
            }
            network <- if (is.null(wt_col)) {
                edges[, c("from", "to"), drop = FALSE]
            } else {
                edges[, c("from", "to", wt_col), drop = FALSE]
            }
        }
        if (!is.data.frame(network) && !is.matrix(network)) {
            stop("network must be a data.frame, matrix, sparseMatrix or mechgraph")
        }
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

build_nsea_result <- function(base_result,
                              gene_sets,
                              rwr_scores,
                              network,
                              mode,
                              iter,
                              p,
                              organism = "UNKNOWN",
                              setType = "UNKNOWN",
                              keytype = "UNKNOWN",
                              params = NULL) {
    if (inherits(base_result, "gseaResult")) {
        result_df <- base_result@result
        gene_sets <- base_result@geneSets
        organism <- base_result@organism
        setType <- base_result@setType
        keytype <- base_result@keytype
        params <- base_result@params
        gene_list <- base_result@geneList
        perm_scores <- base_result@permScores
        gene2symbol <- base_result@gene2Symbol
        readable <- base_result@readable
        termsim <- base_result@termsim
        method <- base_result@method
        dr <- base_result@dr
    } else if (is.data.frame(base_result)) {
        result_df <- base_result
        if (is.null(params)) {
            params <- list(pvalueCutoff = 1.0, pAdjustMethod = "BH")
        }
        gene_list <- rwr_scores
        perm_scores <- matrix(0, nrow = 0, ncol = 0)
        gene2symbol <- character(0)
        readable <- FALSE
        termsim <- matrix(0, nrow = 0, ncol = 0)
        method <- "NSEA"
        dr <- list()
    } else {
        stop("base_result must be a gseaResult or data.frame.")
    }
    
    if (!"Description" %in% colnames(result_df)) {
        result_df$Description <- result_df$ID
    }
    if (!"p.adjust" %in% colnames(result_df)) {
        result_df$p.adjust <- stats::p.adjust(result_df$pvalue, method = params$pAdjustMethod %||% "BH")
    }
    if (!"qvalue" %in% colnames(result_df)) {
        result_df$qvalue <- calculate_qvalue(result_df$pvalue)
    }
    rownames(result_df) <- as.character(result_df$ID)
    
    new("nseaResult",
        result = result_df,
        organism = organism,
        setType = setType,
        geneSets = gene_sets,
        geneList = gene_list,
        keytype = keytype,
        permScores = perm_scores,
        params = params,
        gene2Symbol = gene2symbol,
        readable = readable,
        termsim = termsim,
        method = method,
        dr = dr,
        network = network,
        diffusion_scores = rwr_scores,
        mode = mode,
        iterations = as.integer(iter),
        restart_prob = p)
}

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0) y else x
}

## ---------------------------------------------------------------------------
## Whole-pipeline permutation significance for NSEA
##
## The classical NSEA pipeline ran GSEA's label-permutation test on the
## network-diffused scores. Because diffusion induces strong autocorrelation
## between neighboring genes' scores, that test violates the exchangeability
## assumption of GSEA's permutation null and yields anti-conservative p-values
## (empirically ~2.5x the nominal false-positive rate).
##
## The corrected approach computes the null distribution of the enrichment
## score by re-running the *entire* pipeline under the null: permute the gene
## labels of the input scores, re-diffuse over the network with RWR, and
## recompute the enrichment score. P-values and NES are then derived from this
## whole-pipeline null, which automatically accounts for the smoothing induced
## by the network diffusion.
## ---------------------------------------------------------------------------

# Robust node-name extraction (guards against Matrix accessor quirks on
# serialized sparse matrices).
.nsea_network_nodes <- function(A) {
    nodes <- rownames(A)
    if (is.null(nodes) || length(nodes) == 0) {
        dn <- tryCatch(A@Dimnames, error = function(e) NULL)
        if (is.list(dn) && length(dn) >= 1 && length(dn[[1]]) == nrow(A)) {
            nodes <- as.character(dn[[1]])
        }
    }
    if (is.null(nodes) || length(nodes) == 0) {
        stop("Network matrix has no node names; provide an edge-list data.frame ",
             "or a sparse matrix with dimnames.")
    }
    nodes
}

# One RWR diffusion of the seed values over the network.
# Returns the diffused scores (named, NOT sorted) or NULL when the seed vector
# is degenerate (evidence mode with non-positive total evidence).
# In "signed" mode, positive and negative seeds are propagated in a single
# solve on w = v_up/sum(v_up) - v_down/sum(|v_down|); because RWR is linear in
# the restart vector, RWR(w) == RWR(v_up) - RWR(v_down) exactly.
.nsea_diffuse <- function(A, nodes, geneList, mode, p, threshold, maxIter) {
    n <- length(nodes)
    common <- intersect(names(geneList), nodes)
    if (length(common) == 0) {
        return(NULL)
    }
    
    if (mode == "evidence") {
        v <- rep(0, n)
        names(v) <- nodes
        v[common] <- geneList[common]
        s <- sum(v)
        if (!is.finite(s) || s <= 0) {
            return(NULL)
        }
        v <- v / s
        r <- rwr_eigen_cpp(A, v, restart = p, threshold = threshold, max_iter = maxIter)
        res <- r$score
        names(res) <- nodes
        attr(res, "iterations") <- as.integer(r$iterations)
        return(res)
    }
    
    # signed mode: single linear solve (see note above)
    vals <- geneList[common]
    up <- vals[vals > 0]
    down <- vals[vals < 0]
    if (length(up) == 0 && length(down) == 0) {
        return(NULL)
    }
    w <- rep(0, n)
    names(w) <- nodes
    up_sum <- sum(up)
    dn_sum <- sum(abs(down))
    if (length(up) > 0) {
        w[common][vals > 0] <- vals[vals > 0] / up_sum
    }
    if (length(down) > 0) {
        w[common][vals < 0] <- vals[vals < 0] / dn_sum   # negative weights
    }
    r <- rwr_eigen_cpp(A, w, restart = p, threshold = threshold, max_iter = maxIter)
    res <- r$score
    names(res) <- nodes
    attr(res, "iterations") <- as.integer(r$iterations)
    return(res)
}

# TF-IDF style gene specificity weights over the network nodes (deterministic).
.nsea_specific_weights <- function(gene_sets, nodes) {
    N <- length(gene_sets)
    gene_freq <- table(unlist(gene_sets))
    df_vals <- as.numeric(gene_freq[nodes])
    df_vals[is.na(df_vals)] <- 1 # minimum frequency 1
    w <- log(N / df_vals)
    w[w < 0] <- 0
    names(w) <- nodes
    w
}

# Core of the corrected NSEA significance test.
# Returns a list with the observed ranked diffused scores, per-pathway
# enrichment scores, p-values and NES computed from the whole-pipeline
# permutation null, plus iteration counts.
.nsea_permutation_test <- function(geneList, A, nodes, gene_sets,
                                   mode, p, threshold, maxIter,
                                   specific_weight,
                                   minGSSize, maxGSSize,
                                   nPerm, exponent, scoreType,
                                   seed = NULL, verbose = TRUE) {
    
    common_nodes <- intersect(names(geneList), nodes)
    if (length(common_nodes) == 0) {
        stop("No overlapping genes between geneList and network.")
    }
    
    if (mode == "evidence" && any(geneList < 0)) {
        warning("geneList contains negative values but mode is 'evidence'. ",
                "Negative values will be propagated as is, which might violate RWR ",
                "assumptions. Consider using mode = 'signed'.")
    }
    
    ## observed diffusion
    if (verbose) message("Running Random Walk with Restart (RWR)...")
    obs <- .nsea_diffuse(A, nodes, geneList, mode, p, threshold, maxIter)
    if (is.null(obs)) {
        stop("The total evidence of geneList over the network is not positive. ",
             "In 'evidence' mode the geneList must be non-negative (or have a ",
             "positive total); for signed inputs use mode = 'signed'.")
    }
    iter_obs <- as.integer(attr(obs, "iterations") %||% NA_integer_)
    
    spw <- NULL
    if (specific_weight) {
        if (verbose) message("Applying gene specificity weighting (TF-IDF)...")
        spw <- .nsea_specific_weights(gene_sets, nodes)
        obs <- obs * spw
    }
    obs <- sort(obs, decreasing = TRUE)
    
    ## intersect gene sets with the ranked list and apply size filters
    gs <- lapply(gene_sets, intersect, names(obs))
    sizes <- vapply(gs, length, integer(1))
    keep <- sizes >= minGSSize & sizes <= maxGSSize
    gs <- gs[keep]
    sizes <- sizes[keep]
    if (length(gs) == 0) {
        if (verbose) {
            message("No gene sets have size between ", minGSSize, " and ",
                    maxGSSize, "...")
            message("--> return NULL...")
        }
        return(NULL)
    }
    ids <- names(gs)
    
    es_obs <- gsea_es_all_cpp(obs, unname(gs), exponent, scoreType)
    
    ## whole-pipeline permutations: permute gene labels -> re-diffuse -> ES
    if (!is.null(seed)) {
        set.seed(seed)
    }
    n_set <- length(gs)
    es_perm <- matrix(0, nrow = n_set, ncol = nPerm)
    genes_pool <- names(geneList)
    for (b in seq_len(nPerm)) {
        perm_gl <- setNames(geneList, sample(genes_pool))
        r <- .nsea_diffuse(A, nodes, perm_gl, mode, p, threshold, maxIter)
        if (is.null(r)) {
            next   # degenerate permutation (should not happen for valid input)
        }
        if (specific_weight) {
            r <- r * spw
        }
        r <- sort(r, decreasing = TRUE)
        es_perm[, b] <- gsea_es_all_cpp(r, unname(gs), exponent, scoreType)
        if (verbose && (b %% 100 == 0 || b == nPerm)) {
            message("  permutation ", b, "/", nPerm)
        }
    }
    
    ## p-values from the whole-pipeline null
    if (scoreType == "pos") {
        pv <- (rowSums(es_perm >= es_obs) + 1) / (nPerm + 1)
    } else if (scoreType == "neg") {
        pv <- (rowSums(es_perm <= es_obs) + 1) / (nPerm + 1)
    } else {
        pv <- (rowSums(abs(es_perm) >= abs(es_obs)) + 1) / (nPerm + 1)
    }
    pv <- pmin(pv, 1)
    
    ## NES: normalize ES by the mean absolute null ES (set-size aware)
    mean_abs <- rowMeans(abs(es_perm))
    nes <- ifelse(mean_abs > 0, es_obs / mean_abs, NA_real_)
    
    list(ids = ids, sizes = sizes, gs = gs, es = es_obs, nes = nes,
         pvalue = pv, ranked = obs, iterations = iter_obs,
         nPerm = as.integer(nPerm), seed = seed)
}

# Build the result data.frame for the whole-pipeline NSEA.
.nsea_build_result_df <- function(perm_test, exponent, scoreType,
                                  gsid2name = NULL, pAdjustMethod = "BH") {
    ranked <- perm_test$ranked
    gs <- perm_test$gs
    
    ledge <- lapply(gs, function(g) {
        gsea_leading_edge_details(ranked, g, exponent = exponent, scoreType = scoreType)
    })
    
    df <- data.frame(
        ID = perm_test$ids,
        setSize = perm_test$sizes,
        enrichmentScore = perm_test$es,
        NES = perm_test$nes,
        pvalue = perm_test$pvalue,
        rank = vapply(ledge, `[[`, integer(1), "rank"),
        leading_edge = vapply(ledge, `[[`, character(1), "leading_edge"),
        core_enrichment = vapply(ledge, `[[`, character(1), "core_enrichment"),
        stringsAsFactors = FALSE
    )
    
    if (!is.null(gsid2name)) {
        nm <- gsid2name$name[match(df$ID, gsid2name$gsid)]
        nm[is.na(nm)] <- df$ID[is.na(nm)]
        df$Description <- nm
    } else {
        df$Description <- df$ID
    }
    
    df$p.adjust <- stats::p.adjust(df$pvalue, method = pAdjustMethod)
    df$qvalue <- calculate_qvalue(df$pvalue)
    
    expected <- c("ID", "Description", "setSize", "enrichmentScore", "NES",
                  "pvalue", "p.adjust", "qvalue", "rank", "leading_edge",
                  "core_enrichment")
    df <- df[, c(expected, setdiff(names(df), expected))]
    df <- df[order(abs(df$NES), decreasing = TRUE), ]
    rownames(df) <- as.character(df$ID)
    df
}

#' Network-based Gene Set Enrichment Analysis
#'
#' @param geneList named numeric vector. In "evidence" mode, must be non-negative. In "signed" mode, can contain both positive and negative values.
#' @param network edge list (data.frame) or sparse matrix.
#' @param gene_sets list of gene sets.
#' @param mode character, either "evidence" (default) or "signed". If "signed", the network propagation runs separately for positive and negative values.
#' @param p restart probability for RWR (default is 0.5).
#' @param specific_weight logical, whether to apply gene specificity weighting (TF-IDF style) based on gene frequencies in `gene_sets`. Default is FALSE.
#' @param minGSSize minimal size of each gene set for analyzing. default here is 10.
#' @param maxGSSize maximal size of genes annotated for testing. default here is 500.
#' @param threshold convergence threshold for RWR (default is 1e-9).
#' @param maxIter maximal number of RWR iterations (default is 100).
#' @param significance character, one of "whole_pipeline" (default) or "internal".
#'   "whole_pipeline" computes p-values and NES from a permutation null that
#'   re-runs the full pipeline (permute gene labels -> re-diffuse -> recompute
#'   enrichment scores) and is statistically calibrated. "internal" runs GSEA's
#'   internal label-permutation test on the diffused scores; this is the legacy
#'   behaviour and is **anti-conservative** (it ignores the correlation induced
#'   by network diffusion) — use it only for comparison/debugging.
#' @param nPerm number of whole-pipeline permutations used to estimate p-values
#'   and NES (default: 1000). The smallest estimable p-value is 1/(nPerm + 1).
#' @param seed random seed for the permutation null (default: NULL, use the
#'   current R RNG state). Set a numeric seed for reproducible results.
#' @param exponent weight of each step in the enrichment score (default: 1).
#' @param verbose logical, print messages.
#' @param ... for `significance = "internal"`, other arguments passed to
#'   `gsea()` (e.g. `nPerm`, `method`, `nPermSimple`). For
#'   `significance = "whole_pipeline"` additional arguments are ignored.
#'
#' @return A `nseaResult` object of NSEA results.
#' @export
nsea <- function(geneList,
                 network,
                 gene_sets,
                 mode = c("evidence", "signed"),
                 p = 0.5,
                 specific_weight = FALSE,
                 minGSSize = 10,
                 maxGSSize = 500,
                 threshold = 1e-9,
                 maxIter = 100,
                 significance = c("whole_pipeline", "internal"),
                 nPerm = 1000,
                 seed = NULL,
                 exponent = 1,
                 verbose = TRUE,
                 ...) {
    
    mode <- match.arg(mode)
    significance <- match.arg(significance)
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("geneList must be a named numeric vector")
    }
    
    if (verbose) message("Preparing network...")
    A <- prepare_network(network)
    nodes <- .nsea_network_nodes(A)
    
    ## ---- legacy internal-permutation path (kept for comparison/debugging) ----
    if (significance == "internal") {
        if (mode == "evidence" && any(geneList < 0)) {
            warning("geneList contains negative values but mode is 'evidence'. Negative values will be propagated as is, which might violate RWR assumptions. Consider using mode = 'signed'.")
        }
        
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
            rwr_res <- rwr_eigen_cpp(A, v, restart = p, threshold = threshold, max_iter = maxIter)
            rwr_scores <- rwr_res$score
            names(rwr_scores) <- nodes
            iter <- as.integer(rwr_res$iterations)
            
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
            
            iter_up <- 0
            iter_down <- 0
            
            if (sum(v_up) > 0) {
                r_up <- rwr_eigen_cpp(A, v_up, restart = p, threshold = threshold, max_iter = maxIter)
                rwr_up <- r_up$score
                iter_up <- r_up$iterations
            }
            if (sum(v_down) > 0) {
                r_down <- rwr_eigen_cpp(A, v_down, restart = p, threshold = threshold, max_iter = maxIter)
                rwr_down <- r_down$score
                iter_down <- r_down$iterations
            }
            
            rwr_scores <- rwr_up - rwr_down
            names(rwr_scores) <- nodes
            iter <- as.integer(max(iter_up, iter_down))
        }
        
        if (specific_weight) {
            if (verbose) message("Applying gene specificity weighting (TF-IDF)...")
            N <- length(gene_sets)
            gene_freq <- table(unlist(gene_sets))
            df_vals <- as.numeric(gene_freq[names(rwr_scores)])
            df_vals[is.na(df_vals)] <- 1 # Minimum frequency 1
            w <- log(N / df_vals)
            w[w < 0] <- 0
            rwr_scores <- rwr_scores * w
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
        
        return(build_nsea_result(
            base_result = res,
            gene_sets = gene_sets,
            rwr_scores = rwr_scores,
            network = network,
            mode = mode,
            iter = iter,
            p = p,
            organism = "UNKNOWN",
            setType = "UNKNOWN",
            keytype = "UNKNOWN",
            params = list(
                pvalueCutoff = 1.0,
                pAdjustMethod = "BH",
                minGSSize = minGSSize,
                maxGSSize = maxGSSize
            )
        ))
    }
    
    ## ---- corrected whole-pipeline permutation path (default) ----
    scoreType <- if (mode == "evidence") "pos" else "std"
    
    perm_test <- .nsea_permutation_test(
        geneList = geneList, A = A, nodes = nodes, gene_sets = gene_sets,
        mode = mode, p = p, threshold = threshold, maxIter = maxIter,
        specific_weight = specific_weight,
        minGSSize = minGSSize, maxGSSize = maxGSSize,
        nPerm = nPerm, exponent = exponent, scoreType = scoreType,
        seed = seed, verbose = verbose
    )
    if (is.null(perm_test)) {
        return(NULL)
    }
    
    result_df <- .nsea_build_result_df(perm_test, exponent = exponent,
                                       scoreType = scoreType,
                                       pAdjustMethod = "BH")
    
    build_nsea_result(
        base_result = result_df,
        gene_sets = gene_sets,
        rwr_scores = perm_test$ranked,
        network = network,
        mode = mode,
        iter = perm_test$iterations,
        p = p,
        organism = "UNKNOWN",
        setType = "UNKNOWN",
        keytype = "UNKNOWN",
        params = list(
            pvalueCutoff = 1.0,
            pAdjustMethod = "BH",
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            nPerm = perm_test$nPerm,
            seed = seed,
            significance = significance,
            exponent = exponent
        )
    )
}

#' Network-based GSEA using a GSON object
#'
#' @param geneList named numeric vector. In "evidence" mode, must be non-negative. In "signed" mode, can contain both positive and negative values.
#' @param network edge list (data.frame) or sparse matrix.
#' @param gson a GSON object.
#' @param mode character, either "evidence" (default) or "signed".
#' @param p restart probability for RWR (default is 0.5).
#' @param specific_weight logical, whether to apply gene specificity weighting (TF-IDF style) based on gene frequencies in the GSON object. Default is FALSE.
#' @param minGSSize minimal size of each gene set for analyzing. default here is 10.
#' @param maxGSSize maximal size of genes annotated for testing. default here is 500.
#' @param threshold convergence threshold for RWR (default is 1e-9).
#' @param maxIter maximal number of RWR iterations (default is 100).
#' @param significance character, one of "whole_pipeline" (default) or "internal".
#'   "whole_pipeline" computes p-values and NES from a permutation null that
#'   re-runs the full pipeline (permute gene labels -> re-diffuse -> recompute
#'   enrichment scores) and is statistically calibrated. "internal" runs GSEA's
#'   internal label-permutation test on the diffused scores; this is the legacy
#'   behaviour and is **anti-conservative** (it ignores the correlation induced
#'   by network diffusion) — use it only for comparison/debugging.
#' @param nPerm number of whole-pipeline permutations used to estimate p-values
#'   and NES (default: 1000). The smallest estimable p-value is 1/(nPerm + 1).
#' @param seed random seed for the permutation null (default: NULL, use the
#'   current R RNG state). Set a numeric seed for reproducible results.
#' @param exponent weight of each step in the enrichment score (default: 1).
#' @param pvalueCutoff p-value cutoff applied to both raw and adjusted p-values
#'   (default: 0.05). Set to 1 to retain all tested gene sets.
#' @param pAdjustMethod p-value adjustment method (default: "BH").
#' @param verbose logical, print messages.
#' @param ... for `significance = "internal"`, other arguments passed to
#'   `gsea_gson()`; ignored for `significance = "whole_pipeline"`.
#'
#' @return A `nseaResult` object.
#' @export
nsea_gson <- function(geneList,
                      network,
                      gson,
                      mode = c("evidence", "signed"),
                      p = 0.5,
                      specific_weight = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500,
                      threshold = 1e-9,
                      maxIter = 100,
                      significance = c("whole_pipeline", "internal"),
                      nPerm = 1000,
                      seed = NULL,
                      exponent = 1,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = TRUE,
                      ...) {
    
    mode <- match.arg(mode)
    significance <- match.arg(significance)
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("geneList must be a named numeric vector")
    }
    
    if (!inherits(gson, "GSON")) {
        stop("gson should be a GSON object")
    }
    
    if (verbose) message("Preparing network...")
    A <- prepare_network(network)
    nodes <- .nsea_network_nodes(A)
    
    gsid2gene <- gson@gsid2gene
    if (nrow(gsid2gene) == 0) {
        stop("gson has no gene-set annotation (gsid2gene is empty).")
    }
    geneSets <- split(gsid2gene$gene, gsid2gene$gsid)
    
    ## ---- legacy internal-permutation path ----
    if (significance == "internal") {
        if (mode == "evidence" && any(geneList < 0)) {
            warning("geneList contains negative values but mode is 'evidence'. Negative values will be propagated as is, which might violate RWR assumptions. Consider using mode = 'signed'.")
        }
        
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
            rwr_res <- rwr_eigen_cpp(A, v, restart = p, threshold = threshold, max_iter = maxIter)
            rwr_scores <- rwr_res$score
            names(rwr_scores) <- nodes
            iter <- as.integer(rwr_res$iterations)
            
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
            
            iter_up <- 0
            iter_down <- 0
            
            if (sum(v_up) > 0) {
                r_up <- rwr_eigen_cpp(A, v_up, restart = p, threshold = threshold, max_iter = maxIter)
                rwr_up <- r_up$score
                iter_up <- r_up$iterations
            }
            if (sum(v_down) > 0) {
                r_down <- rwr_eigen_cpp(A, v_down, restart = p, threshold = threshold, max_iter = maxIter)
                rwr_down <- r_down$score
                iter_down <- r_down$iterations
            }
            
            rwr_scores <- rwr_up - rwr_down
            names(rwr_scores) <- nodes
            iter <- as.integer(max(iter_up, iter_down))
        }
        
        if (specific_weight) {
            if (verbose) message("Applying gene specificity weighting (TF-IDF)...")
            gsid2gene <- gson@gsid2gene
            N <- length(unique(gsid2gene$gsid))
            gene_freq <- table(gsid2gene$gene)
            df_vals <- as.numeric(gene_freq[names(rwr_scores)])
            df_vals[is.na(df_vals)] <- 1 # Minimum frequency 1
            w <- log(N / df_vals)
            w[w < 0] <- 0
            rwr_scores <- rwr_scores * w
        }
        
        rwr_scores <- sort(rwr_scores, decreasing = TRUE)
        
        if (verbose) message("Running GSEA...")
        if (mode == "evidence") {
            res <- gsea_gson(geneList = rwr_scores,
                             gson = gson,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             scoreType = "pos",
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             ...)
        } else {
            res <- gsea_gson(geneList = rwr_scores,
                             gson = gson,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             scoreType = "std",
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             ...)
        }
        
        return(build_nsea_result(
            base_result = res,
            gene_sets = NULL,
            rwr_scores = rwr_scores,
            network = network,
            mode = mode,
            iter = iter,
            p = p
        ))
    }
    
    ## ---- corrected whole-pipeline permutation path (default) ----
    scoreType <- if (mode == "evidence") "pos" else "std"
    
    perm_test <- .nsea_permutation_test(
        geneList = geneList, A = A, nodes = nodes, gene_sets = geneSets,
        mode = mode, p = p, threshold = threshold, maxIter = maxIter,
        specific_weight = specific_weight,
        minGSSize = minGSSize, maxGSSize = maxGSSize,
        nPerm = nPerm, exponent = exponent, scoreType = scoreType,
        seed = seed, verbose = verbose
    )
    if (is.null(perm_test)) {
        return(NULL)
    }
    
    gsid2name <- gson@gsid2name
    result_df <- .nsea_build_result_df(perm_test, exponent = exponent,
                                       scoreType = scoreType,
                                       gsid2name = gsid2name,
                                       pAdjustMethod = pAdjustMethod)
    
    # filter by pvalueCutoff (both raw and adjusted), matching gsea_gson
    if (!is.null(pvalueCutoff)) {
        result_df <- result_df[!is.na(result_df$pvalue), ]
        result_df <- result_df[result_df$pvalue <= pvalueCutoff, ]
        result_df <- result_df[result_df$p.adjust <= pvalueCutoff, ]
    }
    if (nrow(result_df) == 0) {
        return(NULL)
    }
    
    build_nsea_result(
        base_result = result_df,
        gene_sets = geneSets,
        rwr_scores = perm_test$ranked,
        network = network,
        mode = mode,
        iter = perm_test$iterations,
        p = p,
        params = list(
            pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            nPerm = perm_test$nPerm,
            seed = seed,
            significance = significance,
            exponent = exponent
        )
    )
}
