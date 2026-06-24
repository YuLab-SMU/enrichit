#' Prepare multi-layer network for repeated propagation
#'
#' @param networks named list of layer-specific networks.
#' @param couplings data.frame of inter-layer edges with columns
#'   `from_layer`, `from_id`, `to_layer`, `to_id`, and optional `weight`.
#' @param directed logical, whether the multi-layer graph is directed.
#' @param intra_normalize one of "column", "row", or "none".
#' @param inter_normalize one of "column", "row", or "none".
#' @param interlayer_strength numeric scalar used to scale all coupling edges.
#' @param layer_order explicit layer order. Defaults to `names(networks)`.
#'
#' @return A `multilayer_network` object.
#' @export
prepare_multilayer_network <- function(networks,
                                       couplings,
                                       directed = FALSE,
                                       intra_normalize = "column",
                                       inter_normalize = "column",
                                       interlayer_strength = 1,
                                       layer_order = names(networks)) {
    intra_normalize <- match.arg(intra_normalize, c("column", "row", "none"))
    inter_normalize <- match.arg(inter_normalize, c("column", "row", "none"))
    
    if (!is.list(networks) || length(networks) == 0) {
        stop("networks must be a non-empty named list.")
    }
    if (is.null(names(networks)) || any(!nzchar(names(networks)))) {
        stop("networks must be a named list with non-empty layer names.")
    }
    if (is.null(layer_order) || !setequal(layer_order, names(networks))) {
        stop("layer_order must contain exactly the same layer names as networks.")
    }
    if (!is.numeric(interlayer_strength) || length(interlayer_strength) != 1 || is.na(interlayer_strength) || interlayer_strength < 0) {
        stop("interlayer_strength must be a single non-negative number.")
    }
    
    layer_order <- as.character(layer_order)
    prepared <- vector("list", length(layer_order))
    names(prepared) <- layer_order
    node_index <- vector("list", length(layer_order))
    offsets <- integer(length(layer_order))
    total_nodes <- 0L
    
    intra_i <- integer(0)
    intra_j <- integer(0)
    intra_x <- numeric(0)
    
    for (k in seq_along(layer_order)) {
        layer <- layer_order[[k]]
        A <- prepare_network(networks[[layer]], directed = directed, normalize = intra_normalize)
        prepared[[layer]] <- A
        nodes <- rownames(A)
        n_nodes <- length(nodes)
        offsets[[k]] <- total_nodes
        
        node_index[[k]] <- data.frame(
            layer = layer,
            node = nodes,
            key = paste(layer, nodes, sep = "::"),
            index = seq_len(n_nodes) + total_nodes,
            stringsAsFactors = FALSE
        )
        
        sm <- Matrix::summary(A)
        if (nrow(sm) > 0) {
            intra_i <- c(intra_i, sm$i + total_nodes)
            intra_j <- c(intra_j, sm$j + total_nodes)
            intra_x <- c(intra_x, sm$x)
        }
        
        total_nodes <- total_nodes + n_nodes
    }
    
    node_index <- do.call(rbind, node_index)
    rownames(node_index) <- NULL
    
    key_to_index <- node_index$index
    names(key_to_index) <- node_index$key
    
    if (missing(couplings) || is.null(couplings)) {
        couplings <- data.frame(
            from_layer = character(0),
            from_id = character(0),
            to_layer = character(0),
            to_id = character(0),
            weight = numeric(0),
            stringsAsFactors = FALSE
        )
    } else {
        couplings <- as.data.frame(couplings, stringsAsFactors = FALSE)
        required_cols <- c("from_layer", "from_id", "to_layer", "to_id")
        if (!all(required_cols %in% colnames(couplings))) {
            stop("couplings must contain columns 'from_layer', 'from_id', 'to_layer', and 'to_id'.")
        }
        if (!"weight" %in% colnames(couplings)) {
            couplings$weight <- 1
        }
        couplings$from_layer <- as.character(couplings$from_layer)
        couplings$from_id <- as.character(couplings$from_id)
        couplings$to_layer <- as.character(couplings$to_layer)
        couplings$to_id <- as.character(couplings$to_id)
        couplings$weight <- as.numeric(couplings$weight)
    }
    
    coupling_matrix <- Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = numeric(0),
        dims = c(total_nodes, total_nodes),
        dimnames = list(node_index$key, node_index$key)
    )
    
    if (nrow(couplings) > 0) {
        from_key <- paste(couplings$from_layer, couplings$from_id, sep = "::")
        to_key <- paste(couplings$to_layer, couplings$to_id, sep = "::")
        
        if (any(!from_key %in% names(key_to_index))) {
            bad <- unique(from_key[!from_key %in% names(key_to_index)])
            stop("Unknown coupling source nodes: ", paste(bad, collapse = ", "))
        }
        if (any(!to_key %in% names(key_to_index))) {
            bad <- unique(to_key[!to_key %in% names(key_to_index)])
            stop("Unknown coupling target nodes: ", paste(bad, collapse = ", "))
        }
        
        ci <- unname(key_to_index[from_key])
        cj <- unname(key_to_index[to_key])
        cx <- couplings$weight
        
        if (!directed) {
            ci_orig <- ci
            cj_orig <- cj
            ci <- c(ci_orig, cj_orig)
            cj <- c(cj_orig, ci_orig)
            cx <- c(cx, cx)
        }
        
        coupling_matrix <- Matrix::sparseMatrix(
            i = ci,
            j = cj,
            x = cx,
            dims = c(total_nodes, total_nodes),
            dimnames = list(node_index$key, node_index$key)
        )
        
        if (inter_normalize == "column") {
            cs <- Matrix::colSums(coupling_matrix)
            cs[cs == 0] <- 1
            coupling_matrix <- coupling_matrix %*% Matrix::Diagonal(x = 1 / cs)
        } else if (inter_normalize == "row") {
            rs <- Matrix::rowSums(coupling_matrix)
            rs[rs == 0] <- 1
            coupling_matrix <- Matrix::Diagonal(x = 1 / rs) %*% coupling_matrix
        }
    }
    
    adjacency <- Matrix::sparseMatrix(
        i = intra_i,
        j = intra_j,
        x = intra_x,
        dims = c(total_nodes, total_nodes),
        dimnames = list(node_index$key, node_index$key)
    )
    adjacency <- adjacency + interlayer_strength * coupling_matrix
    
    structure(
        list(
            adjacency = adjacency,
            intra_matrices = prepared,
            coupling_matrix = coupling_matrix,
            couplings = couplings,
            node_index = node_index,
            layer_order = layer_order,
            directed = directed,
            intra_normalize = intra_normalize,
            inter_normalize = inter_normalize,
            interlayer_strength = interlayer_strength
        ),
        class = "multilayer_network"
    )
}

#' Propagate signals on a multi-layer network
#'
#' @param seed_list named list of named numeric vectors, one per layer.
#' @param network a prepared `multilayer_network` object.
#' @param mode one of "evidence" or "signed".
#' @param p restart probability.
#' @param threshold convergence threshold.
#' @param maxIter maximum number of iterations.
#' @param layer_weights optional named numeric vector of layer weights.
#' @param target_layer optional layer name to focus on downstream.
#'
#' @return A `multilayer_propagation` object.
#' @export
propagate_multilayer <- function(seed_list,
                                 network,
                                 mode = c("evidence", "signed"),
                                 p = 0.5,
                                 threshold = 1e-9,
                                 maxIter = 100,
                                 layer_weights = NULL,
                                 target_layer = NULL) {
    mode <- match.arg(mode)
    
    if (!inherits(network, "multilayer_network")) {
        stop("network must be a multilayer_network object.")
    }
    if (!is.list(seed_list) || length(seed_list) == 0) {
        stop("seed_list must be a non-empty named list.")
    }
    if (is.null(names(seed_list)) || any(!nzchar(names(seed_list)))) {
        stop("seed_list must be a named list with non-empty layer names.")
    }
    
    layer_order <- network$layer_order
    if (!all(names(seed_list) %in% layer_order)) {
        bad <- setdiff(names(seed_list), layer_order)
        stop("Unknown layers in seed_list: ", paste(bad, collapse = ", "))
    }
    
    if (is.null(layer_weights)) {
        layer_weights <- rep(1, length(layer_order))
        names(layer_weights) <- layer_order
    } else {
        if (is.null(names(layer_weights))) {
            if (length(layer_weights) != length(layer_order)) {
                stop("Unnamed layer_weights must match the number of layers.")
            }
            names(layer_weights) <- layer_order
        }
        layer_weights <- layer_weights[layer_order]
        layer_weights[is.na(layer_weights)] <- 1
    }
    
    if (!is.null(target_layer) && !target_layer %in% layer_order) {
        stop("target_layer must be one of the layers in the network.")
    }
    
    node_index <- network$node_index
    node_lookup <- split(node_index, node_index$layer)
    seed_vec <- numeric(nrow(node_index))
    names(seed_vec) <- node_index$key
    
    for (layer in names(seed_list)) {
        s <- seed_list[[layer]]
        if (!is.numeric(s) || is.null(names(s))) {
            stop("Each element in seed_list must be a named numeric vector.")
        }
        idx_df <- node_lookup[[layer]]
        common_nodes <- intersect(names(s), idx_df$node)
        if (length(common_nodes) == 0) {
            next
        }
        match_idx <- idx_df$index[match(common_nodes, idx_df$node)]
        seed_vec[match_idx] <- as.numeric(s[common_nodes]) * layer_weights[[layer]]
    }
    
    if (mode == "evidence" && any(seed_vec < 0)) {
        warning("seed_list contains negative values but mode is 'evidence'. Consider using mode = 'signed'.")
    }
    
    propagation <- .run_multilayer_rwr(
        A = network$adjacency,
        seed = seed_vec,
        mode = mode,
        p = p,
        threshold = threshold,
        maxIter = maxIter
    )
    
    global_scores <- propagation$score
    names(global_scores) <- node_index$key
    layer_scores <- lapply(layer_order, function(layer) {
        idx_df <- node_lookup[[layer]]
        scores <- global_scores[idx_df$key]
        names(scores) <- idx_df$node
        scores
    })
    names(layer_scores) <- layer_order
    
    structure(
        list(
            scores = global_scores,
            layer_scores = layer_scores,
            mode = mode,
            iterations = as.integer(propagation$iterations),
            restart_prob = p,
            layer_weights = layer_weights,
            target_layer = target_layer,
            network = network
        ),
        class = "multilayer_propagation"
    )
}

#' Collapse multi-layer diffusion scores
#'
#' @param x result from `propagate_multilayer()`.
#' @param collapse one of "weighted_mean", "sum", "mean", or "max_abs".
#' @param layer_weights optional named numeric vector used when
#'   `collapse = "weighted_mean"`.
#' @param output_space one of "union" or "gene".
#' @param mapping optional mapping data.frame with `source_id`, `target_id`, and
#'   optional `layer` columns.
#' @param target_layer optional layer name to extract before collapsing.
#'
#' @return A `multilayer_collapsed` object with a `score` vector.
#' @export
collapse_multilayer_scores <- function(x,
                                       collapse = c("weighted_mean", "sum", "mean", "max_abs"),
                                       layer_weights = NULL,
                                       output_space = c("union", "gene"),
                                       mapping = NULL,
                                       target_layer = NULL) {
    collapse <- match.arg(collapse)
    output_space <- match.arg(output_space)
    
    if (!inherits(x, "multilayer_propagation")) {
        stop("x must be a multilayer_propagation object.")
    }
    
    layer_scores <- x$layer_scores
    if (!is.null(target_layer)) {
        if (!target_layer %in% names(layer_scores)) {
            stop("target_layer is not available in the propagation result.")
        }
        layer_scores <- layer_scores[target_layer]
    }
    
    long_df <- do.call(rbind, lapply(names(layer_scores), function(layer) {
        scores <- layer_scores[[layer]]
        data.frame(
            layer = layer,
            source_id = names(scores),
            score = as.numeric(scores),
            stringsAsFactors = FALSE
        )
    }))
    
    if (nrow(long_df) == 0) {
        stop("No layer scores available to collapse.")
    }
    
    if (is.null(mapping)) {
        long_df$target_id <- long_df$source_id
    } else {
        mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
        if (!all(c("source_id", "target_id") %in% colnames(mapping))) {
            stop("mapping must contain 'source_id' and 'target_id' columns.")
        }
        
        if ("layer" %in% colnames(mapping)) {
            merged <- merge(long_df, mapping, by = c("layer", "source_id"), all.x = TRUE, sort = FALSE)
        } else {
            merged <- merge(long_df, mapping, by = "source_id", all.x = TRUE, sort = FALSE)
        }
        merged$target_id[is.na(merged$target_id)] <- merged$source_id[is.na(merged$target_id)]
        long_df <- merged
    }
    
    feature_ids <- unique(long_df$target_id)
    layer_names <- unique(long_df$layer)
    score_mat <- matrix(
        NA_real_,
        nrow = length(feature_ids),
        ncol = length(layer_names),
        dimnames = list(feature_ids, layer_names)
    )
    
    for (feature in feature_ids) {
        for (layer in layer_names) {
            vals <- long_df$score[long_df$target_id == feature & long_df$layer == layer]
            if (length(vals) == 0) {
                next
            }
            score_mat[feature, layer] <- .collapse_vector(vals, method = collapse)
        }
    }
    
    if (is.null(layer_weights)) {
        layer_weights <- x$layer_weights[layer_names]
    } else {
        if (is.null(names(layer_weights))) {
            if (length(layer_weights) != length(layer_names)) {
                stop("Unnamed layer_weights must match the number of layers used for collapsing.")
            }
            names(layer_weights) <- layer_names
        }
        layer_weights <- layer_weights[layer_names]
    }
    layer_weights[is.na(layer_weights)] <- 1
    
    collapsed <- apply(score_mat, 1, function(s) {
        .collapse_vector(s, method = collapse, weights = layer_weights)
    })
    collapsed <- collapsed[!is.na(collapsed)]
    
    structure(
        list(
            score = collapsed,
            score_matrix = score_mat[rownames(score_mat) %in% names(collapsed), , drop = FALSE],
            feature_id = names(collapsed),
            output_space = output_space,
            collapse_method = collapse,
            layer_weights = layer_weights,
            target_layer = target_layer
        ),
        class = "multilayer_collapsed"
    )
}

#' Multi-layer Network-based Gene Set Enrichment Analysis
#'
#' @param seed_list named list of named numeric vectors, one per layer.
#' @param networks named list of layer-specific networks.
#' @param couplings data.frame of inter-layer edges.
#' @param gene_sets list of gene sets.
#' @param mode one of "evidence" or "signed".
#' @param layer_weights optional named numeric vector.
#' @param collapse one of "weighted_mean", "sum", "mean", or "max_abs".
#' @param target_layer optional layer name to export scores from.
#' @param output_space one of "union" or "gene".
#' @param p restart probability.
#' @param interlayer_strength global scaling factor for coupling edges.
#' @param specific_weight logical.
#' @param minGSSize minimal size of each gene set.
#' @param maxGSSize maximal size of genes annotated for testing.
#' @param threshold convergence threshold.
#' @param maxIter maximal number of iterations.
#' @param verbose logical.
#' @param ... additional arguments passed to `gsea()`.
#'
#' @return A `mnseaResult` object.
#' @export
mnsea <- function(seed_list,
                  networks,
                  couplings,
                  gene_sets,
                  mode = c("evidence", "signed"),
                  layer_weights = NULL,
                  collapse = c("weighted_mean", "sum", "mean", "max_abs"),
                  target_layer = NULL,
                  output_space = c("union", "gene"),
                  p = 0.5,
                  interlayer_strength = 1,
                  specific_weight = FALSE,
                  minGSSize = 10,
                  maxGSSize = 500,
                  threshold = 1e-9,
                  maxIter = 100,
                  verbose = TRUE,
                  ...) {
    mode <- match.arg(mode)
    collapse <- match.arg(collapse)
    output_space <- match.arg(output_space)
    
    if (verbose) message("Preparing multi-layer network...")
    ml_net <- prepare_multilayer_network(
        networks = networks,
        couplings = couplings,
        interlayer_strength = interlayer_strength
    )
    
    if (verbose) message("Running multi-layer propagation...")
    prop <- propagate_multilayer(
        seed_list = seed_list,
        network = ml_net,
        mode = mode,
        p = p,
        threshold = threshold,
        maxIter = maxIter,
        layer_weights = layer_weights,
        target_layer = target_layer
    )
    
    collapsed <- collapse_multilayer_scores(
        x = prop,
        collapse = collapse,
        layer_weights = layer_weights,
        output_space = output_space,
        target_layer = target_layer
    )
    
    collapsed_scores <- collapsed$score
    if (length(collapsed_scores) == 0) {
        stop("No collapsed scores available for downstream enrichment.")
    }
    if (specific_weight) {
        if (verbose) message("Applying gene specificity weighting (TF-IDF)...")
        N <- length(gene_sets)
        gene_freq <- table(unlist(gene_sets))
        df_vals <- as.numeric(gene_freq[names(collapsed_scores)])
        df_vals[is.na(df_vals)] <- 1
        w <- log(N / df_vals)
        w[w < 0] <- 0
        collapsed_scores <- collapsed_scores * w
    }
    collapsed_scores <- sort(collapsed_scores, decreasing = TRUE)
    
    if (verbose) message("Running GSEA...")
    if (mode == "evidence") {
        res <- gsea(
            geneList = collapsed_scores,
            gene_sets = gene_sets,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            scoreType = "pos",
            ...
        )
    } else {
        res <- gsea(
            geneList = collapsed_scores,
            gene_sets = gene_sets,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            scoreType = "std",
            ...
        )
    }
    if (is.null(res)) {
        return(NULL)
    }
    
    build_mnsea_result(
        base_result = res,
        gene_sets = gene_sets,
        propagation = prop,
        collapsed = collapsed_scores,
        collapse_method = collapse,
        output_space = output_space,
        p = p,
        params = list(
            pvalueCutoff = 1.0,
            pAdjustMethod = "BH",
            minGSSize = minGSSize,
            maxGSSize = maxGSSize
        )
    )
}

#' Multi-layer NSEA using a GSON object
#'
#' @param seed_list named list of named numeric vectors, one per layer.
#' @param networks named list of layer-specific networks.
#' @param couplings data.frame of inter-layer edges.
#' @param gson a GSON object.
#' @param mode one of "evidence" or "signed".
#' @param layer_weights optional named numeric vector.
#' @param collapse one of "weighted_mean", "sum", "mean", or "max_abs".
#' @param target_layer optional layer name to export scores from.
#' @param output_space one of "union" or "gene".
#' @param p restart probability.
#' @param interlayer_strength global scaling factor for coupling edges.
#' @param specific_weight logical.
#' @param minGSSize minimal size of each gene set.
#' @param maxGSSize maximal size of genes annotated for testing.
#' @param threshold convergence threshold.
#' @param maxIter maximal number of iterations.
#' @param verbose logical.
#' @param ... additional arguments passed to `gsea_gson()`.
#'
#' @return A `mnseaResult` object.
#' @export
mnsea_gson <- function(seed_list,
                       networks,
                       couplings,
                       gson,
                       mode = c("evidence", "signed"),
                       layer_weights = NULL,
                       collapse = c("weighted_mean", "sum", "mean", "max_abs"),
                       target_layer = NULL,
                       output_space = c("union", "gene"),
                       p = 0.5,
                       interlayer_strength = 1,
                       specific_weight = FALSE,
                       minGSSize = 10,
                       maxGSSize = 500,
                       threshold = 1e-9,
                       maxIter = 100,
                       verbose = TRUE,
                       ...) {
    mode <- match.arg(mode)
    collapse <- match.arg(collapse)
    output_space <- match.arg(output_space)
    
    if (verbose) message("Preparing multi-layer network...")
    ml_net <- prepare_multilayer_network(
        networks = networks,
        couplings = couplings,
        interlayer_strength = interlayer_strength
    )
    
    if (verbose) message("Running multi-layer propagation...")
    prop <- propagate_multilayer(
        seed_list = seed_list,
        network = ml_net,
        mode = mode,
        p = p,
        threshold = threshold,
        maxIter = maxIter,
        layer_weights = layer_weights,
        target_layer = target_layer
    )
    
    collapsed <- collapse_multilayer_scores(
        x = prop,
        collapse = collapse,
        layer_weights = layer_weights,
        output_space = output_space,
        target_layer = target_layer
    )
    
    collapsed_scores <- collapsed$score
    if (length(collapsed_scores) == 0) {
        stop("No collapsed scores available for downstream enrichment.")
    }
    if (specific_weight) {
        if (verbose) message("Applying gene specificity weighting (TF-IDF)...")
        gsid2gene <- gson@gsid2gene
        N <- length(unique(gsid2gene$gsid))
        gene_freq <- table(gsid2gene$gene)
        df_vals <- as.numeric(gene_freq[names(collapsed_scores)])
        df_vals[is.na(df_vals)] <- 1
        w <- log(N / df_vals)
        w[w < 0] <- 0
        collapsed_scores <- collapsed_scores * w
    }
    collapsed_scores <- sort(collapsed_scores, decreasing = TRUE)
    
    if (verbose) message("Running GSEA...")
    if (mode == "evidence") {
        res <- gsea_gson(
            geneList = collapsed_scores,
            gson = gson,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            scoreType = "pos",
            ...
        )
    } else {
        res <- gsea_gson(
            geneList = collapsed_scores,
            gson = gson,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            scoreType = "std",
            ...
        )
    }
    if (is.null(res)) {
        return(NULL)
    }
    
    build_mnsea_result(
        base_result = res,
        gene_sets = NULL,
        propagation = prop,
        collapsed = collapsed_scores,
        collapse_method = collapse,
        output_space = output_space,
        p = p,
        params = NULL
    )
}

build_mnsea_result <- function(base_result,
                               gene_sets,
                               propagation,
                               collapsed,
                               collapse_method,
                               output_space,
                               p,
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
        gene_list <- collapsed
        perm_scores <- matrix(0, nrow = 0, ncol = 0)
        gene2symbol <- character(0)
        readable <- FALSE
        termsim <- matrix(0, nrow = 0, ncol = 0)
        method <- "MNSEA"
        dr <- list()
        organism <- "UNKNOWN"
        setType <- "UNKNOWN"
        keytype <- "UNKNOWN"
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
    
    explain_cache <- .build_mnsea_explain_cache(
        result_df = result_df,
        gene_sets = gene_sets,
        layer_scores = propagation$layer_scores
    )
    
    new("mnseaResult",
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
        multilayer_network = propagation$network,
        layer_scores = propagation$layer_scores,
        collapsed_scores = collapsed,
        layer_weights = propagation$layer_weights,
        coupling_table = propagation$network$couplings,
        mode = propagation$mode,
        iterations = as.integer(propagation$iterations),
        restart_prob = p,
        collapse_method = collapse_method,
        target_layer = propagation$target_layer %||% "",
        output_space = output_space,
        pathway_contribution = explain_cache$pathway_contribution,
        feature_contribution = explain_cache$feature_contribution)
}

.run_multilayer_rwr <- function(A, seed, mode, p, threshold, maxIter) {
    if (mode == "evidence") {
        sum_seed <- sum(seed)
        if (sum_seed <= 0) {
            stop("The sum of seed scores in the multi-layer network is zero.")
        }
        v <- seed / sum_seed
        rwr_res <- rwr_eigen_cpp(A, v, restart = p, threshold = threshold, max_iter = maxIter)
        return(list(score = rwr_res$score, iterations = rwr_res$iterations))
    }
    
    v_up <- ifelse(seed > 0, seed, 0)
    v_down <- ifelse(seed < 0, abs(seed), 0)
    if (sum(v_up) > 0) v_up <- v_up / sum(v_up)
    if (sum(v_down) > 0) v_down <- v_down / sum(v_down)
    
    rwr_up <- rep(0, length(seed))
    rwr_down <- rep(0, length(seed))
    iter_up <- 0
    iter_down <- 0
    
    if (sum(v_up) > 0) {
        up_res <- rwr_eigen_cpp(A, v_up, restart = p, threshold = threshold, max_iter = maxIter)
        rwr_up <- up_res$score
        iter_up <- up_res$iterations
    }
    if (sum(v_down) > 0) {
        down_res <- rwr_eigen_cpp(A, v_down, restart = p, threshold = threshold, max_iter = maxIter)
        rwr_down <- down_res$score
        iter_down <- down_res$iterations
    }
    
    list(score = rwr_up - rwr_down, iterations = max(iter_up, iter_down))
}

.collapse_vector <- function(x, method, weights = NULL) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
        return(NA_real_)
    }
    if (method == "sum") {
        return(sum(x))
    }
    if (method == "mean") {
        return(mean(x))
    }
    if (method == "max_abs") {
        return(x[which.max(abs(x))])
    }
    if (method == "weighted_mean") {
        if (is.null(weights)) {
            return(mean(x))
        }
        w <- weights[names(weights) %in% names(x)]
        if (length(w) == 0) {
            return(mean(x))
        }
        return(sum(x * w) / sum(w))
    }
    stop("Unsupported collapse method: ", method)
}

.build_mnsea_explain_cache <- function(result_df, gene_sets, layer_scores) {
    layer_names <- names(layer_scores)
    if (length(layer_names) == 0 || nrow(result_df) == 0) {
        empty_pathway <- data.frame(
            ID = character(0),
            Description = character(0),
            layer = character(0),
            contribution = numeric(0),
            share = numeric(0),
            n_feature = integer(0),
            stringsAsFactors = FALSE
        )
        empty_feature <- data.frame(
            ID = character(0),
            Description = character(0),
            Feature = character(0),
            layer = character(0),
            score = numeric(0),
            abs_score = numeric(0),
            is_core = logical(0),
            stringsAsFactors = FALSE
        )
        return(list(
            pathway_contribution = empty_pathway,
            feature_contribution = empty_feature
        ))
    }
    
    layer_feature_union <- unique(unlist(lapply(layer_scores, names), use.names = FALSE))
    pathway_rows <- vector("list", nrow(result_df))
    feature_rows <- vector("list", nrow(result_df))
    
    for (i in seq_len(nrow(result_df))) {
        path_id <- as.character(result_df$ID[i])
        path_desc <- as.character(result_df$Description[i])
        gene_col <- if ("core_enrichment" %in% colnames(result_df) &&
            !is.na(result_df$core_enrichment[i]) &&
            nzchar(as.character(result_df$core_enrichment[i]))) {
            "core_enrichment"
        } else if ("geneID" %in% colnames(result_df) &&
            !is.na(result_df$geneID[i]) &&
            nzchar(as.character(result_df$geneID[i]))) {
            "geneID"
        } else {
            NA_character_
        }
        
        if (!is.na(gene_col)) {
            selected_features <- unique(unlist(strsplit(as.character(result_df[[gene_col]][i]), "/")))
        } else if (!is.null(gene_sets) && path_id %in% names(gene_sets)) {
            selected_features <- unique(as.character(gene_sets[[path_id]]))
        } else {
            selected_features <- character(0)
        }
        
        selected_features <- intersect(selected_features, layer_feature_union)
        if (length(selected_features) == 0) {
            next
        }
        
        if ("core_enrichment" %in% colnames(result_df) &&
            !is.na(result_df$core_enrichment[i]) &&
            nzchar(as.character(result_df$core_enrichment[i]))) {
            core_features <- unique(unlist(strsplit(as.character(result_df$core_enrichment[i]), "/")))
        } else {
            core_features <- selected_features
        }
        
        layer_contrib <- numeric(length(layer_names))
        names(layer_contrib) <- layer_names
        per_layer_rows <- vector("list", length(layer_names))
        
        for (j in seq_along(layer_names)) {
            layer <- layer_names[[j]]
            scores <- layer_scores[[layer]]
            valid_features <- intersect(selected_features, names(scores))
            if (length(valid_features) == 0) {
                per_layer_rows[[j]] <- NULL
                layer_contrib[[layer]] <- 0
                next
            }
            
            layer_vals <- scores[valid_features]
            layer_contrib[[layer]] <- mean(abs(layer_vals))
            per_layer_rows[[j]] <- data.frame(
                ID = path_id,
                Description = path_desc,
                Feature = valid_features,
                layer = layer,
                score = as.numeric(layer_vals),
                abs_score = abs(as.numeric(layer_vals)),
                is_core = valid_features %in% core_features,
                stringsAsFactors = FALSE
            )
        }
        
        share <- if (sum(layer_contrib) > 0) layer_contrib / sum(layer_contrib) else rep(0, length(layer_contrib))
        pathway_rows[[i]] <- data.frame(
            ID = path_id,
            Description = path_desc,
            layer = layer_names,
            contribution = as.numeric(layer_contrib),
            share = as.numeric(share),
            n_feature = length(selected_features),
            stringsAsFactors = FALSE
        )
        feature_rows[[i]] <- do.call(rbind, per_layer_rows)
    }
    
    pathway_contribution <- do.call(rbind, pathway_rows)
    feature_contribution <- do.call(rbind, feature_rows)
    
    if (is.null(pathway_contribution)) {
        pathway_contribution <- data.frame(
            ID = character(0),
            Description = character(0),
            layer = character(0),
            contribution = numeric(0),
            share = numeric(0),
            n_feature = integer(0),
            stringsAsFactors = FALSE
        )
    }
    if (is.null(feature_contribution)) {
        feature_contribution <- data.frame(
            ID = character(0),
            Description = character(0),
            Feature = character(0),
            layer = character(0),
            score = numeric(0),
            abs_score = numeric(0),
            is_core = logical(0),
            stringsAsFactors = FALSE
        )
    }
    
    list(
        pathway_contribution = pathway_contribution,
        feature_contribution = feature_contribution
    )
}
