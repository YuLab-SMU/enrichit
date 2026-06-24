#' Extract pathway subnetwork data from a `mnseaResult`
#'
#' @param res A `mnseaResult` object.
#' @param pathway_id Optional pathway ID. If `NULL`, the top pathway is used.
#' @param include_couplings Logical, whether to include inter-layer coupling
#'   edges. Default is `TRUE`.
#' @param include_isolated Logical, whether to keep nodes without retained
#'   edges. Default is `TRUE`.
#'
#' @return A list with `pathway`, `layer_contribution`, `nodes`, and `edges`.
#' @export
extract_mnsea_subnetwork <- function(res,
                                     pathway_id = NULL,
                                     include_couplings = TRUE,
                                     include_isolated = TRUE) {
    if (!inherits(res, "mnseaResult")) {
        stop("res must be a mnseaResult object.")
    }
    if (nrow(res@result) == 0) {
        stop("The mnseaResult is empty.")
    }
    
    if (is.null(pathway_id)) {
        pathway_id <- as.character(res@result$ID[1])
        message("pathway_id not provided. Using the top pathway: ", pathway_id)
    }
    if (!pathway_id %in% res@result$ID) {
        stop("pathway_id not found in the mnseaResult.")
    }
    
    pathway_row <- res@result[res@result$ID == pathway_id, , drop = FALSE]
    layer_contribution <- res@pathway_contribution[
        res@pathway_contribution$ID == pathway_id,
        ,
        drop = FALSE
    ]
    nodes <- res@feature_contribution[
        res@feature_contribution$ID == pathway_id,
        ,
        drop = FALSE
    ]
    
    if (nrow(nodes) == 0) {
        return(list(
            pathway = pathway_row,
            layer_contribution = layer_contribution,
            nodes = nodes,
            edges = data.frame(
                from = character(0),
                to = character(0),
                from_layer = character(0),
                to_layer = character(0),
                from_feature = character(0),
                to_feature = character(0),
                weight = numeric(0),
                edge_type = character(0),
                stringsAsFactors = FALSE
            )
        ))
    }
    
    nodes$node_key <- paste(nodes$layer, nodes$Feature, sep = "::")
    nodes$collapsed_score <- unname(res@collapsed_scores[nodes$Feature])
    layer_weight_lookup <- res@layer_weights
    nodes$layer_weight <- unname(layer_weight_lookup[nodes$layer])
    
    edge_tables <- list()
    edge_idx <- 1L
    layer_names <- unique(nodes$layer)
    
    for (layer in layer_names) {
        layer_nodes <- unique(nodes$Feature[nodes$layer == layer])
        A <- res@multilayer_network$intra_matrices[[layer]]
        valid_nodes <- intersect(layer_nodes, rownames(A))
        if (length(valid_nodes) == 0) {
            next
        }
        
        idx <- match(valid_nodes, rownames(A))
        idx <- idx[!is.na(idx)]
        if (length(idx) == 0) {
            next
        }
        subA <- A[idx, idx, drop = FALSE]
        sm <- Matrix::summary(subA)
        if (nrow(sm) == 0) {
            next
        }
        
        edge_tables[[edge_idx]] <- data.frame(
            from = paste(layer, rownames(A)[idx][sm$i], sep = "::"),
            to = paste(layer, colnames(A)[idx][sm$j], sep = "::"),
            from_layer = layer,
            to_layer = layer,
            from_feature = rownames(A)[idx][sm$i],
            to_feature = colnames(A)[idx][sm$j],
            weight = sm$x,
            edge_type = "intra",
            stringsAsFactors = FALSE
        )
        edge_idx <- edge_idx + 1L
    }
    
    if (include_couplings && nrow(res@coupling_table) > 0) {
        coupling_edges <- res@coupling_table
        node_keys <- unique(nodes$node_key)
        keep <- paste(coupling_edges$from_layer, coupling_edges$from_id, sep = "::") %in% node_keys &
            paste(coupling_edges$to_layer, coupling_edges$to_id, sep = "::") %in% node_keys
        coupling_edges <- coupling_edges[keep, , drop = FALSE]
        
        if (nrow(coupling_edges) > 0) {
            edge_tables[[edge_idx]] <- data.frame(
                from = paste(coupling_edges$from_layer, coupling_edges$from_id, sep = "::"),
                to = paste(coupling_edges$to_layer, coupling_edges$to_id, sep = "::"),
                from_layer = coupling_edges$from_layer,
                to_layer = coupling_edges$to_layer,
                from_feature = coupling_edges$from_id,
                to_feature = coupling_edges$to_id,
                weight = coupling_edges$weight,
                edge_type = "coupling",
                stringsAsFactors = FALSE
            )
        }
    }
    
    edges <- do.call(rbind, edge_tables)
    if (is.null(edges)) {
        edges <- data.frame(
            from = character(0),
            to = character(0),
            from_layer = character(0),
            to_layer = character(0),
            from_feature = character(0),
            to_feature = character(0),
            weight = numeric(0),
            edge_type = character(0),
            stringsAsFactors = FALSE
        )
    }
    
    if (!include_isolated && nrow(edges) > 0) {
        keep_nodes <- unique(c(edges$from, edges$to))
        nodes <- nodes[nodes$node_key %in% keep_nodes, , drop = FALSE]
    }
    
    list(
        pathway = pathway_row,
        layer_contribution = layer_contribution,
        nodes = nodes,
        edges = edges
    )
}
