#' Aggregate multi-omics gene/protein-level statistics
#'
#' Aggregate multi-omics or multi-source statistics into a unified object for downstream enrichment analysis.
#'
#' @param x A list of named numeric vectors, a data.frame, or a matrix. Row names (or names for vectors) must represent feature IDs.
#' @param method Character, aggregation method. One of "fisher", "stouffer", "mean", or "max_abs".
#' @param input Character, input type. One of "pvalue" or "signed_score".
#' @param feature_type Character, type of the features (e.g., "gene", "protein"). Default is "gene".
#' @param conflict_policy Character, strategy to handle directional conflicts when input is "signed_score". One of "keep_all" (default, ignore conflicts), "strict" (set to NA if any signs conflict), or "penalty" (divide final score by 2 if signs conflict).
#' @param ... Additional arguments.
#'
#' @return An object of class `omics_aggregated` containing `score`, `pvalue` (if input is "pvalue"), `input_type`, `feature_type`, and `feature_id`.
#' @export
#' @importFrom stats pchisq pnorm qnorm
aggregate_omics <- function(x, method = c("fisher", "stouffer", "brown", "mean", "weighted_mean", "max_abs"), 
                            input = c("pvalue", "signed_score"), feature_type = "gene", 
                            conflict_policy = c("keep_all", "strict", "penalty"), ...) {
    method <- match.arg(method)
    input <- match.arg(input)
    conflict_policy <- match.arg(conflict_policy)
    
    # Convert input to a matrix of features x omics
    if (is.list(x) && !is.data.frame(x)) {
        all_features <- unique(unlist(lapply(x, names)))
        mat <- matrix(NA_real_, nrow = length(all_features), ncol = length(x))
        rownames(mat) <- all_features
        colnames(mat) <- names(x)
        for (i in seq_along(x)) {
            mat[names(x[[i]]), i] <- as.numeric(x[[i]])
        }
    } else if (is.data.frame(x) || is.matrix(x)) {
        mat <- as.matrix(x)
        if (is.null(rownames(mat))) stop("Input must have row names representing feature IDs")
    } else {
        stop("x must be a list of named vectors, a data.frame, or a matrix")
    }
    
    res <- list(
        input_type = input,
        feature_type = feature_type,
        feature_id = rownames(mat)
    )
    
    if (input == "pvalue") {
        if (method == "fisher") {
            # Fisher's method: -2 * sum(ln(p))
            agg_p <- apply(mat, 1, function(p) {
                p <- p[!is.na(p)]
                if (length(p) == 0) return(NA_real_)
                p <- pmax(p, .Machine$double.xmin)
                stat <- -2 * sum(log(p))
                pchisq(stat, df = 2 * length(p), lower.tail = FALSE)
            })
        } else if (method == "stouffer") {
            # Stouffer's method: sum(Z) / sqrt(k)
            agg_p <- apply(mat, 1, function(p) {
                p <- p[!is.na(p)]
                if (length(p) == 0) return(NA_real_)
                p <- pmax(p, .Machine$double.xmin)
                p <- pmin(p, 1 - .Machine$double.eps)
                z <- qnorm(p, lower.tail = FALSE)
                z_stouffer <- sum(z) / sqrt(length(z))
                pnorm(z_stouffer, lower.tail = FALSE)
            })
        } else if (method == "brown") {
            # Brown's method (Empirical Brown's Method)
            args <- list(...)
            if (!is.null(args$cov_matrix)) {
                cov_mat <- args$cov_matrix
            } else {
                # Compute empirical covariance of -2 * log(p) across all features
                logP <- -2 * log(pmax(mat, .Machine$double.xmin))
                cov_mat <- stats::cov(logP, use = "pairwise.complete.obs")
                cov_mat[is.na(cov_mat)] <- 0
            }
            
            agg_p <- apply(mat, 1, function(p) {
                idx <- which(!is.na(p))
                k <- length(idx)
                if (k == 0) return(NA_real_)
                if (k == 1) return(p[idx])
                
                p_vals <- pmax(p[idx], .Machine$double.xmin)
                X <- -2 * sum(log(p_vals))
                
                VX <- sum(cov_mat[idx, idx])
                
                # If empirical variance is <= 4k, fallback to Fisher's independent assumption
                if (is.na(VX) || VX <= 4 * k) {
                    return(pchisq(X, df = 2 * k, lower.tail = FALSE))
                }
                
                c_factor <- VX / (4 * k)
                df <- (8 * k^2) / VX
                
                pchisq(X / c_factor, df = df, lower.tail = FALSE)
            })
        } else {
            stop("For input='pvalue', method must be 'fisher', 'stouffer', or 'brown'")
        }
        res$pvalue <- agg_p
        # Convert aggregated pvalue to score for ranking
        res$score <- -log10(pmax(agg_p, .Machine$double.xmin))
        
    } else if (input == "signed_score") {
        
        has_conflict <- apply(mat, 1, function(s) {
            s <- s[!is.na(s) & s != 0]
            if (length(s) <= 1) return(FALSE)
            return(length(unique(sign(s))) > 1)
        })
        
        if (method == "mean") {
            agg_s <- rowMeans(mat, na.rm = TRUE)
        } else if (method == "weighted_mean") {
            args <- list(...)
            weights <- args$weights
            if (is.null(weights)) {
                warning("method='weighted_mean' requires 'weights' argument. Falling back to simple mean.")
                agg_s <- rowMeans(mat, na.rm = TRUE)
            } else {
                if (length(weights) != ncol(mat)) {
                    stop("Length of 'weights' must match the number of omics columns.")
                }
                agg_s <- apply(mat, 1, function(s) {
                    idx <- !is.na(s)
                    if (sum(idx) == 0) return(NA_real_)
                    sum(s[idx] * weights[idx]) / sum(weights[idx])
                })
            }
        } else if (method == "max_abs") {
            agg_s <- apply(mat, 1, function(s) {
                s <- s[!is.na(s)]
                if (length(s) == 0) return(NA_real_)
                s[which.max(abs(s))]
            })
        } else {
            stop("For input='signed_score', method must be 'mean', 'weighted_mean', or 'max_abs'")
        }
        
        if (conflict_policy == "strict") {
            agg_s[has_conflict] <- NA_real_
        } else if (conflict_policy == "penalty") {
            agg_s[has_conflict] <- agg_s[has_conflict] / 2
        }
        
        res$score <- agg_s
    }
    
    # Remove features with NA scores
    valid_idx <- !is.na(res$score)
    res$feature_id <- res$feature_id[valid_idx]
    res$score <- res$score[valid_idx]
    names(res$score) <- res$feature_id
    
    if (!is.null(res$pvalue)) {
        res$pvalue <- res$pvalue[valid_idx]
        names(res$pvalue) <- res$feature_id
    }
    
    res$original_matrix <- mat[valid_idx, , drop = FALSE]
    
    class(res) <- "omics_aggregated"
    return(res)
}

#' Harmonize feature IDs to a target space
#'
#' Map protein-level or other feature-level statistics to a unified gene-level space.
#'
#' @param x A structured result from `aggregate_omics()`.
#' @param mapping A data.frame with `source_id` and `target_id` columns.
#' @param from Character, source feature type. Default is "protein".
#' @param to Character, target feature type. Default is "gene".
#' @param collapse Character, method to collapse multiple source IDs mapped to a single target ID. One of "max_abs", "mean", or "min_p".
#'
#' @return A harmonized `omics_aggregated` object.
#' @export
harmonize_ids <- function(x, mapping, from = "protein", to = "gene", collapse = c("max_abs", "mean", "min_p")) {
    collapse <- match.arg(collapse)
    if (!inherits(x, "omics_aggregated")) {
        stop("x must be an omics_aggregated object")
    }
    
    if (!all(c("source_id", "target_id") %in% colnames(mapping))) {
        stop("mapping must be a data.frame with 'source_id' and 'target_id' columns")
    }
    
    df <- data.frame(source_id = x$feature_id, score = x$score, stringsAsFactors = FALSE)
    if (!is.null(x$pvalue)) {
        df$pvalue <- x$pvalue
    }
    
    merged <- merge(df, mapping, by = "source_id", all.x = FALSE)
    if (nrow(merged) == 0) {
        warning("No IDs mapped successfully.")
        return(NULL)
    }
    
    res_list <- split(merged, merged$target_id)
    target_ids <- names(res_list)
    scores <- numeric(length(target_ids))
    pvals <- if (!is.null(x$pvalue)) numeric(length(target_ids)) else NULL
    
    for (i in seq_along(res_list)) {
        sub_df <- res_list[[i]]
        if (collapse == "max_abs") {
            idx <- which.max(abs(sub_df$score))
            scores[i] <- sub_df$score[idx]
            if (!is.null(pvals)) pvals[i] <- sub_df$pvalue[idx]
        } else if (collapse == "mean") {
            scores[i] <- mean(sub_df$score, na.rm = TRUE)
            if (!is.null(pvals)) pvals[i] <- mean(sub_df$pvalue, na.rm = TRUE)
        } else if (collapse == "min_p") {
            if (is.null(pvals)) stop("collapse='min_p' requires pvalue input in the aggregated object")
            idx <- which.min(sub_df$pvalue)
            scores[i] <- sub_df$score[idx]
            pvals[i] <- sub_df$pvalue[idx]
        }
    }
    
    names(scores) <- target_ids
    if (!is.null(pvals)) names(pvals) <- target_ids
    
    res <- list(
        input_type = x$input_type,
        feature_type = to,
        feature_id = target_ids,
        score = scores
    )
    if (!is.null(pvals)) res$pvalue <- pvals
    
    class(res) <- "omics_aggregated"
    return(res)
}

#' Select features for ORA
#'
#' Convert continuous aggregated statistics into a discrete list of genes and a universe for Over-Representation Analysis.
#'
#' @param x A structured result from `aggregate_omics()` or `harmonize_ids()`.
#' @param cutoff Numeric, the threshold to apply.
#' @param by Character, metric to apply the threshold on. One of "pvalue" or "score".
#' @param ... Additional arguments.
#'
#' @return A list containing `gene` (the selected feature IDs) and `universe` (all feature IDs).
#' @export
select_features_for_ora <- function(x, cutoff = 0.05, by = c("pvalue", "score"), ...) {
    by <- match.arg(by)
    if (!inherits(x, "omics_aggregated")) {
        stop("x must be an omics_aggregated object")
    }
    
    if (by == "pvalue") {
        if (is.null(x$pvalue)) stop("x does not contain pvalue. Use by='score' instead.")
        gene <- x$feature_id[x$pvalue < cutoff]
    } else if (by == "score") {
        gene <- x$feature_id[x$score > cutoff]
    }
    
    list(
        gene = gene,
        universe = x$feature_id
    )
}

#' Aggregate multiple enrichment results (Late Fusion)
#'
#' Combine pathway-level enrichment results from multiple omics or independent analyses.
#' P-values of identical pathways are merged using statistical methods (e.g., Brown's method).
#'
#' @param res_list A named list of enrichment result objects (e.g., `enrichResult`, `gseaResult`, `nseaResult`).
#' @param method Character, aggregation method for p-values. One of "brown", "fisher", or "stouffer".
#' @param ... Additional arguments passed to `aggregate_omics` (e.g., `cov_matrix` for Brown's method).
#'
#' @return An `enrichResult` object containing the aggregated p-values, FDR, and combined gene lists.
#' @export
aggregate_enrichment <- function(res_list, method = c("brown", "fisher", "stouffer"), ...) {
    method <- match.arg(method)
    
    if (!is.list(res_list) || length(res_list) < 2) {
        stop("res_list must be a list of at least two enrichment result objects.")
    }
    
    if (is.null(names(res_list))) {
        names(res_list) <- paste0("Omics_", seq_along(res_list))
    }
    
    # Extract results to data.frames
    df_list <- lapply(res_list, function(x) {
        if (inherits(x, c("enrichResult", "gseaResult", "nseaResult"))) {
            return(as.data.frame(x))
        } else if (is.data.frame(x)) {
            return(x)
        } else {
            stop("Elements in res_list must be enrichResult, gseaResult, nseaResult, or data.frame.")
        }
    })
    
    # Get all unique pathway IDs
    all_ids <- unique(unlist(lapply(df_list, function(df) df$ID)))
    
    # Build p-value matrix (rows = pathways, cols = omics)
    p_mat <- matrix(NA_real_, nrow = length(all_ids), ncol = length(df_list))
    rownames(p_mat) <- all_ids
    colnames(p_mat) <- names(res_list)
    
    # Build description map and combined gene lists
    desc_map <- character(length(all_ids))
    names(desc_map) <- all_ids
    
    gene_list_map <- lapply(all_ids, function(id) character(0))
    names(gene_list_map) <- all_ids
    
    for (i in seq_along(df_list)) {
        df <- df_list[[i]]
        if (nrow(df) == 0) next
        
        idx <- match(df$ID, all_ids)
        p_mat[idx, i] <- df$pvalue
        
        # Map descriptions (take the first non-NA/non-empty encountered)
        unmapped <- desc_map[idx] == "" | is.na(desc_map[idx])
        if (any(unmapped)) {
            desc_map[idx[unmapped]] <- df$Description[unmapped]
        }
        
        # Combine genes from 'geneID' (ORA) or 'core_enrichment' (GSEA)
        gene_col <- if ("geneID" %in% colnames(df)) "geneID" else if ("core_enrichment" %in% colnames(df)) "core_enrichment" else NULL
        if (!is.null(gene_col)) {
            for (j in seq_len(nrow(df))) {
                id <- df$ID[j]
                genes <- unlist(strsplit(as.character(df[[gene_col]][j]), "/"))
                gene_list_map[[id]] <- unique(c(gene_list_map[[id]], genes))
            }
        }
    }
    
    # Aggregate p-values using the underlying feature-level aggregator
    # This perfectly reuses Brown's/Fisher's/Stouffer's logic on the pathway level!
    agg_res <- aggregate_omics(p_mat, method = method, input = "pvalue", ...)
    combined_p <- agg_res$pvalue
    
    # Drop pathways with NA aggregated p-value
    valid_idx <- !is.na(combined_p)
    combined_p <- combined_p[valid_idx]
    all_ids <- names(combined_p)
    
    padj <- stats::p.adjust(combined_p, method = "BH")
    
    # Safely compute qvalue if function exists, else NA
    qval <- tryCatch({
        calculate_qvalue(combined_p)
    }, error = function(e) rep(NA_real_, length(combined_p)))
    
    combined_genes <- vapply(all_ids, function(id) paste(gene_list_map[[id]], collapse = "/"), character(1))
    counts <- vapply(all_ids, function(id) length(gene_list_map[[id]]), integer(1))
    
    res_df <- data.frame(
        ID = all_ids,
        Description = desc_map[all_ids],
        pvalue = combined_p,
        p.adjust = padj,
        qvalue = qval,
        geneID = combined_genes,
        Count = counts,
        stringsAsFactors = FALSE
    )
    
    res_df <- res_df[order(res_df$pvalue), ]
    rownames(res_df) <- res_df$ID
    
    # Create a generic enrichResult to hold the late fusion output
    methods::new("enrichResult",
        result = res_df,
        pvalueCutoff = 1.0,
        pAdjustMethod = "BH",
        qvalueCutoff = 1.0,
        gene = character(0),
        universe = character(0),
        geneSets = list(),
        organism = "UNKNOWN",
        keytype = "UNKNOWN",
        ontology = "Multi-omics Late Fusion",
        readable = FALSE
    )
}
