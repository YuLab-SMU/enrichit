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
aggregate_omics <- function(x, method = c("fisher", "stouffer", "mean", "max_abs"), 
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
        } else {
            stop("For input='pvalue', method must be 'fisher' or 'stouffer'")
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
        } else if (method == "max_abs") {
            agg_s <- apply(mat, 1, function(s) {
                s <- s[!is.na(s)]
                if (length(s) == 0) return(NA_real_)
                s[which.max(abs(s))]
            })
        } else {
            stop("For input='signed_score', method must be 'mean' or 'max_abs'")
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
