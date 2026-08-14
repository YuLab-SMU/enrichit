#' Gene Set Enrichment Analysis (GSEA)
#'
#' Perform Gene Set Enrichment Analysis (GSEA) using a ranked gene list.
#'
#' @inheritParams enrichit_params
#' @param weight A named numeric vector of weights for genes. The names should match the names of geneList. If provided, the geneList will be multiplied by the weight and resorted before GSEA (default: NULL).
#' @param eps Epsilon for multilevel methods (default: 1e-10). Sets the smallest p-value that can be estimated.
#' @param sampleSize Sample size for multilevel methods (default: 101).
#' @param seed Random seed for reproducibility (default: FALSE). Set to a number to make
#'   the result reproducible; set to TRUE to use a fixed default seed. If FALSE (default),
#'   a seed is drawn from R's random number generator on each run, so results may vary
#'   between runs. Note that the C++ engine uses its own RNG seeded with this value:
#'   passing an explicit numeric seed (or calling \code{set.seed()} before \code{gsea()}
#'   when \code{seed = FALSE}) makes the result reproducible.
#' @param nPermSimple Number of permutations for the simple method (default: 1000).
#' @param scoreType Type of enrichment score calculation: "std", "pos", "neg" (default: "std").
#'
#' @return A data.frame with columns:
#' - **ID**: Gene set name
#' - **enrichmentScore**: Enrichment Score
#' - **NES**: Normalized Enrichment Score
#' - **pvalue**: Empirical p-value from permutation test
#' - **setSize**: Size of the gene set (number of genes found in geneList)
#' - **nPerm**: (adaptive mode only) Actual number of permutations used
#' - **rank**: Rank at which the maximum enrichment score is attained
#' - **leading_edge**: Leading edge statistics (tags, list, signal)
#' - **core_enrichment**: Genes in the leading edge, separated by '/'
#'
#' @examples
#' # Example data
#' stats <- rnorm(1000)
#' names(stats) <- paste0("Gene", 1:1000)
#' stats <- sort(stats, decreasing = TRUE)
#' 
#' gs1 <- paste0("Gene", 1:50)
#' gs2 <- paste0("Gene", 500:550)
#' gene_sets <- list(Pathway1 = gs1, Pathway2 = gs2)
#' 
#' # Use default fixed permutation method
#' result <- gsea(geneList=stats, gene_sets=gene_sets, nPerm=100)
#' 
#' # Use adaptive permutation for more accurate p-values
#' \donttest{
#' result_adaptive <- gsea(geneList=stats, gene_sets=gene_sets, adaptive=TRUE)
#' }
#'
#' @export
gsea <- function(geneList, gene_sets, 
                 weight = NULL,
                 minGSSize = 10,
                 maxGSSize = 500,
                 nPerm = 1000, 
                 exponent = 1.0, 
                 method = "multilevel",
                 adaptive = FALSE, 
                 minPerm = 101, 
                 maxPerm = 100000, 
                 pvalThreshold = 0.1, 
                 eps = 1e-10,
                 sampleSize = 101,
                 seed = FALSE,
                 nPermSimple = 1000,
                 scoreType = "std",
                 verbose = TRUE) {
    
    gene_sets <- validate_gene_sets(gene_sets)
    
    if (!is.null(weight)) {
        if (!is.numeric(weight) || is.null(names(weight))) {
            stop("weight must be a named numeric vector")
        }
        
        common_genes <- intersect(names(geneList), names(weight))
        if (length(common_genes) == 0) {
            stop("No common genes found between geneList and weight")
        }
        if (verbose && length(common_genes) < length(geneList)) {
            message("Weight provided: ", length(common_genes), " out of ", length(geneList), " genes are weighted.")
        }
        
        geneList <- geneList[common_genes]
        weight_aligned <- weight[common_genes]
        
        if (any(weight_aligned < 0, na.rm = TRUE)) {
            warning("Negative weights detected. Taking absolute values.")
            weight_aligned <- abs(weight_aligned)
        }
        
        # Apply weight
        geneList <- geneList * weight_aligned
        # Re-sort to maintain preranked requirement
        geneList <- sort(geneList, decreasing = TRUE)
    }
    
    method <- match.arg(method, c("sample", "permute", "multilevel"))
    scoreType <- match.arg(scoreType, c("std", "pos", "neg"))

    prepared <- prepare_gsea_inputs(geneList, scoreType, exponent)
    geneList <- prepared$geneList
    multilevelRanks <- prepared$multilevelRanks
    sampleSize <- normalize_multilevel_sample_size(sampleSize)
    # seed semantics: a numeric (or TRUE -> fixed default seed) makes the result
    # reproducible; FALSE draws a fresh seed from R's RNG on each run, so results
    # may differ between runs (inherent to permutation-based GSEA). Note that the
    # C++ engine seeds its own RNG with this value, so passing an explicit seed
    # (or set.seed() before the call when seed = FALSE) fixes the result.
    if (isTRUE(seed)) {
        seed <- 12345L  # fixed default seed, consistent with the C++ default
    } else if (isFALSE(seed)) {
        seed <- sample.int(1e9, 1)
    }
    
    # Intersect with geneList to get overlapping genes;
    # minGSSize/maxGSSize should constrain the overlap, not the raw set size
    # see https://github.com/YuLab-SMU/clusterProfiler/issues/824
    gene_sets <- lapply(gene_sets, intersect, names(geneList))

    # Filter by size
    idx <- get_geneSet_index(gene_sets, minGSSize, maxGSSize)
    if (sum(idx) == 0) {
        if (verbose) {
            msg <- paste("No gene sets have size between", minGSSize, "and", maxGSSize, "...")
            message(msg)
            message("--> return NULL...")
        }
        return (NULL)
    }
    gene_sets <- gene_sets[idx]   

    gene_set_names <- names(gene_sets)
    
    # Call appropriate C++ function
    if (method == "multilevel") {
        result <- gsea_multilevel_cpp(geneList = multilevelRanks, gene_sets = gene_sets, gene_set_names = gene_set_names,
                                 minPerm = minPerm, maxPerm = maxPerm, pvalThreshold = pvalThreshold,
                                 exponent = 1.0, method = method, eps = eps, sampleSize = sampleSize, seed = seed,
                                 nPermSimple = nPermSimple, scoreType = scoreType)
    } else if (adaptive) {
        result <- gsea_adaptive_cpp(geneList, gene_sets, gene_set_names, 
                                    minPerm, maxPerm, pvalThreshold, exponent, method, seed)
    } else {
        result <- gsea_cpp(geneList, gene_sets, gene_set_names, nPerm, exponent, method, seed)
    }
    
    # Rename columns to standard names
    if (!"ID" %in% names(result) && "GeneSet" %in% names(result)) names(result)[names(result) == "GeneSet"] <- "ID"
    if (!"enrichmentScore" %in% names(result) && "ES" %in% names(result)) names(result)[names(result) == "ES"] <- "enrichmentScore"
    if (!"pvalue" %in% names(result) && "PValue" %in% names(result)) names(result)[names(result) == "PValue"] <- "pvalue"
    if (!"setSize" %in% names(result) && "Size" %in% names(result)) names(result)[names(result) == "Size"] <- "setSize"

    # Add setSize if missing
    if (!"setSize" %in% names(result)) {
        # Calculate size based on gene_sets (filtered)
        # Note: gene_sets order must match result ID order
        # Assuming result ID matches gene_set_names order passed to C++
        # If C++ preserves order (it does), we can map
        set_sizes <- sapply(gene_sets, length)
        result$setSize <- set_sizes[result$ID]
    }

    # Add NES if missing (placeholder)
    if (!"NES" %in% names(result)) {
        result$NES <- NA
    }

    if (!all(c("rank", "leading_edge", "core_enrichment") %in% names(result))) {
        ledge_rank <- integer(nrow(result))
        ledge_str <- character(nrow(result))
        ledge_core <- character(nrow(result))
        for (i in seq_len(nrow(result))) {
            sid <- result$ID[[i]]
            gs <- gene_sets[[sid]]
            ledge <- gsea_leading_edge_details(geneList, gs, exponent = exponent, scoreType = scoreType)
            ledge_rank[[i]] <- ledge$rank
            ledge_str[[i]] <- ledge$leading_edge
            ledge_core[[i]] <- ledge$core_enrichment
        }

        if (!"rank" %in% names(result)) result$rank <- ledge_rank
        if (!"leading_edge" %in% names(result)) result$leading_edge <- ledge_str
        if (!"core_enrichment" %in% names(result)) result$core_enrichment <- ledge_core
    }

    # Sort by absolute NES (descending) or pvalue (ascending)
    if (nrow(result) > 0) {
        if (!all(is.na(result$NES))) {
            result <- result[order(abs(result$NES), decreasing = TRUE), ]
        } else {
             result <- result[order(result$pvalue, decreasing = FALSE), ]
        }
    }

    if (method == "multilevel" && any(is.na(result$pvalue))) {
        warning(
            "There were ", sum(is.na(result$pvalue)),
            " pathways for which P-values were not calculated properly due to unbalanced gene-level statistic values. ",
            "For such pathways pvalue, NES and log2err are set to NA. You can try to increase nPermSimple."
        )
    }
    if (method == "multilevel" && any(!is.na(result$pvalue) & result$pvalue == eps & is.na(result$log2err))) {
        warning(
            "For some pathways, in reality P-values are less than ", eps,
            ". You can set the eps argument to zero for better estimation."
        )
    }
    
    return(result)
}

normalize_multilevel_sample_size <- function(sampleSize) {
    sampleSize <- as.integer(sampleSize[[1]])
    if (!is.finite(sampleSize)) {
        stop("sampleSize must be a finite integer")
    }
    sampleSize <- max(sampleSize, 3L)
    if (sampleSize %% 2L == 0L) {
        sampleSize <- sampleSize + 1L
    }
    sampleSize
}

prepare_gsea_inputs <- function(geneList, scoreType, exponent) {
    if (!is.numeric(geneList) || is.null(names(geneList))) {
        stop("geneList must be a named numeric vector")
    }
    if (any(!is.finite(geneList))) {
        stop("Not all stats values are finite numbers")
    }
    if (is.unsorted(rev(geneList))) {
        warning("geneList is not sorted in descending order. Sorting it now.")
        geneList <- sort(geneList, decreasing = TRUE)
    }

    ties <- sum(duplicated(geneList[geneList != 0]))
    if (ties != 0) {
        warning(
            "There are ties in the preranked stats (",
            round(ties * 100 / length(geneList), digits = 2),
            "% of the list). The order of those tied genes will be arbitrary, which may produce unexpected results."
        )
    }
    if (all(geneList > 0) && scoreType == "std") {
        warning(
            "All values in the stats vector are greater than zero and scoreType is \"std\", maybe you should switch to scoreType = \"pos\"."
        )
    }

    multilevelRanks <- scale_fgsea_ranks(geneList, exponent)

    list(
        geneList = geneList,
        multilevelRanks = multilevelRanks
    )
}

scale_fgsea_ranks <- function(geneList, exponent) {
    ranks <- abs(geneList)^exponent

    # Mirror fgsea::prepareStats(): scale to integer-like values based on the
    # total weight rather than a fixed multiplier, which keeps the multilevel
    # backend numerically closer to the long-used fgsea implementation.
    scaleCoeff <- 2^30 / sum(ranks)
    if (scaleCoeff >= 1) {
        scaleCoeff <- floor(scaleCoeff)
    }

    ranks <- round(ranks * scaleCoeff)
    storage.mode(ranks) <- "integer"
    structure(ranks, names = names(geneList))
}


#' generic function for gene set enrichment analysis
#'
#'
#' @title gsea_gson
#' @inheritParams enrichit_params
#' @param weight A named numeric vector of weights for genes.
#' @param seed Random seed for reproducibility. See \code{gsea()} for details
#'   (default: FALSE).
#' @param pvalueCutoff P-value cutoff applied to both the raw p-value and the
#'   adjusted p-value (\code{p.adjust}), consistent with the historical
#'   clusterProfiler/DOSE behavior (default: 0.05).
#' @param ... Additional parameters passed to gsea()
#' @return gseaResult object
#' @author Guangchuang Yu
#' @export
gsea_gson <- function(geneList,
                 gson,
                 weight = NULL,
                 nPerm = 1000,
                 exponent = 1.0,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 method = "multilevel",
                 adaptive = FALSE,
                 minPerm = 101,
                 maxPerm = 100000,
                 pvalThreshold = 0.1,
                 seed = FALSE,
                 verbose = TRUE,
                 ...) {

    if (!inherits(gson, "GSON")) {
        stop("gson should be a GSON object")
    }

    # Ensure geneList is sorted
    if (is.unsorted(rev(geneList))) {
        if (verbose) {
            warning("geneList is not sorted in descending order. Sorting it now.")
        }
        geneList <- sort(geneList, decreasing = TRUE)
    }

    ## query external ID to Term ID
    gene <- names(geneList)
    
    # Extract gene sets from GSON
    gsid2gene <- gson@gsid2gene
    
    # ID Match Check
    if (!check_gene_id(gene, gsid2gene)) {
        return(NULL)
    }

    # Prepare Gene Sets
    geneSets <- split(gsid2gene$gene, gsid2gene$gsid)
    
    gsea_res <- gsea(geneList = geneList, 
                     gene_sets = geneSets, 
                     weight = weight,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     nPerm = nPerm, 
                     exponent = exponent, 
                     method = method,
                     adaptive = adaptive,
                     minPerm = minPerm,
                     maxPerm = maxPerm,
                     pvalThreshold = pvalThreshold,
                     seed = seed,
                     verbose = verbose,
                     ...)
                     
    if (is.null(gsea_res) || nrow(gsea_res) == 0) {
        return(NULL)
    }

    # Add Description
    gsid2name <- gson@gsid2name
    if (!is.null(gsid2name) && "ID" %in% names(gsea_res)) {
        description <- gsid2name$name[match(gsea_res$ID, gsid2name$gsid)]
        na_idx <- is.na(description)
        description[na_idx] <- gsea_res$ID[na_idx]
        gsea_res$Description <- description
    } else {
        if (!"Description" %in% names(gsea_res)) {
            gsea_res$Description <- gsea_res$ID
        }
    }

    # Calculate p.adjust
    gsea_res$p.adjust <- p.adjust(gsea_res$pvalue, method=pAdjustMethod)
    
    # Calculate qvalue
    gsea_res$qvalue <- calculate_qvalue(gsea_res$pvalue)
    
    # Filter by pvalueCutoff: both the raw p-value and the adjusted p-value
    # must pass the cutoff, matching the historical DOSE/clusterProfiler
    # behavior (see DOSE::GSEA_fgsea: res[res$pvalue <= pvalueCutoff, ] followed
    # by res[res$p.adjust <= pvalueCutoff, ]).
    if (!is.null(pvalueCutoff)) {
        gsea_res <- gsea_res[!is.na(gsea_res$pvalue), ]
        gsea_res <- gsea_res[gsea_res$pvalue <= pvalueCutoff, ]
        gsea_res <- gsea_res[gsea_res$p.adjust <= pvalueCutoff, ]
    }
    
    if (nrow(gsea_res) == 0) {
        return(NULL)
    }
    
    # Reorder columns
    expected_cols <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")
    
    # Ensure all expected columns exist to avoid error
    missing_cols <- setdiff(expected_cols, names(gsea_res))
    if (length(missing_cols) > 0) {
        warning("Missing columns in GSEA result: ", paste(missing_cols, collapse=", "), 
                ". Filling with NA. (Note: 'qvalue' may be missing if p-value distribution does not allow q-value calculation)")
        for (col in missing_cols) {
            gsea_res[[col]] <- NA
        }
    }

    if ("qvalues" %in% names(gsea_res) && !"qvalue" %in% names(gsea_res)) {
        gsea_res$qvalue <- gsea_res$qvalues
    }

    other_cols <- setdiff(names(gsea_res), expected_cols)
    gsea_res <- gsea_res[, c(expected_cols, other_cols)]
    
    # Set row names
    # Ensure IDs are robust (handle NA and duplicates)
    ids <- as.character(gsea_res$ID)
    
    # Handle NA IDs (convert to string "NA")
    if (any(is.na(ids))) {
        warning("NA values detected in gene set IDs. Replacing with string 'NA'.")
        ids[is.na(ids)] <- "NA"
    }
    
    # Handle Duplicate IDs
    if (any(duplicated(ids))) {
        dup_ids <- unique(ids[duplicated(ids)])
        warning("Duplicate gene set IDs detected: ", paste(head(dup_ids), collapse=", "), 
                "... (Total ", length(dup_ids), "). Unique suffixes added.")
        ids <- make.unique(ids)
    }
    
    rownames(gsea_res) <- ids
    
    params <- list(pvalueCutoff = pvalueCutoff,
                   nPerm = nPerm,
                   pAdjustMethod = pAdjustMethod,
                   exponent = exponent,
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize)
                   
    res <- new("gseaResult",
               result = gsea_res,
               organism = if (!is.null(gson@species)) gson@species else "UNKNOWN",
               setType = if (!is.null(gson@gsname)) gsub(".*;", "", gson@gsname) else "UNKNOWN",
               geneSets = geneSets,
               geneList = geneList,
               keytype = if (!is.null(gson@keytype)) gson@keytype else "UNKNOWN",
               permScores = matrix(), 
               params = params,
               gene2Symbol = character(), 
               readable = FALSE
              )
              
    return(res)
}


