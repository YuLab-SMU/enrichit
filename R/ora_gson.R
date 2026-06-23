#' interal method for enrichment analysis
#'
#' using the hypergeometric model
#' @title ora-gson
#' @param gene a vector of entrez gene id.
#' @param qvalueCutoff cutoff of qvalue
#' @param universe background genes, default is the intersection of the 'universe' with genes that have annotations. 
#' Users can set `options(enrichment_force_universe = TRUE)` to force the 'universe' untouched.
#' @param weight A named numeric vector of weights for background genes. If provided, Weighted ORA will be performed.
#' @inheritParams enrichit_params
#' @return  A `enrichResult` instance.
#' @importClassesFrom methods data.frame
#' @importFrom methods new
#' @importFrom stats p.adjust
#' @keywords manip
#' @author Guangchuang Yu <https://yulab-smu.top>
#' @export
ora_gson <- function(gene,
                              pvalueCutoff,
                              pAdjustMethod="BH",
                              universe = NULL,
                              weight = NULL,
                              minGSSize=10,
                              maxGSSize=500,
                              qvalueCutoff=0.2,
                              gson){

    if (!inherits(gson, "GSON")) {
        stop("gson should be a GSON object")
    }

    ## query external ID to Term ID
    gene <- as.character(unique(gene))
    
    # Extract gene sets from GSON
    gsid2gene <- gson@gsid2gene
    
    # ID Match Check
    if (!check_gene_id(gene, gsid2gene)) {
        return(NULL)
    }

    # Handle universe
    extID <- unique(gsid2gene$gene)
    if (missing(universe))
        universe <- NULL
    if(!is.null(universe)) {
        if (is.character(universe)) {
            force_universe <- getOption("enrichment_force_universe", FALSE)
            if (force_universe) {
                extID <- universe
            } else {
                extID <- intersect(extID, universe)
            }
        } else {
            ## https://github.com/YuLab-SMU/clusterProfiler/issues/217
            message("`universe` is not in character and will be ignored...")
        }
    }

    # Prepare Gene Sets
    geneSets <- split(gsid2gene$gene, gsid2gene$gsid)
    # Intersect with universe
    geneSets <- lapply(geneSets, intersect, extID)
    
    # Filter by size
    idx <- get_geneSet_index(geneSets, minGSSize, maxGSSize)
    if (sum(idx) == 0) {
        msg <- paste("No gene sets have size between", minGSSize, "and", maxGSSize, "...")
        message(msg)
        message("--> return NULL...")
        return (NULL)
    }
    geneSets <- geneSets[idx]

    ora_res <- ora(gene, geneSets, universe = extID, weight = weight)

    if (is.null(ora_res) || nrow(ora_res) == 0) {
        return(NULL)
    }

    # Rename columns to match expectation
    if ("PValue" %in% names(ora_res)) {
        names(ora_res)[names(ora_res) == "PValue"] <- "pvalue"
    }
    if ("GeneSet" %in% names(ora_res)) {
        names(ora_res)[names(ora_res) == "GeneSet"] <- "ID"
    }
    if ("DEInSet" %in% names(ora_res)) {
        names(ora_res)[names(ora_res) == "DEInSet"] <- "Count"
    }

    # Only gene sets with at least one matched input gene should contribute to
    # multiple testing correction. Keeping zero-overlap rows here inflates the
    # number of hypotheses and shifts adjusted statistics relative to the
    # historical DOSE/clusterProfiler behavior.
    if ("Count" %in% names(ora_res)) {
        ora_res <- ora_res[ora_res$Count > 0, , drop = FALSE]
    }

    if (nrow(ora_res) == 0) {
        return(NULL)
    }
    
    # Calculate ratios
    if (all(c("Count", "DESize") %in% names(ora_res))) {
        ora_res$GeneRatio <- paste0(ora_res$Count, "/", ora_res$DESize)
    }
    if (all(c("SetSize", "UniverseSize") %in% names(ora_res))) {
        ora_res$BgRatio <- paste0(ora_res$SetSize, "/", ora_res$UniverseSize)
    }

    # Calculate RichFactor and FoldEnrichment
    if (all(c("Count", "SetSize", "DESize", "UniverseSize") %in% names(ora_res))) {
        ora_res$RichFactor <- ora_res$Count / ora_res$SetSize
        ora_res$FoldEnrichment <- (ora_res$Count / ora_res$DESize) / (ora_res$SetSize / ora_res$UniverseSize)
    }

    p <- suppressWarnings(as.numeric(ora_res$pvalue))
    invalid <- is.na(p) | !is.finite(p) | p < 0 | p > 1
    if (any(invalid)) {
        ids <- as.character(ora_res$ID)
        bad_i <- which(invalid)
        head_i <- bad_i[seq_len(min(length(bad_i), 10))]
        bad_preview <- paste0(ids[head_i], "=", format(p[head_i], digits = 4, scientific = TRUE), collapse = ", ")
        warning(
            "Invalid p-values detected in ORA result (",
            length(bad_i),
            "/",
            length(p),
            "). Examples: ",
            bad_preview
        )
        too_small <- !is.na(p) & is.finite(p) & p < 0
        too_large <- !is.na(p) & is.finite(p) & p > 1
        p[too_small] <- 0
        p[too_large] <- 1
        ora_res$pvalue <- p
    } else if (!isTRUE(all.equal(p, ora_res$pvalue))) {
        ora_res$pvalue <- p
    }

    # Calculate p.adjust
    ora_res$p.adjust <- p.adjust(ora_res$pvalue, method=pAdjustMethod)

    # Calculate qvalue
    ora_res$qvalue <- calculate_qvalue(ora_res$pvalue)

    # Calculate zScore
    # Need k, M, n, N    
    N <- length(extID)
    n <- length(intersect(gene, extID))
    k <- ora_res$Count
    
    # M = size of gene set (in universe)
    # Map ID to M
    # ora_res$ID should match names(geneSets)
    M <- sapply(geneSets[ora_res$ID], length)
    
    mu <- M * n / N
    sigma <- mu * (N - n) * (N - M) / N / (N - 1)
    zScore <- (k - mu) / sqrt(sigma)
    ora_res$zScore <- zScore

    # Add Description
    gsid2name <- gson@gsid2name
    if (!is.null(gsid2name) && "ID" %in% names(ora_res)) {
        description <- gsid2name$name[match(ora_res$ID, gsid2name$gsid)]
        na_idx <- is.na(description)
        description[na_idx] <- ora_res$ID[na_idx]
        ora_res$Description <- description
    } 

    if (!"Description" %in% names(ora_res)) {
        ora_res$Description <- ora_res$ID
    }

    # Reorder columns
    expected_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    
    # Ensure all expected columns exist
    missing_cols <- setdiff(expected_cols, names(ora_res))
    if (length(missing_cols) > 0) {
        warning("Missing columns in ORA result: ", paste(missing_cols, collapse=", "), 
                ". Filling with NA. (Note: 'qvalue' may be missing if p-value distribution does not allow q-value calculation)")
        for (col in missing_cols) {
            ora_res[[col]] <- NA
        }
    }

    #other_cols <- setdiff(names(ora_res), expected_cols)
    ora_res <- ora_res[, expected_cols]
    
    # Sort by pvalue
    ora_res <- ora_res[order(ora_res$pvalue), ]
    
    # Set row names
    # Ensure IDs are robust (handle NA and duplicates)
    ids <- as.character(ora_res$ID)
    
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
    
    rownames(ora_res) <- ids

    x <- new("enrichResult",
             result         = ora_res,
             pvalueCutoff   = pvalueCutoff,
             pAdjustMethod  = pAdjustMethod,
             qvalueCutoff   = qvalueCutoff,
             gene           = as.character(gene),
             universe       = extID,
             geneSets       = geneSets,
             organism       = if (!is.null(gson@species)) gson@species else "UNKNOWN",
             keytype        = if (!is.null(gson@keytype)) gson@keytype else "UNKNOWN",
             ontology       = if (!is.null(gson@gsname)) gsub(".*;", "", gson@gsname) else "UNKNOWN",
             readable       = FALSE
             )
             
    return(x)
}


