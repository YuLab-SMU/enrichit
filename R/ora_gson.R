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
#' @return A `enrichResult` instance. The `result` slot carries the three
#'   significance columns `pvalue`, `p.adjust`, and `qvalue`:
#'   raw nominal p-values, BH-style multiple-testing adjusted p-values
#'   (via `stats::p.adjust()`; default `pAdjustMethod = "BH"`), and
#'   Storey q-values (via `qvalue::qvalue()`). If q-value estimation
#'   fails, `qvalue` remains `NA`. See `?as_enrichResult` for the full
#'   column semantics.
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
    annotatedGenes <- unique(gsid2gene$gene)
    extID <- annotatedGenes
    if (missing(universe))
        universe <- NULL
    if(!is.null(universe)) {
        ## A numeric/integer/factor universe (e.g. plain Entrez IDs supplied as
        ## numbers) is a perfectly meaningful background and should not be
        ## silently discarded -- coerce it to character the same way `gene`
        ## is coerced a few lines above. Only genuinely non-vector input
        ## (list, data.frame, S4 object, ...) is still rejected.
        ## see https://github.com/YuLab-SMU/clusterProfiler/issues/217
        ##     https://github.com/YuLab-SMU/clusterProfiler/issues/737
        if (!is.character(universe)) {
            if (is.atomic(universe) && !is.null(universe)) {
                message(
                    "`universe` is not in character, converting via as.character()..."
                )
                universe <- as.character(universe)
            } else {
                message("`universe` is not in character and will be ignored...")
                universe <- NULL
            }
        }

        if (!is.null(universe)) {
            force_universe <- getOption("enrichment_force_universe", FALSE)
            if (force_universe) {
                extID <- universe
            } else {
                extID <- intersect(extID, universe)
            }
        }
    }

    # Prepare Gene Sets
    geneSets <- split(gsid2gene$gene, gsid2gene$gsid)
    ## Intersect the gene sets with the background only when the background was
    ## actually RESTRICTED. An intersection with the full annotated gene list is
    ## the identity -- the sets come from that same table -- but it is not free:
    ## base `intersect()` hashes `extID` once per set, measured at 20 s for the
    ## 14,209 human GO BP sets (PBMC 10k, 2026-10-18), and ORA runs once per
    ## cluster in compareCluster(). Note that clusterProfiler passes the full
    ## annotated list EXPLICITLY when the caller supplies no universe, so
    ## `universe != NULL` alone is not the right test; `extID` is shorter than
    ## the annotated set only when a universe actually removed genes.
    if (!is.null(universe) && length(extID) < length(annotatedGenes)) {
        geneSets <- filter_gene_sets(geneSets, extID)
    }
    
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

    ## Fisher's exact odds ratio: what fisher.test()$estimate reports for the 2x2
    ## table below. FoldEnrichment is a ratio of proportions, whereas the odds
    ## ratio is the effect size the test itself is built on, so it is worth
    ## exposing next to it.
    ##
    ##            in gene set   not in gene set
    ##  in list        k            n - k
    ##  not in list  M - k     N - M - n + k
    ##
    ## A zero cell gives 0 or Inf, exactly as fisher.test() does.
    ora_res$oddsRatio <- (k * (N - M - n + k)) / ((n - k) * (M - k))

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
    expected_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment", "oddsRatio", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count")

    # Ensure all expected columns exist
    missing_cols <- setdiff(expected_cols, names(ora_res))
    if (length(missing_cols) > 0) {
        warning("Missing columns in ORA result: ", paste(missing_cols, collapse=", "),
                ". Filling with NA. (Note: 'qvalue' may be missing if p-value distribution does not allow q-value calculation)")
        for (col in missing_cols) {
            ora_res[[col]] <- NA
        }
    }

    ## Keep the denominators the engine already computed. `DESize` is the
    ## GeneRatio denominator (query genes inside the universe) and
    ## `UniverseSize` the BgRatio denominator (universe ∩ annotated genes);
    ## both were calculated here and then dropped, which is why four users in a
    ## row (clusterProfiler#703, #523, #618, #723) had to ask what the
    ## denominators of the "20/207" strings actually are. `SetSize` stays too:
    ## RichFactor divides by it, so all three ratios are now checkable from the
    ## table itself. GSEA's gsea_gson() has always kept its extra columns.
    extra_cols <- intersect(c("SetSize", "DESize", "UniverseSize"), names(ora_res))
    ora_res <- ora_res[, c(expected_cols, extra_cols)]
    
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
