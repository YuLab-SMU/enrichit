#' interal method for enrichment analysis
#'
#' using the hypergeometric model
#' @title ora-gson
#' @param gene a vector of entrez gene id.
#' @param qvalueCutoff cutoff of qvalue
#' @param universe background genes, default is the intersection of the 'universe' with genes that have annotations. 
#' Users can set `options(enrichment_force_universe = TRUE)` to force the 'universe' untouched.
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

    ora_res <- ora(gene, geneSets, universe = extID)

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
    # other_cols <- setdiff(names(ora_res), expected_cols)
    # ora_res <- ora_res[, c(expected_cols, other_cols)]
    ora_res <- ora_res[, expected_cols]

    # Sort by pvalue
    ora_res <- ora_res[order(ora_res$pvalue), ]
    
    # Set row names
    row.names(ora_res) <- ora_res$ID

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


