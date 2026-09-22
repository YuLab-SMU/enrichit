## Converters: build enrichResult / gseaResult objects from external
## enrichment-analysis result tables (enrichr, g:Profiler, WebGestalt,
## fgsea, ...).  Importers for specific tools live in 'enrichplot' and
## call these constructors.

##' Convert a result table to an \code{enrichResult} object
##'
##' Generic constructor for over-representation analysis results produced
##' by other tools.  The input must follow the canonical column schema
##' documented in \code{enrichplot::fortify()}; a minimal set of columns
##' (\code{ID}, \code{pvalue} and at least one of \code{geneID},
##' \code{Count} or \code{GeneRatio}) is required, the remaining canonical
##' columns are derived when possible.
##' @title as_enrichResult
##' @param x result table (data.frame) with the canonical ORA columns, or
##' a table with common aliases (e.g. \code{PValue}, \code{term_id});
##' tool-specific naming should be mapped to the canonical schema by the
##' per-tool importers in 'enrichplot'
##' @param ... additional arguments passed to methods
##' @return An \code{enrichResult} object
##' @export
as_enrichResult <- function(x, ...) {
    UseMethod("as_enrichResult")
}

##' @rdname as_enrichResult
##' @param geneSets gene sets as a named list, a two-column
##' data.frame (term, gene), or a \code{GSON} object. If \code{NULL}, gene
##' sets are rebuilt from the \code{geneID} column (overlap genes only,
##' sufficient for \code{cnetplot()}/\code{heatplot()}).
##' @param gene query gene vector used in the analysis. If \code{NULL},
##' inferred as the union of \code{geneID} entries.
##' @param universe background gene vector. If \code{NULL}, inferred as the
##' union of \code{geneSets} when available.
##' @param ontology,organism,keytype metadata stored in the object slots.
##' @param pAdjustMethod method passed to \code{stats::p.adjust} when
##' \code{p.adjust} is missing.
##' @export
as_enrichResult.default <- function(
    x,
    geneSets = NULL,
    gene = NULL,
    universe = NULL,
    ontology = "UNKNOWN",
    organism = "UNKNOWN",
    keytype = "UNKNOWN",
    pAdjustMethod = "BH",
    ...
) {
    df <- .prepare_result_df(x)
    df <- .rename_canonical_ora(df)

    df <- .finalize_ids(df)

    if (is.null(df$Description)) {
        df$Description <- df$ID
    }

    ## validate and normalize pvalue / p.adjust / qvalue early
    df <- .finalize_pvalues(df, pAdjustMethod)

    ## geneID: normalize separators
    if (!is.null(df$geneID)) {
        df$geneID <- gsub("[;,]\\s*", "/", trimws(as.character(df$geneID)))
        genes_list <- strsplit(df$geneID, "/", fixed = TRUE)
    } else {
        genes_list <- NULL
    }

    gene_provided <- !is.null(gene)
    if (is.null(gene)) {
        if (!is.null(genes_list)) {
            gene <- unique(unlist(genes_list))
        } else if (!is.null(geneSets)) {
            gene <- unique(unlist(.normalize_geneSets(geneSets)))
        }
    }
    gene <- as.character(gene)

    ## Count
    if (is.null(df$Count)) {
        if (!is.null(genes_list)) {
            df$Count <- lengths(genes_list)
        } else if (!is.null(df$GeneRatio)) {
            df$Count <- as.numeric(sub("/.*$", "", as.character(df$GeneRatio)))
        } else {
            stop("at least one of 'geneID', 'Count' or 'GeneRatio' is required")
        }
    }
    df$Count <- as.numeric(df$Count)

    ## geneSets: keep user-provided sets separate from the rebuilt slot value,
    ## so derived statistics never use overlap-only rebuilt sets
    geneSets_norm <- .normalize_geneSets(geneSets)

    user_ratio <- df$GeneRatio

    ## GeneRatio
    if (is.null(df$GeneRatio)) {
        if (length(gene) > 0) {
            df$GeneRatio <- paste0(df$Count, "/", length(gene))
            if (!gene_provided) {
                warning(
                    "query genes not provided; 'GeneRatio' uses the union of ",
                    "overlap genes as the denominator and may be inaccurate. ",
                    "Pass 'gene' for exact ratios."
                )
            }
        } else {
            df$GeneRatio <- NA_character_
        }
    } else if (is.numeric(df$GeneRatio) && gene_provided && length(gene) > 0) {
        ## rebuild "k/n" so downstream ratio parsing and cluster labeling work
        df$GeneRatio <- paste0(df$Count, "/", length(gene))
    }
    ## otherwise keep as-is: character "k/n" or numeric ratios both supported

    ## BgRatio (only from user-provided gene sets)
    if (is.null(df$BgRatio) && length(geneSets_norm) > 0) {
        if (is.null(universe)) {
            universe <- unique(unlist(geneSets_norm))
        }
        df$BgRatio <- vapply(
            df$ID,
            function(i) paste0(length(geneSets_norm[[i]]), "/", length(universe)),
            character(1)
        )
    }
    if (!is.null(df$BgRatio)) {
        df$BgRatio <- as.character(df$BgRatio)
    }
    if (is.null(universe) && length(geneSets_norm) > 0) {
        universe <- unique(unlist(geneSets_norm))
    }
    universe <- as.character(universe)

    ## optional derived statistics (only when all components are known)
    N_universe <- if (length(universe) > 0) length(universe) else NA_integer_
    M_set <- .set_sizes(df, geneSets_norm)

    if (!"RichFactor" %in% names(df) && !anyNA(M_set)) {
        df$RichFactor <- df$Count / M_set
    }
    n_query <- if (gene_provided && length(gene) > 0) {
        length(gene)
    } else if (!is.null(user_ratio)) {
        .n_from_ratio(user_ratio)
    } else {
        NULL
    }
    if (!is.na(N_universe) && !is.null(n_query) && !anyNA(M_set)) {
        k <- df$Count
        if (!"FoldEnrichment" %in% names(df)) {
            df$FoldEnrichment <- (k / n_query) / (M_set / N_universe)
        }
        if (!"zScore" %in% names(df)) {
            mu <- M_set * n_query / N_universe
            sigma <- mu * (N_universe - n_query) * (N_universe - M_set) /
                N_universe / (N_universe - 1)
            sigma[sigma <= 0] <- NA_real_
            df$zScore <- (k - mu) / sqrt(sigma)
        }
    }

    ## geneSets slot: rebuild from geneID when not provided
    if (length(geneSets_norm) == 0 && !is.null(genes_list)) {
        geneSets_norm <- split(unlist(genes_list), rep(df$ID, lengths(genes_list)))
        geneSets_norm <- lapply(geneSets_norm, unique)
    }

    df <- .order_columns_ora(df)
    rownames(df) <- df$ID

    methods::new(
        "enrichResult",
        result = df,
        pvalueCutoff = 1,
        pAdjustMethod = pAdjustMethod,
        qvalueCutoff = 1,
        organism = organism,
        ontology = ontology,
        gene = gene,
        keytype = keytype,
        universe = universe,
        gene2Symbol = character(),
        geneSets = geneSets_norm,
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list()
    )
}

##' Convert a result table to a \code{gseaResult} object
##'
##' Generic constructor for pre-ranked GSEA results produced by other tools
##' (e.g. 'fgsea', Broad GSEA reports).  The ranked gene statistics must be
##' supplied separately as \code{geneList}; the result table alone cannot be
##' converted because most GSEA visualizations consume the ranked list.
##' @title as_gseaResult
##' @param x result table (data.frame) with GSEA columns
##' (\code{ID}, \code{enrichmentScore}, \code{pvalue}; commonly also
##' \code{NES}, \code{p.adjust}, \code{setSize}, \code{core_enrichment};
##' fgsea-style \code{pathway}/\code{ES}/\code{pval}/\code{padj}/\code{size}
##' are recognized)
##' @param ... additional arguments passed to methods
##' @return A \code{gseaResult} object
##' @export
as_gseaResult <- function(x, ...) {
    UseMethod("as_gseaResult")
}

##' @rdname as_gseaResult
##' @param geneList named numeric vector of ranked statistics, sorted in
##' descending order (sorted automatically with a warning if not).
##' @param geneSets gene sets as a named list, a two-column data.frame
##' (term, gene), or a \code{GSON} object. If \code{NULL}, gene sets are
##' rebuilt from \code{core_enrichment} (leading edge genes only) and
##' running-score plots will only be approximate.
##' @param setType,organism,keytype metadata stored in the object slots.
##' @param exponent,scoreType parameters used to recompute missing
##' \code{rank}/\code{leading_edge}/\code{core_enrichment} columns.
##' @param pAdjustMethod method passed to \code{stats::p.adjust} when
##' \code{p.adjust} is missing.
##' @export
as_gseaResult.default <- function(
    x,
    geneList,
    geneSets = NULL,
    setType = "UNKNOWN",
    organism = "UNKNOWN",
    keytype = "UNKNOWN",
    exponent = 1,
    scoreType = "std",
    pAdjustMethod = "BH",
    ...
) {
    prep <- prepare_gsea_inputs(geneList, scoreType, exponent)
    geneList <- prep$geneList

    df <- .prepare_result_df(x)
    df <- .rename_canonical_gsea(df)
    df <- .finalize_ids(df)

    if (is.null(df$Description)) {
        df$Description <- df$ID
    }

    if (is.null(df$enrichmentScore)) {
        stop("the result table must contain 'enrichmentScore' (or 'ES')")
    }

    ## validate and normalize pvalue / p.adjust / qvalue early
    df <- .finalize_pvalues(df, pAdjustMethod)

    ## core_enrichment: prefer leadingEdge list column (fgsea) if present
    if (!is.null(df$.leadingEdge_list)) {
        le <- vapply(
            df$.leadingEdge_list,
            function(g) paste0(unique(as.character(g)), collapse = "/"),
            character(1)
        )
        if (is.null(df$core_enrichment)) {
            df$core_enrichment <- le
        }
        df$.leadingEdge_list <- NULL
    }
    if (!is.null(df$core_enrichment)) {
        df$core_enrichment <- gsub("[;,]\\s*", "/", trimws(as.character(df$core_enrichment)))
    }

    geneSets <- .normalize_geneSets(geneSets)
    if (length(geneSets) == 0) {
        if (is.null(df$core_enrichment)) {
            stop("either 'geneSets' or a 'core_enrichment' column is required")
        }
        warning(
            "geneSets is NULL: rebuilding gene sets from 'core_enrichment'; ",
            "running-score plots will only be approximate"
        )
        cores <- strsplit(as.character(df$core_enrichment), "/", fixed = TRUE)
        geneSets <- lapply(seq_len(nrow(df)), function(i) unique(cores[[i]]))
        names(geneSets) <- df$ID
    }
    geneSets <- lapply(geneSets, function(gs) as.character(unique(gs)))

    ## fill missing rank / leading_edge / core_enrichment
    missing_details <- is.null(df$rank) || is.null(df$leading_edge) ||
        is.null(df$core_enrichment)
    if (missing_details) {
        details <- lapply(seq_len(nrow(df)), function(i) {
            gsea_leading_edge_details(
                geneList,
                geneSets[[df$ID[i]]],
                exponent = exponent,
                scoreType = scoreType
            )
        })
        if (is.null(df$rank)) {
            df$rank <- vapply(details, `[[`, integer(1), "rank")
        }
        if (is.null(df$leading_edge)) {
            df$leading_edge <- vapply(details, `[[`, character(1), "leading_edge")
        }
        if (is.null(df$core_enrichment)) {
            df$core_enrichment <- vapply(details, `[[`, character(1), "core_enrichment")
        }
    }

    if (is.null(df$setSize)) {
        df$setSize <- vapply(seq_len(nrow(df)), function(i) {
            length(intersect(geneSets[[df$ID[i]]], names(geneList)))
        }, numeric(1))
    }

    if (!is.null(df$NES) && all(!is.na(df$NES))) {
        df <- df[order(abs(df$NES), decreasing = TRUE), , drop = FALSE]
    } else {
        df <- df[order(df$pvalue, decreasing = FALSE), , drop = FALSE]
    }
    rownames(df) <- df$ID
    df$rank <- as.integer(df$rank)

    methods::new(
        "gseaResult",
        result = df,
        organism = organism,
        setType = setType,
        geneSets = geneSets,
        geneList = geneList,
        keytype = keytype,
        permScores = matrix(0, 0, 0),
        params = list(exponent = exponent, scoreType = scoreType),
        gene2Symbol = character(),
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list()
    )
}

## ---- internal helpers -------------------------------------------------

.prepare_result_df <- function(x) {
    df <- tryCatch(as.data.frame(x), error = function(e) NULL)
    if (is.null(df)) {
        stop("x cannot be converted to a data.frame")
    }
    if (nrow(df) == 0) {
        stop("the result table is empty")
    }
    df
}

.rename_first <- function(df, canonical, aliases) {
    if (!is.null(df[[canonical]])) {
        return(df)
    }
    hit <- intersect(aliases, names(df))
    if (length(hit) == 0) {
        return(df)
    }
    names(df)[names(df) == hit[1]] <- canonical
    df
}

.rename_canonical_ora <- function(df) {
    df <- .rename_first(df, "ID", c("GeneSet", "term_id", "term", "pathway", "geneSet"))
    df <- .rename_first(df, "Description", c("description", "term_name", "name"))
    df <- .rename_first(df, "pvalue", c("PValue", "p_value", "P.value"))
    df <- .rename_first(
        df, "p.adjust",
        c("padj", "p_adj", "FDR", "Adjusted.P.value", "Benjamini")
    )
    df <- .rename_first(df, "qvalue", c("q_values", "q.value"))
    df <- .rename_first(df, "geneID", c("gene_id", "Genes", "overlapId", "intersection"))
    df <- .rename_first(df, "Count", c("count", "intersection_size", "overlap"))
    df <- .rename_first(df, "GeneRatio", c("geneRatio"))
    df <- .rename_first(df, "BgRatio", c("bgRatio"))
    df
}

.rename_canonical_gsea <- function(df) {
    df <- .rename_first(df, "ID", c("pathway", "GeneSet", "term_id", "term", "geneSet"))
    df <- .rename_first(df, "Description", c("description", "term_name", "name"))
    df <- .rename_first(df, "enrichmentScore", c("ES"))
    df <- .rename_first(df, "pvalue", c("pval", "PValue", "p_value", "P.value"))
    df <- .rename_first(df, "p.adjust", c("padj", "FDR", "Adjusted.P.value"))
    df <- .rename_first(df, "setSize", c("size", "Size", "set_size", "term_size"))
    df <- .rename_first(df, "core_enrichment", c("core_enrichment_genes"))
    ## keep fgsea list column aside; handled after geneSets normalization
    if (!is.null(df$leadingEdge) && is.list(df$leadingEdge)) {
        df$.leadingEdge_list <- df$leadingEdge
        df$leadingEdge <- NULL
    }
    df
}

.finalize_ids <- function(df) {
    if (is.null(df$ID)) {
        stop("the result table must contain an 'ID' column")
    }
    ids <- as.character(df$ID)
    if (anyNA(ids)) {
        warning("NA values detected in IDs. Replacing with string 'NA'.")
        ids[is.na(ids)] <- "NA"
    }
    if (any(duplicated(ids))) {
        dups <- unique(ids[duplicated(ids)])
        warning(
            length(dups), " duplicated ID(s) detected: ",
            paste(utils::head(dups, 5), collapse = ", "),
            ". Unique suffixes added."
        )
        ids <- make.unique(ids)
    }
    df$ID <- ids
    df
}

.finalize_pvalues <- function(df, pAdjustMethod) {
    if (is.null(df$pvalue)) {
        stop("the result table must contain a 'pvalue' column")
    }
    p <- suppressWarnings(as.numeric(df$pvalue))
    invalid <- !is.na(p) & (!is.finite(p) | p < 0 | p > 1)
    if (any(invalid)) {
        warning(
            sum(invalid), " invalid p-value(s) clamped to [0, 1]."
        )
        p[!is.na(p) & p < 0] <- 0
        p[!is.na(p) & p > 1] <- 1
        p[!is.finite(p) & !is.na(p)] <- 1
    }
    if (all(is.na(p))) {
        stop("no valid p-values in the result table")
    }
    df$pvalue <- p

    if (is.null(df$p.adjust)) {
        df$p.adjust <- stats::p.adjust(df$pvalue, method = pAdjustMethod)
    } else {
        df$p.adjust <- suppressWarnings(as.numeric(df$p.adjust))
    }

    if (is.null(df$qvalue)) {
        q <- calculate_qvalue(df$pvalue)
        ## keep the column non-empty when qvalue estimation fails
        q[is.na(q)] <- df$p.adjust[is.na(q)]
        df$qvalue <- q
    } else {
        df$qvalue <- suppressWarnings(as.numeric(df$qvalue))
    }
    df
}

.n_from_ratio <- function(ratio) {
    ratio <- as.character(ratio)
    if (!any(grepl("/", ratio))) {
        return(NULL)
    }
    n <- suppressWarnings(as.numeric(sub("^.*/", "", ratio)))
    if (all(is.na(n))) {
        return(NULL)
    }
    n
}

.set_sizes <- function(df, geneSets) {
    if (length(geneSets) > 0 && all(df$ID %in% names(geneSets))) {
        return(vapply(df$ID, function(i) length(geneSets[[i]]), numeric(1)))
    }
    m <- suppressWarnings(as.numeric(sub("/.*$", "", as.character(df$BgRatio))))
    if (all(is.na(m))) {
        return(rep(NA_real_, nrow(df)))
    }
    m
}

.normalize_geneSets <- function(geneSets) {
    if (is.null(geneSets)) {
        return(list())
    }
    if (inherits(geneSets, "GSON")) {
        geneSets <- split(
            as.character(geneSets@gsid2gene$gene),
            as.character(geneSets@gsid2gene$gsid)
        )
        return(geneSets)
    }
    if (is.data.frame(geneSets) || is.matrix(geneSets)) {
        geneSets <- as.data.frame(geneSets)
        if (ncol(geneSets) < 2) {
            stop("a geneSets data.frame needs two columns: term, gene")
        }
        geneSets <- split(
            as.character(geneSets[[2]]),
            as.character(geneSets[[1]])
        )
        return(geneSets)
    }
    if (!is.list(geneSets) || is.null(names(geneSets))) {
        stop(
            "geneSets must be a named list, a two-column data.frame ",
            "(term, gene), or a GSON object"
        )
    }
    lapply(geneSets, function(gs) as.character(unique(gs)))
}

.order_columns_ora <- function(df) {
    expected <- c(
        "ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
        "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue",
        "geneID", "Count"
    )
    keep <- c(intersect(expected, names(df)), setdiff(names(df), expected))
    df[, keep, drop = FALSE]
}
