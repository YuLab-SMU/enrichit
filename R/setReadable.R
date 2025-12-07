#' mapping geneID to gene Symbol
#'
#'
#' @title setReadable
#' @param x enrichResult Object
#' @param OrgDb OrgDb
#' @param keyType keyType of gene
#' @return enrichResult Object
#' @author Guangchuang Yu
# @importFrom yulab.utils load_OrgDb
#' @export
setReadable <- function(x, OrgDb, keyType="auto") {
    OrgDb <- load_OrgDb(OrgDb)
    if (!'SYMBOL' %in% AnnotationDbi::columns(OrgDb)) {
        warning("Fail to convert input geneID to SYMBOL since no SYMBOL information available in the provided OrgDb...")
    }

    if (!(is(x, "enrichResult") || is(x, "groupGOResult") || is(x, "gseaResult") || is(x,"compareClusterResult")))
        stop("input should be an 'enrichResult' , 'gseaResult' or 'compareClusterResult' object...")

    isGSEA <- FALSE
    isCompare <- FALSE
    if (is(x, 'gseaResult'))
        isGSEA <- TRUE

    if (is(x, 'compareClusterResult'))
        isCompare <- TRUE

    if (keyType == "auto") {
        keyType <- x@keytype
        if (keyType == 'UNKNOWN') {
            stop("can't determine keyType automatically; need to set 'keyType' explicitly...")
        }
    }

    if (x@readable)
        return(x)

    gc <- geneInCategory(x)
    if (isGSEA) {
        genes <- names(x@geneList)
    } else if (isCompare) {
        if ("core_enrichment" %in% colnames(as.data.frame(x))) {
            geneslist <- x@geneClusters
            names(geneslist) <- NULL
            genes <- unique(names(unlist(geneslist)))
        } else {
            genes <- unique(unlist(x@geneClusters))
        }     
    } else {
        genes <- x@gene
    }

    gn <- EXTID2NAME(OrgDb, genes, keyType)


    if(isCompare) {
        gc2 <- list()
        k <- 1
        for(i in seq_len(length(gc))) {
            for(j in seq_len(length(gc[[i]]))) {
                gc2[[k]] <- gc[[i]][[j]]
                names(gc2)[k] <- paste(names(gc)[[i]], names(gc[[i]])[j], sep="-")
                k <- k + 1
            }
        }
        gc <- gc2
        gc <- lapply(gc, function(i) gn[i])
        res <- x@compareClusterResult
        gc <- gc[paste(res$Cluster, res$ID, sep= "-")]
    } else {
        gc <- lapply(gc, function(i) gn[i])
        res <- x@result
        gc <- gc[as.character(res$ID)]
    }

    ## names(gc) should be identical to res$ID

    ## gc <- gc[as.character(res$ID)]


    geneID <- sapply(gc, paste0, collapse="/")
    # if (isGSEA) {
    if ("core_enrichment" %in% colnames(as.data.frame(x))) {
        res$core_enrichment <- unlist(geneID)
    } else {
        res$geneID <- unlist(geneID)
    }
    x@gene2Symbol <- gn
    x@keytype <- keyType
    x@readable <- TRUE
    if(isCompare){
        x@compareClusterResult <- res
    } else {
        x@result <- res
    }


    return(x)
}

#' mapping gene ID to gene Symbol
#'
#'
#' @title EXTID2NAME
#' @param OrgDb OrgDb
#' @param geneID entrez gene ID
#' @param keytype keytype
#' @return gene symbol
# @importFrom yulab.utils load_OrgDb
#' @export
#' @author Guangchuang Yu \url{https://yulab-smu.top}
EXTID2NAME <- function(OrgDb, geneID, keytype) {
    OrgDb <- load_OrgDb(OrgDb)
    kt <- AnnotationDbi::keytypes(OrgDb)
    if (! keytype %in% kt) {
        stop("keytype is not supported...")
    }

    gn.df <- suppressMessages(
        AnnotationDbi::select(OrgDb, keys=geneID, keytype=keytype, columns="SYMBOL")
    )

    gn.df <- unique(gn.df)
    colnames(gn.df) <- c("GeneID", "SYMBOL")

    unmap_geneID <- geneID[!geneID %in% gn.df$GeneID]
    if (length(unmap_geneID) != 0) {
        unmap_geneID.df = data.frame(GeneID = unmap_geneID,
                                     SYMBOL = unmap_geneID)
        gn.df <- rbind(gn.df, unmap_geneID.df)
    }

    gn <- gn.df$SYMBOL
    names(gn) <- gn.df$GeneID
    return(gn)
}


# to remove and imported from yulab.utils
load_OrgDb <- function(OrgDb) {
    #if (is(OrgDb, "character")) {
    #    require(OrgDb, character.only = TRUE)
    #    OrgDb <- eval(parse(text=OrgDb))
    #}
    if (is(OrgDb, "character")) {
        OrgDb <- utils::getFromNamespace(OrgDb, OrgDb)
    } 
    
    return(OrgDb)
}
