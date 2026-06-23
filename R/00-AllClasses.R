#' Class "compareClusterResult"
#' This class represents the comparison result of gene clusters by GO
#' categories at specific level or GO enrichment analysis.
#'
#'
#' @name compareClusterResult-class
#' @aliases compareClusterResult-class show,compareClusterResult-method summary,compareClusterResult-method plot,compareClusterResult-method
#' @docType class
#' @slot compareClusterResult cluster comparing result
#' @slot geneClusters a list of genes
#' @slot fun one of groupGO, enrichGO and enrichKEGG
#' @slot gene2Symbol gene ID to Symbol
#' @slot keytype Gene ID type
#' @slot readable logical flag of gene ID in symbol or not.
#' @slot .call function call
#' @slot termsim Similarity between term
#' @slot method method of calculating the similarity between nodes
#' @slot dr dimension reduction result
#' @slot organism organism
#' @exportClass compareClusterResult
#' @author Guangchuang Yu \url{https://yulab-smu.top}
#' @exportClass compareClusterResult
#' @seealso 
#'   \code{\linkS4class{enrichResult}}
#' @keywords classes
setClass("compareClusterResult",
         representation = representation(
             compareClusterResult = "data.frame",
             geneClusters = "list",
             fun = "character",
             gene2Symbol    = "character",
             keytype        = "character",
             readable       = "logical",
             .call          = "call",
             termsim        = "matrix",
             method         = "character",
             dr             = "list",
             organism       = "character"
         )
         )

#' Class "enrichResult"
#' This class represents the result of enrichment analysis.
#'
#'
#' @name enrichResult-class
#' @aliases enrichResult-class show,enrichResult-method summary,enrichResult-method
#'
#' @docType class
#' @slot result enrichment analysis
#' @slot pvalueCutoff pvalueCutoff
#' @slot pAdjustMethod pvalue adjust method
#' @slot qvalueCutoff qvalueCutoff
#' @slot organism only "human" supported
#' @slot ontology biological ontology
#' @slot gene Gene IDs
#' @slot keytype Gene ID type
#' @slot universe background gene
#' @slot gene2Symbol mapping gene to Symbol
#' @slot geneSets gene sets
#' @slot readable logical flag of gene ID in symbol or not.
#' @slot termsim Similarity between term
#' @slot method method of calculating the similarity between nodes
#' @slot dr dimension reduction result
#' @exportClass enrichResult
#' @author Guangchuang Yu \url{https://yulab-smu.top}
#' @keywords classes
setClass("enrichResult",
         representation=representation(
             result         = "data.frame",
             pvalueCutoff   = "numeric",
             pAdjustMethod  = "character",
             qvalueCutoff   = "numeric",
             organism       = "character",
             ontology       = "character",
             gene           = "character",
             keytype        = "character",
             universe       = "character",
             gene2Symbol    = "character",
             geneSets       = "list",
             readable       = "logical",
             termsim        = "matrix",
             method         = "character",
             dr             = "list"
             ),
         prototype=prototype(readable = FALSE)
         )



#' Class "gseaResult"
#' This class represents the result of GSEA analysis
#'
#'
#' @name gseaResult-class
#' @aliases gseahResult-class show,gseaResult-method summary,gseaResult-method
#'
#' @docType class
#' @slot result GSEA anaysis
#' @slot organism organism
#' @slot setType setType
#' @slot geneSets geneSets
#' @slot geneList order rank geneList
#' @slot keytype ID type of gene
#' @slot permScores permutation scores
#' @slot params parameters
#' @slot gene2Symbol gene ID to Symbol
#' @slot readable whether convert gene ID to symbol
#' @slot dr dimension reduction result
#' @exportClass gseaResult
#' @author Guangchuang Yu \url{https://yulab-smu.top}
#' @keywords classes
setClass("gseaResult",
         representation   = representation(
             result          = "data.frame",
             organism        = "character",
             setType         = "character",
             geneSets        = "list",
             geneList        = "numeric",
             keytype         = "character",
             permScores      = "matrix",
             params          = "list",
             gene2Symbol     = "character",
             readable        = "logical",
             termsim         = "matrix",
             method          = "character",
             dr              = "list"
         )
         )


#' Class "nseaResult"
#' This class represents the result of Network-based Set Enrichment Analysis (NSEA).
#'
#' @name nseaResult-class
#' @aliases nseaResult-class
#' @docType class
#' @slot result enrichment analysis
#' @slot organism organism label for the enrichment result
#' @slot setType gene set collection type
#' @slot geneSets gene sets
#' @slot geneList order rank geneList
#' @slot keytype ID type of gene
#' @slot permScores permutation score matrix inherited from `gseaResult`
#' @slot gene2Symbol gene ID to symbol mapping
#' @slot readable logical flag of gene ID in symbol or not.
#' @slot termsim Calculation matrix of termsim.
#' @slot method Method of termsim.
#' @slot params parameters
#' @slot dr dimension reduction result
#' @slot network sparse matrix or data.frame representing the underlying network.
#' @slot diffusion_scores numeric vector of RWR diffusion scores for each node.
#' @slot mode character, "evidence" or "signed", describing the RWR propagation mode.
#' @slot iterations integer, the actual number of iterations RWR took to converge.
#' @slot restart_prob numeric, the restart probability used in RWR.
#' @exportClass nseaResult
#' @author Guangchuang Yu \url{https://yulab-smu.top}
setClass("nseaResult",
         contains = "gseaResult",
         representation = representation(
             network = "ANY",
             diffusion_scores = "numeric",
             mode = "character",
             iterations = "integer",
             restart_prob = "numeric"
         )
)

