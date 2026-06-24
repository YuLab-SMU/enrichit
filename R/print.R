#' show method for `gseaResult` instance
#'
#' @name show
#' @docType methods
#' @rdname show-methods
#'
#' @title show method
#' @return message
#' @importFrom methods show
#' @exportMethod show
#' @usage show(object)
#' @author Guangchuang Yu <https://yulab-smu.top>
setMethod("show", signature(object="gseaResult"),
          function (object){
              params <- object@params
              cat("#\n# Gene Set Enrichment Analysis\n#\n")
              .print_common_info(object)
              
              cat("#...@geneList", "\t")
              str(object@geneList)
              cat("#...nPerm", "\t", params$nPerm, "\n")
              cat(sprintf("#...pvalues adjusted by '%s' with cutoff < %s\n", 
                          params$pAdjustMethod, params$pvalueCutoff))
              cat(sprintf("#...%d enriched terms found\n", nrow(object@result)))
              str(object@result)
              cat("#...Citation\n")
              print_citation_msg(object@setType)
          }
)

#' show method for `nseaResult` instance
#'
#' @name show
#' @aliases show,nseaResult-method
#' @docType methods
#' @rdname show-methods
#'
#' @title show method
#' @param object A `nseaResult` instance.
#' @return message
#' @importFrom utils str
#' @importFrom methods show
#' @exportMethod show
#' @usage show(object)
#' @author Guangchuang Yu <https://yulab-smu.top>
setMethod("show", signature(object="nseaResult"),
          function (object){
              params <- object@params
              cat("#\n# Network-based Set Enrichment Analysis (NSEA)\n#\n")
              .print_common_info(object)
              
              cat("#...@mode", "\t", object@mode, "\n")
              cat(sprintf("#...@network\t%d nodes, %d edges\n", 
                          nrow(object@network), sum(object@network > 0)))
              cat(sprintf("#...@rwr\t%d iterations to converge (restart_prob = %s)\n", 
                          object@iterations, object@restart_prob))
              
              cat("#...@geneList (input)", "\t")
              str(object@geneList)
              cat("#...@diffusion_scores", "\t")
              str(object@diffusion_scores)
              
              cat(sprintf("#...pvalues adjusted by '%s' with cutoff < %s\n", 
                          params$pAdjustMethod, params$pvalueCutoff))
              cat(sprintf("#...%d enriched terms found\n", nrow(object@result)))
              str(object@result)
              cat("#...Citation\n")
              print_citation_msg(object@setType)
          }
)

#' show method for `mnseaResult` instance
#'
#' @name show
#' @aliases show,mnseaResult-method
#' @docType methods
#' @rdname show-methods
#'
#' @title show method
#' @param object A `mnseaResult` instance.
#' @return message
#' @importFrom utils str
#' @importFrom methods show
#' @exportMethod show
#' @usage show(object)
#' @author Guangchuang Yu <https://yulab-smu.top>
setMethod("show", signature(object="mnseaResult"),
          function (object){
              params <- object@params
              cat("#\n# Multi-layer Network-based Set Enrichment Analysis (MNSEA)\n#\n")
              .print_common_info(object)
              
              cat("#...@mode", "\t", object@mode, "\n")
              cat("#...@collapse", "\t", object@collapse_method, "\n")
              cat("#...@output_space", "\t", object@output_space, "\n")
              if (nzchar(object@target_layer)) {
                  cat("#...@target_layer", "\t", object@target_layer, "\n")
              }
              cat(sprintf("#...@layers\t%d layers, %d coupling edges\n",
                          length(object@layer_scores), nrow(object@coupling_table)))
              cat(sprintf("#...@rwr\t%d iterations to converge (restart_prob = %s)\n",
                          object@iterations, object@restart_prob))
              cat(sprintf("#...@explain\t%d pathway rows, %d feature rows cached\n",
                          nrow(object@pathway_contribution), nrow(object@feature_contribution)))
              
              cat("#...@collapsed_scores", "\t")
              str(object@collapsed_scores)
              cat(sprintf("#...pvalues adjusted by '%s' with cutoff < %s\n",
                          params$pAdjustMethod, params$pvalueCutoff))
              cat(sprintf("#...%d enriched terms found\n", nrow(object@result)))
              str(object@result)
              cat("#...Citation\n")
              print_citation_msg(object@setType)
          }
)


#' show method for `enrichResult` instance
#'
#' @name show
#' @docType methods
#' @rdname show-methods
#'
#' @title show method
#' @param object A `enrichResult` instance.
#' @return message
#' @importFrom utils str
#' @importFrom methods show
#' @exportMethod show
#' @usage show(object)
#' @author Guangchuang Yu <https://yulab-smu.top>
setMethod("show", signature(object="enrichResult"),
        function (object){
              
              cat("#\n# over-representation test\n#\n")
              .print_common_info(object)
              
              cat("#...@gene", "\t")
              str(object@gene)
              cat(sprintf("#...pvalues adjusted by '%s' with cutoff < %s\n", 
                          object@pAdjustMethod, object@pvalueCutoff))
              
              object <- get_enriched(object)
              n <- nrow(object@result)
              cat(sprintf("#...%d enriched terms found\n", n))
              if (n > 0) str(object@result)
              cat("#...Citation\n")
              print_citation_msg(object@ontology)
        }
)

#' @importFrom methods .hasSlot
.print_common_info <- function(object) {
    cat("#...@organism", "\t", object@organism, "\n")
    
    # Handle ontology or setType depending on object type
    if (.hasSlot(object, "ontology")) {
        cat("#...@ontology", "\t", object@ontology, "\n")
    } else if (.hasSlot(object, "setType")) {
        cat("#...@setType", "\t", object@setType, "\n")
    }
    
    kt <- object@keytype
    if (kt != "UNKNOWN") {
        cat("#...@keytype", "\t", kt, "\n")
    }
}


print_citation_msg <- function(ontology) {
    refs <- yulab.utils:::ref_knownledge()

    if (ontology == "HDO" || ontology == "NCG") {
        citation_msg <- refs["DOSE"]
    } else if (ontology == "Reactome") {
        citation_msg <- refs["ReactomePA"]
    } else if (ontology == "MeSH") {
        citation_msg <- refs["meshes"]
    } else {
        citation_msg <- refs["clusterProfiler_NP"]
    }
    cat(citation_msg, "\n\n")
}

