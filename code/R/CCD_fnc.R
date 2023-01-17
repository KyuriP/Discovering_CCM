## ============================================================================
## It contains five function: `ccdKP`, `loadContinuousData`, `loadDiscreteData`,
## `extractTetradNodes`, `extractTetradEdges`.
## The main function is `ccdKP`, which runs the CCD algorithm,
## and it is dependent on the other four functions.
##
## Purpose: To run the CCD algorithm.
##
## Note: This is a wrapper function for the ccd algorithm in the Tetrad software.
## Please check Tetrad for more detailed information
## on its CCD function (https://sites.google.com/view/tetradcausal).
## ============================================================================

## ====================
## preparation
## ====================
# install.packages("stringr")
# install.packages("rJava")
# install.packages("devtools")
# install.packages("DOT")
# install_github("bd2kccd/r-causal")
library(rJava)
library(usethis)
library(devtools)
library(rcausal)
library(DOT)
library(pcalg)



ccdKP <- function (df, dataType = "continuous", numCategoriesToDiscretize = 4,
                  depth = 3, alpha = 0.05, numBootstrap = -1, ensembleMethod = "Highest",
                  java.parameters = NULL, priorKnowledge = NULL)
{
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - df: data
  ## - dataType: data type either continuous or discrete
  ## - numCategoriesToDiscretize: number of categories when wanting to discretize
  ## - depth: search depth (default = 3)
  ## - alpha: alpha value for conditional independence test
  ## - numBootstrap: number of bootstrapping performed (default = -1 : no bootstrap)
  ## - ensembleMethod: type of ensemble method to use
  ## - java.parameters: when specify tetrad java parameters
  ## - priorKnowledge: when incorporate additional knowledge

  ## ----------------------------------------------------------------------
  ## Value:
  ## ccd model object
  ## ----------------------------------------------------------------------
  params <- list(NULL)
  if (!is.null(java.parameters)) {
    options(java.parameters = java.parameters)
    params <- c(java.parameters = java.parameters)
  }
  indTest <- NULL
  if (dataType == "continuous") {
    tetradData <- loadContinuousData(df)
    if (numBootstrap < 1) {
      indTest <- .jnew("edu/cmu/tetrad/search/IndTestFisherZ",
                       tetradData, alpha)
    }
    else {
      indTest <- .jnew("edu/cmu/tetrad/algcomparison/independence/FisherZ")
    }
  }
  else if (dataType == "discrete") {
    tetradData <- loadDiscreteData(df)
    if (numBootstrap < 1) {
      indTest <- .jnew("edu/cmu/tetrad/search/IndTestChiSquare",
                       tetradData, alpha)
    }
    else {
      indTest <- .jnew("edu/cmu/tetrad/algcomparison/independence/ChiSquare")
    }
  }
  else {
    tetradData <- loadMixedData(df, numCategoriesToDiscretize)
    if (numBootstrap < 1) {
      indTest <- .jnew("edu/cmu/tetrad/search/IndTestConditionalGaussianLRT",
                       tetradData, alpha)
    }
    else {
      indTest <- .jnew("edu/cmu/tetrad/algcomparison/independence/ConditionalGaussianLRT")
    }
  }
  ccd <- list()
  class(ccd) <- "ccd"
  ccd$datasets <- deparse(substitute(df))
  cat("Datasets:\n")
  cat(deparse(substitute(df)), "\n\n")
  ccd_instance <- NULL
  if (numBootstrap < 1) {
    indTest <- .jcast(indTest, "edu/cmu/tetrad/search/IndependenceTest")
    ccd_instance <- .jnew("edu/cmu/tetrad/search/Ccd", indTest)
    .jcall(ccd_instance, "V", "setDepth", as.integer(depth))
  }
  else {
    indTest <- .jcast(indTest, "edu/cmu/tetrad/algcomparison/independence/IndependenceWrapper")
    algorithm <- .jnew("edu/cmu/tetrad/algcomparison/algorithm/oracle/pag/Ccd",
                       indTest)
    algorithm <- .jcast(algorithm, "edu/cmu/tetrad/algcomparison/algorithm/Algorithm")
    parameters_instance <- .jnew("edu/cmu/tetrad/util/Parameters")
    obj_depth <- .jnew("java/lang/Integer", as.integer(depth))
    parameter_instance <- .jcast(obj_depth, "java/lang/Object")
    parameters_instance$set("depth", parameter_instance)
    obj_alpha <- .jnew("java/lang/Double", alpha)
    parameter_instance <- .jcast(obj_alpha, "java/lang/Object")
    parameters_instance$set("alpha", parameter_instance)
    #obj_verbose <- .jnew("java/lang/Boolean", verbose)
    #parameter_instance <- .jcast(obj_verbose, "java/lang/Object")
    #parameters_instance$set("verbose", parameter_instance)
    ccd_instance <- .jnew("edu/pitt/dbmi/algo/bootstrap/GeneralBootstrapTest",
                          tetradData, algorithm, numBootstrap)
    edgeEnsemble <- .jfield("edu/pitt/dbmi/algo/bootstrap/BootstrapEdgeEnsemble",
                            name = ensembleMethod)
    ccd_instance$setEdgeEnsemble(edgeEnsemble)
    ccd_instance$setParameters(parameters_instance)
  }
  #ccd_instance$setVerbose(verbose)
  if (!is.null(priorKnowledge)) {
    .jcall(ccd_instance, "V", "setKnowledge", priorKnowledge)
  }
  params <- c(params, dataType = dataType) # why as.integer(dataType)??
  params <- c(params, depth = as.integer(depth))
  params <- c(params, alpha = alpha)
  if (numBootstrap > 0) {
    params <- c(params, numBootstrap = as.integer(numBootstrap))
    params <- c(params, ensembleMethod = ensembleMethod)
  }
  params <- c(params)
  if (!is.null(priorKnowledge)) {
    params <- c(params, prior = priorKnowledge)
  }
  ccd$parameters <- params
  cat("Graph Parameters:\n")
  cat("dataType = ", dataType, "\n")
  cat("depth = ", as.integer(depth), "\n")
  cat("alpha = ", as.numeric(alpha), "\n")
  if (numBootstrap > 0) {
    cat("numBootstrap = ", as.integer(numBootstrap), "\n")
    cat("ensembleMethod = ", ensembleMethod, "\n")
  }
  #cat("verbose = ", verbose, "\n")
  tetrad_graph <- .jcall(ccd_instance, "Ledu/cmu/tetrad/graph/Graph;",
                         "search", check = FALSE)
  if (!is.null(e <- .jgetEx())) {
    .jclear()
    ccd$nodes <- NULL
    ccd$edges <- NULL
    print("Java exception was raised")
    print(e)
  }
  else {
    ccd$graph <- tetrad_graph
    V <- extractTetradNodes(tetrad_graph)
    ccd$nodes <- V
    ccd_edges <- extractTetradEdges(tetrad_graph)
    ccd$edges <- ccd_edges
  }
  return(ccd)
}


############################################################
loadContinuousData <- function(df){
  node_names <- colnames(df)
  node_list <- .jnew("java/util/ArrayList")
  for (i in 1:length(node_names)){
    nodname <- .jnew("java/lang/String", node_names[i])
    nodi <- .jnew("edu/cmu/tetrad/data/ContinuousVariable", nodname)
    node_list$add(nodi)
  }
  node_list <- .jcast(node_list, "java/util/List")
  mt <- as.matrix(df)
  mat <- .jarray(mt, dispatch=TRUE)

  data <- .jnew("edu/cmu/tetrad/data/DoubleDataBox", mat)
  data <- .jcast(data, "edu/cmu/tetrad/data/DataBox")
  boxData <- .jnew("edu/cmu/tetrad/data/BoxDataSet", data, node_list)
  boxData <- .jcast(boxData, "edu/cmu/tetrad/data/DataSet")
  return(boxData)
}

############################################################
loadDiscreteData <- function(df){
  node_names <- colnames(df)
  node_list <- .jnew("java/util/ArrayList")
  for (i in 1:length(node_names)){
    nodname <- .jnew("java/lang/String", node_names[i])
    cat("node_names: ", node_names[i],"\n")
    cate <- unique(df[[node_names[i]]])
    cate <- sort(cate)
    cat("value: ")
    print(cate)
    cat("\n")
    cate_list <- .jnew("java/util/ArrayList")
    for(j in 1:length(cate)){
      cate_list$add(as.character(cate[j]))
    }
    cate_list <- .jcast(cate_list, "java/util/List")
    nodi <- .jnew("edu/cmu/tetrad/data/DiscreteVariable",
                  nodname, cate_list)
    node_list$add(nodi)

    # Substitute a new categorial value
    cate <- data.frame(cate)
    new_col <- sapply(df[,i],function(x,cate)
      as.integer(which(cate[,1] == x)),cate=cate)
    new_col = as.integer(new_col - 1)
    df[,i] <- (data.frame(new_col))[,1]
  }
  node_list <- .jcast(node_list, "java/util/List")
  mt <- as.matrix(df)
  mat <- .jarray(t(mt), dispatch=TRUE)
  data <- .jnew("edu/cmu/tetrad/data/VerticalIntDataBox", mat)
  data <- .jcast(data, "edu/cmu/tetrad/data/DataBox")
  boxData <- .jnew("edu/cmu/tetrad/data/BoxDataSet", data, node_list)
  boxData <- .jcast(boxData, "edu/cmu/tetrad/data/DataSet")
  return(boxData)
}

############################################################
extractTetradNodes <- function(resultGraph){
  nods <- resultGraph$getNodes()
  V <- sapply(as.list(nods), with, toString())
  return(V)
}

############################################################
extractTetradEdges <- function(resultGraph){
  eds <- resultGraph$getEdges()
  fgs_edges <- c()
  if(!is.null(eds)){
    fgs_edges <- sapply(as.list(eds), .jrcall, "toString")
  }
  return(fgs_edges)
}
