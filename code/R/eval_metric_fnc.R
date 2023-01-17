## ============================================================================
## It contains five functions: `GGMdensity`, `DCGdensity`,`truemoddensity`,
## `GGMdegree`, and `DCGdegree`.
##
## Purpose: It computes density and degrees per node for
## true cyclic models, GGMs, and DCGs. ============================================================================

#' Compute density for GGM
#' @param qgraph.object the qgraph object
#'
#' @return density value of GGM
GGMdensity <- function(qgraph.object){
  # extract number of nodes
  p <- qgraph.object[["graphAttributes"]][["Graph"]][["nNodes"]]
  # compute overall density
  density <- sum(round(abs(qgraph.object$Edgelist$weight), 2) > 0) / ((p * (p-1))/2)
  return(density)
}



#' Compute average density for equivalent class of DCGs (directed cyclic graph)
#' @param equivclass the list of equivalence class of DCG matrices
#'
#' @return the average density of equivalent class of DCGs
DCGdensity <- function (equivclass){
  # storage for density
  avg_density <- c()
  for(i in 1:length(equivclass)){
    # extract number of nodes
    p <- ncol(equivclass[[i]])
    # compute the number of edges such that
    # if there is a cycle between a pair, count it as one.
    x <- equivclass[[i]] + t(equivclass[[i]])
    x <- ifelse(x == 2, 1, x)
    edgenumber <- sum(x[lower.tri(x, diag = FALSE)])
    # compute the overall density for each of DCGs
    avg_density[i] <- edgenumber/ ((p * (p-1))/2)
  }
  # average all the obtained density
  return(mean(avg_density))
}



#' Compute density of true model (directed cyclic graph)
#' @param B regression matrix of true model
#'
#' @return the density of true DCG
truemoddensity <- function (B){
  # get the number of edges
  p <- ncol(B)
  # convert it to an adjency matrix
  x <- ifelse(B==0, 0, 1)
  # compute the number of edges such that
  # if there is a cycle between a pair, count it as one.
  x <- x + t(x)
  x <- ifelse(x == 2, 1, x)
  edgenumber <- sum(x[lower.tri(x, diag = FALSE)])
  # compute the density
  density <- edgenumber/ ((p * (p-1))/2)
  return(density)
}



#' Compute degree centrality for GGM
#' @param qgraph.object the qgraph object
#'
#' @return the dataframe of degrees per node
GGMdegree <- function(qgraph.object){
  # get the existing edge indices
  ind <- round(abs(qgraph.object$Edgelist$weight), 2) >0
  # how many edges per node
  tb <- table(c(qgraph.object$Edgelist$from[ind], qgraph.object$Edgelist$to[ind]))
  # extract the node labels for table
  names(tb) <- qgraph.object$graphAttributes$Nodes$names[as.numeric(names(tb))]
  degree <- as.data.frame(tb)
  colnames(degree) <- c("node", "degree")
  return(degree)
}



#' Compute average degree for equivalent class of DCGs (directed cyclic graph)
#' @param equivclass the list of equivalence class of DCG matrices
#'
#' @return the dataframe of average degree per node in equivalent class of DCGs
DCGdegree <- function (equivclass){
  # storage for overall degree per node
  overall_degrees <- list()
  for(i in 1:length(equivclass)){
    # compute the overall degree
    outdegree <- colSums(equivclass[[i]])
    indegree <- rowSums(equivclass[[i]])
    overalldegree <- outdegree + indegree
    overall_degrees[[i]] <- overalldegree
  }
  avg_deg <- do.call(rbind, overall_degrees)
  # compute the average degree
  avg_degree <- as.data.frame(apply(avg_deg, 2, mean))
  avg_degree$node <- rownames(avg_degree)
  colnames(avg_degree) <- c("average_degree", "node")
  rownames(avg_degree) <- NULL
  return(avg_degree[,c("node","average_degree")])
}

