## ============================================================================
## It contains two functions `DCGdensities` and `DCGdegrees`.
##
## Purpose: To compute the density and degree for an equivalence class of DCGs
## per model. It is used in `variation_DCG.R` where the corresponding figures 
## are created.
## ============================================================================



#' Density for equivalence class of DCGs
#'
#' @description `DCGdensities` compute density for equivalence class of DCGs (directed cyclic graph)
#' @param equivclass the list of equivalence class of DCG matrices
#'
#' @return the list containing the size of equivalence class, density of all equivalent DCGs, average density of DCGs, and variation in the set of densities
DCGdensities <- function (equivclass){
  # storage for densities
  density <- c()
  # for all DCGs in the equivalence class
  for(i in 1:length(equivclass)){
    # p: number of nodes
    p <- ncol(equivclass[[i]])
    # only count the edge once when there is a cycle between a pair
    x <- equivclass[[i]] + t(equivclass[[i]])
    x <- ifelse(x == 2, 1, x)
    # get the count for the number of edges
    edgenumber <- sum(x[lower.tri(x, diag = FALSE)])
    # compute the overall density
    density[i] <- edgenumber/ ((p * (p-1))/2)
  }
  return(list(class_size =length(equivclass), densities = density, avg_density = mean(density), variation = var(density)))
}



#' Degree for equivalence class of DCGs
#'
#' @description `DCGdegrees` compute degrees per node for equivalence class of DCGs (directed cyclic graph)
#' @param equivclass the list of equivalence class of DCG matrices
#'
#' @return the dataframe containing the list of node names, average degree per node, and standard deviation of degrees per node
DCGdegrees <- function (equivclass){
  # storage for the degrees per node
  overall_degrees <- list()
  # for all DCGs in the equivalence class
  for(i in 1:length(equivclass)){
    # compute out-degree per node
    outdegree <- colSums(equivclass[[i]])
    # compute in-degree per node
    indegree <- rowSums(equivclass[[i]])
    # sum out-degree & in-degree to get the overall degree
    overalldegree <- outdegree + indegree
    overall_degrees[[i]] <- overalldegree
  }
  # combine the overall degree for all nodes
  overall_deg <- do.call(rbind, overall_degrees)
  # compute the standard deviation in degrees per node
  deg_var <- as.data.frame(apply(overall_deg, 2, sd))
  colnames(deg_var) <- "sd"
  # compute the average degree per node
  avg_degree <- as.data.frame(apply(overall_deg, 2, mean))
  avg_degree$node <- rownames(avg_degree)
  colnames(avg_degree) <- c("average_degree", "node")
  rownames(avg_degree) <- NULL
  return(data.frame(node = avg_degree$node, avg_degree= avg_degree$average_degree, deg_sd = deg_var$sd))
}
