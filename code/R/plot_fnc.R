## ============================================================================
## It contains four functions:
## `CreateAdjMat`, `extracttriples`, `extractnode`, and `plotPAG`.
##
## Purpose: It plots a partial ancestral graph (PAG).
## `plotPAG` is specifically designed for the ccd object.
## `plotAG` plots a PAG given an (ancestral) adjacency matrix.
## ============================================================================

## load necessary package
library(Rgraphviz)



#' Create adjacency matrix
#'
#' @param ccd.obj the resulting object of ccdKP function
#' @param p the number of nodes
#'
#' @details
#' "0": no edge; "1": circle; "2": arrow; "3": tail
#'
#' @return an adjacency matrix of PAG
CreateAdjMat <- function(ccd.obj, p){
  # extract the edge list
  ccd.edges <- ccd.obj$edges
  # storage matrix
  mat <- matrix(0, p, p, dimnames = list(ccd.obj$nodes, ccd.obj$nodes))

  # for now, I assume the names of nodes are numbers.
  for (i in 1:length(ccd.edges)){
    nodes <- unlist(stringr::str_split(ccd.edges[i], " "))
    row <- nodes[1]
    column <- nodes[3]
    # extract the edge type
    row_edgetype <- substr(nodes[2], 1, 1)
    column_edgetype <- substr(nodes[2], 3, 3)
    # fill the matrix
    mat[row, column] <- column_edgetype
    mat[column, row] <- row_edgetype
  }
  # encode each type with the corresponding number
  mat[mat=="o"] <- 1
  mat[mat==">"|mat=="<"] <- 2
  mat[mat=="-"] <- 3
  # convert the character matrix to numeric matrix
  class(mat) <- "numeric"

  return(mat)
}



#' Extract triples
#'
#' @param triples the triples in PAG from a ccd object
#'
#' @return the first, middle, last node in each triple
extracttriples <- function(triples) {
  # extract all triples from ccd object
  triplesets <- unlist(stringr::str_extract_all(triples, "(?<=\\<).+?(?=\\>)"))
  # separate each triple
  triplesets <- stringr::str_split(triplesets, ",")
  # extract first node
  firstnode <- stringr::str_trim(sapply(triplesets, function(x) x[[1]]))
  # extract second node
  middlenode <- stringr::str_trim(sapply(triplesets, function(x) x[[2]]))
  # extract third node
  lastnode <- stringr::str_trim(sapply(triplesets, function(x) x[[3]]))
  return(data.frame(fistnode = firstnode, middlenode= middlenode, lastnode=lastnode))
}



#' Extract middle node in a triple
#'
#' @param triples the triples in PAG
#'
#' @return the middle node in the triples
extractnode <- function(triples) {
  # extract all triples from ccd object
  triplesets <- unlist(stringr::str_extract_all(triples, "(?<=\\<).+?(?=\\>)"))
  # separate each triple
  triplesets <- stringr::str_split(triplesets, ",")
  # extract only the middle node
  middlenode <- stringr::str_trim(sapply(triplesets, function(x) x[[2]]))
  return(middlenode)
}



#' #' Extract dotted-underlined triplet nodes
#' #'
#' #' @param triples the triples in PAG
#' #'
#' #' @return the first,second,third node in the dotted-underlined triples
#' #'
#' extractdottednodes <- function(dottedtriples) {
#'   triplesets <- unlist(stringr::str_extract_all(dottedtriples,  "(?<=\\<).+?(?=\\>)"))
#'   triplesets <- stringr::str_split(triplesets, ",")
#'   firstnode <- stringr::str_trim(sapply(triplesets, function(x) x[[1]]))
#'   middlenode <- stringr::str_trim(sapply(triplesets, function(x) x[[2]]))
#'   lastnode <- stringr::str_trim(sapply(triplesets, function(x) x[[3]]))
#'   return(data.frame(fistnode = firstnode, middlenode= middlenode, lastnode=lastnode))
#' }



#' Plot PAG (partial ancestral graph) for CCD algorithm
#'
#' @param amat the adjacency matrix of PAG
#' @param ccd.obj the resulting object of ccdKP function
#'
#' @details
#' "0": no edge; "1": circle; "2": arrow; "3": tail
#'
#' @return a PAG graph of graphNEL class
plotPAG <- function(ccd.obj, amat)
{
  # get the underline triples
  underline <- ccd.obj$graph$getUnderLines()$toString()
  # get the dotted-underline triples
  dottedunderline <- ccd.obj$graph$getDottedUnderlines()$toString()
  # get underline node
  underline_node <- extractnode(underline)
  # get  dotted-underline node
  dottedunderline_node <- extractnode(dottedunderline)
  # convert it to a graphNEL object
  g <- as(amat, "graphNEL")
  # extract node info
  nn <- nodes(g)
  p <- numNodes(g)
  # extract edge info
  n.edges <- numEdges(g)
  ah.list <- at.list <- vector("list", n.edges)
  l.names <- character(n.edges)
  # assign a shape for each edge type
  amat[amat == 1] <- "odot"
  amat[amat == 2] <- "vee"
  amat[amat == 3] <- "none"
  iE <- 0
  for (i in 1:(p-1)) {
    x <- nn[i]
    for (j in (i+1):p) {
      y <- nn[j]
      if (amat[x,y] != 0) {
        iE <- iE + 1
        ah.list[[iE]] <- amat[x,y]
        at.list[[iE]] <- amat[y,x]
        l.names[[iE]] <- paste0(x,"~",y)
      }
    }
  }
  names(ah.list) <- names(at.list) <- l.names
  # create a graph object
  graph <- layoutGraph(g) # layoutType="neato"
  # for the solid underlines, color them in blue
  fill <- rep("lightblue", length(underline_node))
  names(fill) <- underline_node
  nodeRenderInfo(graph) <- list(fill = fill)
  # for the dotted underlines, make the outlines dashed
  lty <- rep(2, length(dottedunderline_node))
  names(lty) <- dottedunderline_node
  nodeRenderInfo(graph) <- list(lty = lty)

  edgeRenderInfo(graph) <- list(arrowhead = ah.list, arrowtail = at.list)
  # global features
  graph.par(list(nodes=list(cex = 1)))
  # plot the pag
  pag <- Rgraphviz::renderGraph(graph)
  # also show the underlining information
  cat(paste("Dotted-underlined triples:", dottedunderline, "\nUnderlined triples:", underline, "\n"))
  return(pag)
}



## =======================================
## plotAG function (edited slightly from `pcalg` package)
## =======================================
#' Plot PAG (partial ancestral graph) for FCI and CCI algorithms
#'
#' @param amat adjacency matrix of the resulting estimated graph
#'
#' @details
#' "0": no edge; "1": circle; "2": arrow; "3": tail
#'
#' @return a PAG graph of graphNEL class
plotAG <- function(amat)
{
  g <- as(amat,"graphNEL")
  nn <- nodes(g)
  p <- numNodes(g)
  n.edges <- numEdges(g)
  ah.list <- at.list <- vector("list", n.edges)
  l.names <- character(n.edges)
  amat[amat == 1] <- "odot"
  amat[amat == 2] <- "vee"
  amat[amat == 3] <- "none"
  iE <- 0
  for (i in 1:(p-1)) {
    x <- nn[i]
    for (j in (i+1):p) {
      y <- nn[j]
      if (amat[x,y] != 0) {
        iE <- iE + 1
        ah.list[[iE]] <- amat[x,y]
        at.list[[iE]] <- amat[y,x]
        l.names[[iE]] <- paste0(x,"~",y)
      }
    }
  }
  names(ah.list) <- names(at.list) <- l.names
  
  edgeRenderInfo(g) <- list(arrowhead = ah.list, arrowtail = at.list)
  # global features
  graph.par(list(nodes=list(cex = 1)))
  Rgraphviz::renderGraph(Rgraphviz::layoutGraph(g))
}

