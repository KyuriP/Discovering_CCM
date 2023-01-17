## ============================================================================
## It contains a function `dsep_KP`.
##
## Purpose: It tests if the set "a" and the set "b" are d-separated given the set "S".
## This function is used when finding all equivalence class of DCGs
## in `semiequiv_dcg` function from `equivset_fnc.R`.
## ============================================================================

dsep_KP <- function(a,b, S = NULL, graph)
{
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - a,b,S: vectors of node names
  ## - graph: pag or cdgset
  ## ----------------------------------------------------------------------
  ## Value:
  ## Boolean decision
  ## ----------------------------------------------------------------------

  amatTmp <- graph ## i->j if amatTmp[j,i]!=0
  amatTmp[amatTmp != 0] <- 1
  ## create a graphNEL object
  g <- as(t(amatTmp), "graphNEL")
  p <- numNodes(g)

  ## build node union of a,b,S
  nodeUnion <- if(length(S) > 0) c(a,b,S) else c(a,b)
  my.nodes <- nodes(g) # check the nodes

  ## find ancestor graph of node union
  anc.set <- NULL
  for (i in seq_len(p)) {
    desc.nodes <- my.nodes[possDe(m = amatTmp, x = i, possible = TRUE, ds = FALSE, type = "pdag") ]
    if (any(desc.nodes %in% nodeUnion)) anc.set <- c(anc.set, my.nodes[i])
  } ## for (i in 1:p)
  gS <- subGraph(anc.set, g)

  ## Moralize in amatM
  ## !!! in the following line:
  ## i->j if amat[i,j], i.e. different than default coding !!!
  ## (*)
  amat <- wgtMatrix(gS, transpose = FALSE)
  if(all(a0 <- amat == 0))
    ## if no edge in graph, nodes are d-separated
    return( TRUE )
  ## else :
  amat[!a0] <- 1
  amatM <- amat
  ind <- which(amat == 1,arr.ind = TRUE)

  for (i in seq_len(nrow(ind))) {
    ## input is guaranteed to be directed
    x <- ind[i,1]
    y <- ind[i,2] ## x -> y
    ## using different coding, see (*) -> OK
    allZ <- setdiff(which(amat[y,] == 0 & amat[,y] == 1), x) ## x -> y <- z
    for (z in allZ)
      if (amat[x,z] == 0 && amat[z,x] == 0)
        amatM[x,z] <- 1 ## moralize
  } ## for (i in seq_len(nrow(ind)))

  ## make undirected graph
  ## (up to now, there is NO undirected edge -> just add t(amat))
  gSM <- as(amatM+t(amatM),"graphNEL")

  if (length(S) > 0) { ## check separation
    RBGL::separates(a,b,S,gSM)
  } else {
    b %nin% (if(is.list(bfs. <- bfs(gSM,a))) bfs.[[1]] else bfs.)
  }
} ## {dsep}
