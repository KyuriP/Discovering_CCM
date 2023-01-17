## ============================================================================
## It contains a function `semiequiv_dcg`.
##
## Purpose: It computes all directed cyclic graph (DCGs) in an equivalence class
## based on a PAG estimated by the CCD algorithm.
##
## (Note: I've tested this function and it has performed well so far.
## But I am still not 100% sure if I've covered all possible cases.
## I am currently working on making a more thorough set of tests).
## ============================================================================


#' Search for equivalent class of DCGs
#' @param ccdobj the resulting object of ccdKP function
#' @param pag the resulting object of CreateAdjMat function
#'
#' @return list of equivalent class of DCGs
semiequiv_dcg <- function(ccdobj, pag, verbose = FALSE)
  {
  # possible dcg matrix starting with fully-connected:
  dcgset <- ifelse(pag==0, 0, 1)
  num_node <- ncol(dcgset)
  ## STEP 1) first removing graph that doesn't satisfy ancestral relationships
  equivset <- list() # storage for the set
  # size of target vector
  trg_vec <- sum(dcgset!=0)
  # if there are more than 1 target vector
  if(trg_vec != 0 && trg_vec > 1){
    # all possible vectors from permutation
    possible_vec <- gtools::permutations(2, trg_vec, v = 0:1, repeats.allowed = T)
    for (j in 1:nrow(possible_vec)){
      tmp <- dcgset # copy in tmp
      tmp[which(tmp!=0)] <- possible_vec[j,]
      checktruefalse <- c()
      for(i in seq_len(num_node)){
        # think again about using "type = cpdag" (what's the implication)
        # check if all definite ancestors with "-" arrow tail is included in the pool
        checktruefalse[i] <- (all(searchAM_KP(pag, i, type="an") %in% possAn(tmp, i, type="pdag")) &&
        # check if the ancestors in cpdag are part of all possible ancestors considering circles
                                all(possAn(tmp, i, type="pdag") %in% searchAM_KP(pag, i, type="ant")))
      }
      # if ancestral relations are identical for all nodes
      if(all(checktruefalse)){
        # then save it
        equivset[[j]] <- tmp
        # otherwise
      } else {
        # assign NULL
        equivset[[j]] <- NULL
      }
    }
  }
  # if there is only 1 target
  if(trg_vec == 1){
    # possible vectors
    tmp[which(tmp!=0)] <- 0
    # do the same checks as before
    checktruefalse <- c()
    for(k in seq_len(num_node)){
      checktruefalse[i] <- (all(searchAM_KP(pag, i, type="an") %in% possAn(tmp, i, type="cpdag")) &&
                              all(possAn(tmp, i, type="cpdag") %in% searchAM_KP(pag, i, type="ant")))
    }

    if(all(checktruefalse)){
      equivset[[k]] <- tmp
    } else {
      equivset[[k]] <- NULL
    }
  }
  # exclude null values
  ## output after checking ancestral relationship (step1)
  equivset <- equivset[!sapply(equivset, is.null)]


  ## STEP 2) if solid-underline exists, then
  # A- _B_ - C : B has to be an ancestor of at least either A or C
  if(nrow(extracttriples(ccdobj$graph$getUnderLines()$toString()))!=0){
    # extract first, second, third nodes in the triples
    solidline <- extracttriples(ccdobj$graph$getUnderLines()$toString())
    # get the nodes name
    vars <- ccdobj$nodes
    removelist <- list()
    for(e in 1:nrow(solidline)){
      # extract first, second, third nodes in the triple
      fn <- solidline[e,]$fistnode
      mn <- solidline[e,]$middlenode
      ln <- solidline[e,]$lastnode

      # go through the current equivsets
      ancestorno <- c()
      for (q in 1:length(equivset)){
        equivgraph <- equivset[[q]]
        # find the ancestors
        ind_fn <- which(vars %in% fn)
        ind_ln <- which(vars %in% ln)
        fn_an <- possAn(equivgraph, ind_fn, type="pdag")
        ln_an <- possAn(equivgraph, ind_ln, type="pdag")
          # if mn is ancestor of either of them
          if (!(mn %in% names(equivgraph[,fn]) | mn %in% names(equivgraph[,ln]) )){
            # then record the index
            ancestorno[q] <- q
          }
        removelist[[e]] <- ancestorno
      }
    }
    if(length(removelist)!= 0){
    # remove the DCGs that are in the removelist
    removelist <- na.omit(unique(unlist(removelist)))
    equivset <- equivset[-removelist]
    }
  }


  ## STEP 3) if dotted-underline exists, then
  # A- ..B.. - C : B cannot be a descendant of common child of A and C
  if(nrow(extracttriples(ccdobj$graph$getDottedUnderlines()$toString()))!=0){
    # extract first, second, third nodes in the triples
    dotted <- extracttriples(ccdobj$graph$getDottedUnderlines()$toString())
    # get the nodes name
    vars <- ccdobj$nodes
    removelist <- list()
    for(w in 1:nrow(dotted)){
      # extract first, second, third nodes in the triple
      fn <- dotted[w,]$fistnode
      mn <- dotted[w,]$middlenode
      ln <- dotted[w,]$lastnode
      # go through the current equivsets
      descendtyesno <- c()
      for (p in 1:length(equivset)){
        equivgraph <- equivset[[p]]
        # find the common child
        fn_child <- names(equivgraph[,fn][equivgraph[,fn] == 1])
        ln_child <- names(equivgraph[,ln][equivgraph[,ln] == 1])
        commonchild <- intersect(fn_child, ln_child)
        if (length(commonchild) != 0){
          # get the index of common child
          ind_ch <- which(vars %in% commonchild) # vars are the nodes defined earlier
          # check if the middle node is the descendant of common child
          de_ch <- list()
          for (u in 1:length(ind_ch)){
            de_ch[[u]] <- vars[possDe(equivgraph, ind_ch[u], type = "pdag")] ## what is the implication of using "type=pdag" here?
          }
          # descendants of common child
          de_ch <- unique(unlist(de_ch))
          # if it is then remove the corresponding graph
          if (mn %in% de_ch){
            descendtyesno[p] <- p
          }
        }
        removelist[[w]] <- descendtyesno
      }
    }
    # remove the DCGs that are in the removelist
    removelist <- na.omit(unique(unlist(removelist)))
    equivset <- equivset[-removelist]
  }

  if(length(equivset)==0) stop("nothing left in equivalnce class after step2 (checking underlinined nodes).")

  ## STEP 4) if there is an edge in PAG A *-*B --> A and B are d-connected in all DCGs
  # check if the p-adjacency is preserved as d-connected in circle-circle edges
  circleind <- which(pag==1, arr.ind=T)
  if(length(circleind)!=0){
    circlecircle <- list() # storage for circle-circle edge
    for(o in 1:nrow(circleind)){
      # see if [i,j] == [j,i] == 1
      if(any(apply(apply(circleind, 1, function(x) x == circleind[o,][c("col", "row")]), 2, sum) == 2)){
        circlecircle[[o]] <- circleind[o,]
      }
    }
    # indices for circle-circle
    ind <- do.call(rbind, circlecircle)
    # go through the equivset to check circle-circle edges are d-connected
    for(g in 1:length(equivset)){  # change the iterative symbol
      candidate <- equivset[[g]]
      vars <- ccdobj$nodes
      for(h in 1:nrow(ind)){
        # x-y needs to be d-connected
        x <- vars[ind[h,1]]
        y <- vars[ind[h,2]]
        # rest of vars except for x and y
        rest <- setdiff(vars, c(x, y))
        # possible sepsets
        for(len in 1:length(rest)){
          sep <- arrangements::combinations(rest, k = len, freq = table(rest))
          # if there is a dsep set then remove the graph
          for (r in 1:nrow(sep)){
            # if there is a dsep set then remove the graph
            if(dsep_KP(x, y, S=sep[r,], graph = candidate)){
              equivset[g] <- NA
            }
          }
        }
      }
    }
    equivset <- equivset[!is.na(equivset)]
    if(length(equivset)==0) stop("nothing left in equivalnce class after step3 (checking p-adjacency == d-adjacency).")
  }
  return(equivset)
}



######################
## TO-DOs in future ##
######################
# - using "cpdag" or "pdag" when searching ancestors/descendants... --> think of the implications!!
# - p-adjacency preserved in all DCGs... is it sufficient to check only circle-circle edges. think about it again. (the intuition is when there is no circle, and the endpoints are definitely specified, the ancestor information seems to be sufficient to obtain the equivalence class. But when the circle occurs, suddenly the ancestor information becomes much less informative).
