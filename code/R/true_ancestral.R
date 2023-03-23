
#' Create true ancestral graph adjacency matrix
#'
#' @param trueadj true directed (cyclic) graph adjacency matrix ([i,j] = 1 means i -> j)
#' @param dat data generated given `trueadj`
#' @param indepTest type of tests used for conditional independence
#' @param alpha level of type-I error rate (default = 1e-10)
#' @param p number of variables (default = nrow of `trueadj`)
#'
#' @details
#' "0": no edge; 2": arrow; "3": tail
#'
#' @return an adjacency matrix of the true ancestral graph for the given directed (cyclic) graph
#' 
true_ancestral <- function(trueadj, dat, indepTest, alpha = 1e-10, p = nrow(trueadj)){
  # 1) estimate a skeleton 
  skel.fit <- CCI.KP::skeleton_new(suffStat = list(C = cor(dat), n = nrow(dat)), indepTest=indepTest,
                           alpha = alpha, p = p)
  skel <- as(skel.fit@graph, "matrix")
  
  # 2) extract ancestral relations
  truemat <- skel
  for (i in 1:p){
    # ancestors
    ancs <- possAn(trueadj, i, possible = FALSE, ds = TRUE,
                   type = c("dag"))
    # descendants
    desc <- possDe(trueadj, i, possible = FALSE, ds = TRUE,
                   type = c("dag"))
    truemat[i, ancs] <- 2 
    truemat[i, desc] <- 3
  }
  # if LV, keep the bidirected edge
  # lv_ind <- which(trueadj == 2, arr.ind = TRUE)
  # truemat[lv_ind] <- 2
  
  # set to zero when no edges in skeleton
  ind <- which(skel == 0, arr.ind = TRUE)
  truemat[ind] <- 0
  dimnames(truemat) <- list(1:p, 1:p)
  
  return(truemat)
}

