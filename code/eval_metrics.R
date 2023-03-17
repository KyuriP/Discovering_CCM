
precision_recall <- function(truepag, estimatedpag){
  ## compare to the true pag
  # create a confusion matrix
  cm <- table(true = truepag, est = estimatedpag, exclude="1") # exclude circle
  # compute precision & recall for each category
  # (0 = null, 1 = circle, 2 = arrow head, 3 = arrow tail)
  if(0 %in% estimatedpag) {
    pre_null = cm["0","0"]/sum(cm[,"0"])
    rec_null = cm["0","0"]/sum(cm["0",])
  } else {
    pre_null = 0; rec_null = 0;
    #print("No null endpoint occurred in the estimated graph.")
  }
  
  if(2 %in% estimatedpag) {
    pre_head = cm["2","2"]/sum(cm[,"2"])
    rec_head = cm["2","2"]/sum(cm["2",])
  } else {
    pre_head = 0; rec_head = 0;
    #print("No tail endpoint occurred in the estimated graph.")
  }
  
  if(3 %in% estimatedpag) {
    pre_tail = cm["3","3"]/sum(cm[,"3"])
    rec_tail = cm["3","3"]/sum(cm["3",])
  } else {
    pre_tail = 0; rec_tail = 0;
    #print("No head endpoint occurred in the estimated graph.")
  }
  # average precision and recall
  avg_pre = sum(pre_null, pre_tail, pre_head) / 3
  avg_rec = sum(rec_null, rec_tail, rec_head) / 3
  
  return(list(precision_null = pre_null,
              precision_tail = pre_tail,
              precision_head = pre_head,
              average_precision = avg_pre,
              recall_null = rec_null,
              recall_tail = rec_tail,
              recall_head = rec_head,
              average_recall = avg_rec))
}





uncertainty <- function(estimatedpag){
  p <- ncol(estimatedpag)
  n_circ <- sum(estimatedpag == 1) 
  uncertainty_rate <- n_circ / (p * (p-1)) # total number of possible edgeendpoint = choose(p,2) * 2
  # uncertainty_rate <- n_circ  # total number of possible edgeendpoint = choose(p,2) * 2 
  return(uncertainty_rate)
}





SHD <- function(truepag, estimatedpag){
  # g1 = estimated graph (graph matrix1: m1)
  # g2 = true graph (graph matrix2: m2)
  m1 <- estimatedpag
  m2 <- truepag
  
  # adjacency matrix
  m1[m1 != 0] <- 1
  m2[m2 != 0] <- 1
  
  # shd starts with 0
  shd <- 0
  
  ## Remove superfluous edges from g1
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  # get indices for superfluous edges
  ind <- which(ds > 0)
  # remove them from g1
  m1[ind] <- 0
  # add the number of superfluous edges to the SHD value
  shd <- shd + length(ind)/2
  
  ## Add missing edges to g1
  # get indices for missing edges
  ind <- which(ds < 0)
  # add them to g1
  m1[ind] <- m2[ind]
  # add the number of missing edges to the SHD value
  shd <- shd + length(ind)/2
  
  ## Compare orientation
  ind <- which(m1==0)
  # match the skeleton of g1 to that of true graph (g2)
  estimatedpag[ind] <- 0
  # check how many edge-endpoints do not match to the true graph
  d <- sum(estimatedpag != truepag)
  # add the number of unmatching edge-endpoints to the SHD value
  shd <- shd + d
  
  return(shd)
}


#' Create an adjacency matrix consisting of the highest frequency values
#'
#' @param amat list of adjacency matrices
#' @param p number of variables (dimension of amat = p by p)
high_freq <- function(amat, p){
  # highest freq value storage
  freq <- matrix(NA, p, p)
  # loop through every row & column
  for(i in 1:p){
    for(j in 1:p){
      freq[i, j] <- amat %>% 
        purrr::map(~.x[i,j]) %>% 
        unlist %>%  
        table  %>%    
        # find the highest freq val
        which.max %>%  
        names %>% as.numeric
    }
  }
  dimnames(freq) <- list(paste0("X", 1:p), paste0("X", 1:p))
  return(freq)
}


#' Compute the proportion of correct estimation per cell 
#'
#' @param amat list of adjacency matrices
#' @param truemat adjacency matrix of the true model
#' @param p number of variables (dimension of amat = p by p)
prop_correct <- function(amat, truemat, p){
  # highest freq value storage
  correct <-  amat %>% 
        purrr::map(~.x == truemat) %>% 
        # sum them all and divide by the number of amats
        Reduce("+", .) / length(amat)
  return(correct)
}




#' Compute the proportion of circle (uncertainty) estimation per cell 
#'
#' @param amat list of adjacency matrices
#' @param p number of variables (dimension of amat = p by p)
prop_uncertain <- function(amat, p){
  uncertain <- amat %>% 
    purrr::map(~.x  == 1) %>% 
    # sum them all and divide by the number of amats
    Reduce("+", .) / length(amat)
  return(uncertain)
}



