## =============================================================================
## Description
#
# This script contains six functions for computing evaluation metrics 
# to assess algorithm performance: 
# `precision_recall`, `uncertainty`,`SHD`,
# `high_freq`, `prop_correct`, and `prop_uncertain`.
# 
# Each function computes the corresponding evaluation metrics.
# See below for the detailed description of each function.
## =============================================================================

# load package
library(magrittr) # for pipes used in functions


#' Compute precision and recall
#' @param truepag adjacency matrix of true ancestral graph
#' @param estimatedpag estimated pag adjacency matrix
#'
#' @return precision and recall value for all three categories of edge-endpoint 
#' (e.g., arrowhead, arrowtail, null)
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
  }
  
  if(2 %in% estimatedpag) {
    pre_head = cm["2","2"]/sum(cm[,"2"])
    rec_head = cm["2","2"]/sum(cm["2",])
  } else {
    pre_head = 0; rec_head = 0;
  }
  
  if(3 %in% estimatedpag) {
    pre_tail = cm["3","3"]/sum(cm[,"3"])
    rec_tail = cm["3","3"]/sum(cm["3",])
  } else {
    pre_tail = 0; rec_tail = 0;
  }
  # average precision and recall
  avg_pre = sum(pre_null, pre_tail, pre_head) / 3
  avg_rec = sum(rec_null, rec_tail, rec_head) / 3
  # return precision and recall value for each edge-endpoint along with the average
  return(list(precision_null = pre_null,
              precision_tail = pre_tail,
              precision_head = pre_head,
              average_precision = avg_pre,
              recall_null = rec_null,
              recall_tail = rec_tail,
              recall_head = rec_head,
              average_recall = avg_rec))
}




#' Compute uncertainty rate
#' @param estimatedpag estimated pag adjacency matrix
#'
#' @return uncertainty rate (proportion of circle edge-endpoint in the corresponding graph)
uncertainty <- function(estimatedpag){
  # p = number of variables
  p <- ncol(estimatedpag)
  # get the number of circle endpoint
  n_circ <- sum(estimatedpag == 1) 
  # compute the uncertainty rate = number of circle / total number of possible edge-endpoint
  uncertainty_rate <- n_circ / (p * (p-1)) # total number of possible edge-endpoint = choose(p,2) * 2
  # number of circle edge-endpoints
  # uncertainty_rate <- n_circ  
  return(uncertainty_rate)
}




#' Compute structural Hamming distance (SHD)
#' @param truepag adjacency matrix of true ancestral graph
#' @param estimatedpag estimated pag adjacency matrix
#' 
#' @return structural Hamming distance (int)
SHD <- function(truepag, estimatedpag){
  # g1 = estimated graph (graph matrix1: m1)
  # g2 = true graph (graph matrix2: m2)
  m1 <- estimatedpag
  m2 <- truepag
  
  # adjacency matrix
  m1[m1 != 0] <- 1
  m2[m2 != 0] <- 1
  
  # container for the shd value: starts from 0
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




#' Create an adjacency matrix consisting of the most frequent edge-endpoint
#'
#' @param amat list of adjacency matrices
#' @param p number of variables (dimension of amat = p by p)
#' 
#' @return adjacency matrix consist of the most frequent values (edge-endpoint)
high_freq <- function(amat, p){
  # highest freq value storage
  freq <- matrix(NA, p, p)
  # loop through every row & column
  for(i in 1:p){
    for(j in 1:p){
      # compute the frequency
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
#' 
#' @return matrix containing the proportion of correct estimation per each edge-endpoint
prop_correct <- function(amat, truemat, p){
  correct <-  amat %>% 
        # compare to the true adjacency matrix 
        purrr::map(~.x == truemat) %>% 
        # count the number of "TRUE"s and divide by the total number of amats
        Reduce("+", .) / length(amat)
  # return the proportion of correct estimation out of the whole estimated sets
  return(correct)
}




#' Compute the proportion of circle (uncertainty) estimation per cell 
#'
#' @param amat list of adjacency matrices
#' @param p number of variables (dimension of amat = p by p)
#' 
#' @return matrix containing the proportion of uncertainty per each edge-endpoint
prop_uncertain <- function(amat, p){
  uncertain <- amat %>% 
    # check if it is a circle endpoint
    purrr::map(~.x  == 1) %>% 
    # count the number of "TRUE"s and divide by the total number of amats
    Reduce("+", .) / length(amat)
  # return the proportion of circle endpoint (i.e., uncertainty rate) out of the whole estimated sets
  return(uncertain)
}




## ============================================================================
# Extra: two functions below `precision2` and `recall2` are essentially the same 
# as `precision_recall` function but only return average precision and recall 
# value, respectively. They are used in the extended simulation study 
# where we randomly sample B matrices. 
## ============================================================================

#' Compute only average precision 
#' @param truepag adjacency matrix of true ancestral graph
#' @param estimatedpag estimated pag adjacency matrix
#'
#' @return average precision 
#' (e.g., arrowhead, arrowtail, null)
precision2 <- function(truepag, estimatedpag){
  ## compare to the true pag
  # create a confusion matrix
  cm <- table(true = truepag, est = estimatedpag, exclude="1") # exclude circle
  # compute precision & recall for each category
  # (0 = null, 1 = circle, 2 = arrow head, 3 = arrow tail)
  if(0 %in% estimatedpag) {
    pre_null = cm["0","0"]/sum(cm[,"0"])
  } else {
    pre_null = 0;
  }
  
  if(2 %in% estimatedpag) {
    pre_head = cm["2","2"]/sum(cm[,"2"])
  } else {
    pre_head = 0; 
  }
  
  if(3 %in% estimatedpag) {
    pre_tail = cm["3","3"]/sum(cm[,"3"])
  } else {
    pre_tail = 0; 
  }
  # if there is any NaN, replace with zero (it is due to the zero in the denominator)
  tmp <- c(pre_null, pre_head, pre_tail)
  tmp[is.nan(tmp)] <- 0
  pre_null <- tmp[1]; pre_head <- tmp[2]; pre_tail <- tmp[3]
  # average precision and recall
  avg_pre = sum(pre_null, pre_tail, pre_head) / 3
  # return average precision value
  return(average_precision = avg_pre)
}




#' Compute only average recall 
#' @param truepag adjacency matrix of true ancestral graph
#' @param estimatedpag estimated pag adjacency matrix
#'
#' @return average recall 
#' (e.g., arrowhead, arrowtail, null)
recall2 <- function(truepag, estimatedpag){
  ## compare to the true pag
  # create a confusion matrix
  cm <- table(true = truepag, est = estimatedpag, exclude="1") # exclude circle
  # compute precision & recall for each category
  # (0 = null, 1 = circle, 2 = arrow head, 3 = arrow tail)
  if(0 %in% estimatedpag) {
    rec_null = cm["0","0"]/sum(cm["0",])
  } else {
    rec_null = 0;
  }
  
  if(2 %in% estimatedpag) {
    rec_head = cm["2","2"]/sum(cm["2",])
  } else {
    rec_head = 0;
  }
  
  if(3 %in% estimatedpag) {
    rec_tail = cm["3","3"]/sum(cm["3",])
  } else {
    rec_tail = 0;
  }
  # if there is any NaN, replace with zero (it is due to the zero in the denominator)
  tmp <- c(rec_null, rec_head, rec_tail)
  tmp[is.nan(tmp)] <- 0
  rec_null <- tmp[1]; rec_head <- tmp[2]; rec_tail <- tmp[3]
  
  # average recall
  avg_rec = sum(rec_null, rec_tail, rec_head) / 3
  # return average recall value 
  return(average_recall = avg_rec)
}
