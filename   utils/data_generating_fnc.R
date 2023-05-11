## ============================================================================
## Description
#
# It contains five functions: `equilibrium_check`, `gen_dat`, `generatesimdat`,
# `sampleranB2`, and `sampleranB2_pos`.
#
# - `equilibrium_check`: checks if the specified regression matrix (B) reaches 
#    an equilibrium state, which is necessary for a cyclic model to converge.
# - `gen_dat`: generates data according to the specified B matrix.
# - `generatesimdat`: generates *n* datasets (dependent on gen_dat function).
# - `sampleranB2`: generates data given randomly sampled weights of B matrices.
# - `sampleranB2_pos`: generates data given randomly sampled *positive* weights 
#    of B matrices.
## ============================================================================



#' Check equilibrium state
#'
#' @param B square regression matrix
#'
#' @return message whether equilibrium is met or not.
equilibrium_check <- function(B){
  # if the absolute values of all eigen values of B are smaller than 1,
  # then it reaches an equilibrium state
  if(all(abs(Re(eigen(B)$values)) < 1)){
    print("Equilibirum state has been reached.")
  } else {
    print("Equilibirum state has NOT been reached.")
  }
}



#' Generate data from multivariate standard normal distribution
#'
#' @param B square regression matrix
#' @param N sample size (default = 1e6)
#' @param seed seed number
#'
#' @return matrix of simulated data (N by variable numbers)
gen_dat <- function(B, N = 1e6, seed = NULL){
  dimension <- ncol(B)
  # I - B inverse
  inverse <- solve(diag(dimension) - B)
  # sample errors (std. normal)
  if(!is.null(seed)) set.seed(seed)
  errors <- MASS::mvrnorm(n = N, mu = rep(0, dimension), Sigma = diag(dimension))
  # generate data: (I_B)*errors
  data <- t(apply(errors, 1, function(x) inverse %*% x))
  colnames(data) <-  colnames(B)
  return(data)
}



#' Generate n datasets from multivariate standard normal distribution 
#' (used in sensitivity analyses)
#'
#' @param B square regression matrix
#' @param N sample size  
#' @param LV index of latent variable
#' @param n number of replication (default = 500)
#' @param seed seed number (default = 123)
#'
#' @return list of matrices of simulated data 
generatesimdat <- function(B, N, LV = NULL, n = 500, seed=123){
  # generate data n times for each N
  simdat <- N %>% future_map(function(z) {
    replicate(n = n,
              expr = if(is.null(LV)) gen_dat(B, N = z) else gen_dat(B, N = z)[,-LV],
              simplify = FALSE)
    # set the seed
  }, .options = furrr_options(seed=seed) 
  ) %>% 
    # label it with N
    rlang::set_names(.,  N)
  return(simdat)
}



#' Generate data using randomly sampled parameters of B
#'
#' @param B square regression matrix
#' @param LV index of latent variable
#' @param n number of replication (default = 500)
#' @param seed seed number (default = 123)
#'
#' @return list of matrices of simulated data
sampleranB2 <- function(B, LV=NULL, n = 500, seed=123){
  # set the seed
  set.seed(seed)
  # variable number
  p <- ncol(B)
  # data storage
  simdat <- list()
  # number of diff sample size
  for(i in 1:length(N)){
    simdat[[i]] <- list()
    # number of repetition
    for(j in 1:n){
      # sample random weights for B from uniform distribution (unif[-0.9, -0.1] & unif[0.1, 0.9])
      ranB <- matrix((0.8*runif(p^2)+0.1)*sample(c(-1,1), p^2,replace=TRUE), p, p)
      # assign 0 where B = 0
      ind <- which(B == 0, arr.ind = TRUE)
      ranB[ind] <- 0
      # generate data
      if(is.null(LV)){
        simdat[[i]][[j]] <- gen_dat(ranB, N = N[i])
        colnames(simdat[[i]][[j]]) <- colnames(B)
      } else {
        simdat[[i]][[j]] <- gen_dat(ranB, N = N[i])[,-LV]
        colnames(simdat[[i]][[j]]) <- colnames(B)[-LV]
      }
    }
  }
  return(simdat)
}



#' Generate data using randomly sampled *positive* parameters of B
#'
#' @param B square regression matrix
#' @param LV index of latent variable
#' @param n number of replication (default = 500)
#' @param seed seed number (default = 123)
#'
#' @return list of matrices of simulated data

sampleranB2_pos <- function(B, LV=NULL, n = 500, seed=123){
  # set the seed
  set.seed(seed)
  # variable number
  p <- ncol(B)
  # data storage
  simdat <- list()
  # number of diff sample size
  for(i in 1:length(N)){
    simdat[[i]] <- list()
    # number of repetition
    for(j in 1:n){
      # sample random weights for B from uniform distribution: unif[0.1, 0.8]
      ranB <- matrix((0.7*runif(p^2)+0.1), p, p)
      # assign 0 where B = 0
      ind <- which(B == 0, arr.ind = TRUE)
      ranB[ind] <- 0
      # generate data
      if(is.null(LV)){
        simdat[[i]][[j]] <- gen_dat(ranB, N = N[i])
        colnames(simdat[[i]][[j]]) <- colnames(B)
      } else {
        simdat[[i]][[j]] <- gen_dat(ranB, N = N[i])[,-LV]
        colnames(simdat[[i]][[j]]) <- colnames(B)[-LV]
      }
    }
  }
  return(simdat)
}