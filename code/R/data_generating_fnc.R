## ============================================================================
## It contains two functions: `equilibrium_check` and `gen_dat`.
##
## Purpose: It checks if the specified regression matrix (B) reaches
## an equilibrium state, which is necessary for a cyclic model to converge.
## Then, it generates the data according to the specified B matrix.
## ============================================================================

#' Check equilibrium state
#'
#' @param B the square regression matrix
#'
#' @return print out message whether equilibrium is met.
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
#' @param B the square regression matrix
#' @param N the sample size (default = 1e6)
#' @param seed set the specific seed
#'
#' @return the matrix of data (N by variable numbers)
gen_dat <- function(B, N = 1e6, seed = NULL){
  dimension <- ncol(B)
  # I - B inverse
  inverse <- solve(diag(dimension) - B)
  # sample errors (std. normal)
  errors <- MASS::mvrnorm(n = N, mu = rep(0, dimension), Sigma = diag(dimension))
  # generate data: (I_B)*errors
  data <- t(apply(errors, 1, function(x) inverse %*% x))
  colnames(data) <-  colnames(B)
  return(data)
}
