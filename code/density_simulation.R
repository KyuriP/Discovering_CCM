## =============================================================================
## Description
#
# This script contains code for simulating the high density conditions 
# to see what happens with the model in terms of the uncertainty rate and 
# equilibrium state. 
# It includes one function `isCyclic` which is used in the simulation.
## =============================================================================


## ========================
## Preparation
## ========================

## source the necessary functions
source("code/R/CCD_fnc.R")
source("code/R/data_generating_fnc.R")

## load necessary packages
library(magrittr)
library(purrr)
library(dplyr)
library(qgraph)
library(ggplot2)
library(pcalg)
library(svMisc) # for the progress bar

## slightly modified CCI package
# remotes::install_github("KyuriP/CCI_KP")
library(CCI.KP)

# remotes::install_github("bd2kccd/r-causal")
library(rcausal)


#' Check graph is cyclic
#' @param graph adjacency matrix of a graph
#'
#' @return TRUE/FALSE
isCyclic <- function(graph)
{
  p = nrow(graph);
  visited = rep(FALSE, p);
  recStack = rep(FALSE, p);
  
  for (i in 1:p){
    if (isCyclicUtil(i, visited, recStack, graph)){
      return(TRUE);}
  }
  
  return(FALSE);
}


## ========================
## Simulation
## ========================

# set the seed
# set.seed(123) # you can try out many other seeds

# specify the conditions
p <- 10 # 10 variable model
alpha <- 0.05 
storage <- list()
equilibrium <- c()

for (i in 1:(2*(p-1))){ # 2*(p-1) is the max number of expected neighborhood size E(N)
  en <- i
  DCG=list();
  adjmat=matrix(0,p,p);
  N = p*p - p;
  # Till it is cyclic, keep randomly sampling the graph structure
  while(!isCyclic(adjmat>0)){
    # sample from Bern(en/2(p-1))
    samplesB = rbinom(N, 1, en/(2*(p-1)));
    
    adjmat=matrix(0,p,p);
    # fill in the matrix
    adjmat[upper.tri(adjmat, diag=FALSE)] <- samplesB[1:(N/2)];
    adjmat[lower.tri(adjmat, diag=FALSE)] <- samplesB[(N/2+1):N];
  }
  # sample weights from unif[-0.6, -0.1] & unif[0.1, 0.6]
  weights = matrix((0.5*runif(p^2)+0.1)*sample(c(-1,1),p^2,replace=TRUE),p,p)
  # assign weights and complete B matrix
  adjmat[adjmat != 0] <- weights[adjmat != 0] 
  # specify dim names
  dimnames(adjmat) <- list(paste0("X",1:p), paste0("X",1:p))
  # complete B matrix
  B <- adjmat
  # check equilibrium condition
  equilibrium[i] <- all(abs(Re(eigen(B)$values)) < 1)
  # generate data 
  Bdat <- gen_dat(B, N = 1e5) # N should be large enough to pick up on small edges
  # run algorithms
  CCD  <- ccdKP(df = Bdat, dataType = "continuous", alpha = alpha) %>% 
    CreateAdjMat(length(.$nodes))
  
  FCI <- fci(list(C = cor(Bdat), n = nrow(Bdat)), indepTest=gaussCItest,
             alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
             labels = colnames(Bdat)) %>% .@amat
  
  CCI <- cci(list(C = cor(Bdat), n = nrow(Bdat)), gaussCItest, alpha = alpha, 
             labels = colnames(Bdat), p = ncol(Bdat)) %>% .$maag 
  # update the storage with the number of circles
  storage[[i]] <- list(CCD, FCI, CCI) %>% 
    map_dbl(~sum(.x == 1))
  if(! equilibrium[i]) storage[[i]] <- NA
  # progress bar 
  progress(i, 2*(p-1))
}

# results
storage %>% 
  do.call("rbind",.) %>% 
  as_tibble() %>% 
  mutate(equilibrium) %>% 
  rename(CCD = V1, FCI = V2, CCI = V3)

