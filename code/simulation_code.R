## =============================================================================
## Description
#
# This script contains code for generating data from different models 
# that are used in the main simulation study.
# There are in total 8 models considered: 5nodes-sparse, 5nodes-dense, 10nodes-sparse,
# 10nodes-dense, 5nodes-sparse with latent variables (LV), 5nodes-dense with LVs,
# 10nodes-sparse with LVs, 10nodes-dense with LVs.
#
# We generate 500 datasets from each model and estimate partial ancestral graphs (PAGs)
# using CCD, FCI, and CCI algorithm.
#
# Running algorithms (except for CCD) and plotting the resulting PAGs are commented out
# due to a long running time. if interested, simply uncomment them and run the code.
## =============================================================================


## ========================
## Preparation
## ========================

## source all the necessary functions
source("code/R/CCD_fnc.R")
source("code/R/plot_fnc.R")
source("code/R/searchAM_KP_fnc.R")
source("code/R/data_generating_fnc.R")
source("code/R/eval_metrics.R")
source("code/R/true_ancestral.R")

## load necessary packages
library(magrittr)
library(purrr)
library(furrr)
library(dplyr)
library(qgraph)
library(ggplot2)
library(ggpubr)
library(pcalg)
library(rJava)
library(usethis)
library(devtools)

## slightly modified CCI package
# remotes::install_github("KyuriP/CCI_KP")
library(CCI.KP)

# remotes::install_github("bd2kccd/r-causal")
library(rcausal)

# bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rgraphviz")
BiocManager::install("graph")
BiocManager::install("RBGL")


## set the seed
set.seed(123)

## allow parallel processing
plan(multisession) 

## ====================
## 5p - sparse
## ====================
# specify B matrix
B5sparse = matrix(c(0, 0, 0, 0, 0,
                    1, 0, 0.8, 0, 0,
                    0, 0, 0, 0.9, 0,
                    0, 0.7, 0, 0, 1.5,
                    0, 0, 0, 0, 0), 5, 5, byrow = T)
colnames(B5sparse) <- c("X1", "X2", "X3", "X4", "X5")

# specify layout
layout5 = matrix(c(0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1),5,2,byrow = T)

par(mfrow=c(1,2))
# true 5p sparse DCG
true5psparse <- qgraph(t(B5sparse), layout=layout5, labels = colnames(B5sparse), theme="colorblind")
# equilibrium check
equilibrium_check(B5sparse)

## Data generating
# specify the sample sizes
N <- c(50, 150, 500, 1000, 1500, 2000, 3000, 4000, 5000, 10000)
# specify replication number
n <- 500
# specify alpha level
alpha <- 0.05

# generate data
simdata_5psparse <- N %>% future_map(function(z) {
  # later increase n again to 1000 or sth
  replicate(n = n,
            expr = gen_dat(B5sparse, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))

## True ancestral graph (AG)
dcg_5psparse <- matrix(c(0,1,0,0,0,
                         0,0,0,1,0,
                         0,1,0,0,0,
                         0,0,1,0,0,
                         0,0,0,1,0), 5,5,byrow=T)
trueag_5psparse <- true_ancestral(dcg_5psparse, gen_dat(B5sparse), gaussCItest)
dimnames(trueag_5psparse) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
# plot true AG
plotAG(trueag_5psparse)


## Run CCD algorithm
ccd_5psparse <- simdata_5psparse %>%
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
  )
mat_5psparse <- ccd_5psparse %>%
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

######### CCD object is not intact when saving'em as .Rdata   
######### We save adj matrices but not the ccd objects themselves...
# save(ccd_5psparse, file="data/ccd_5psparse.RData")
# save(mat_5psparse, file="data/mat_5psparse.RData")
# load("data/mat_5psparse.RData")
# can we have a more elegant way to combine "graphNEL" plots? ? ? ? ? ?

# plot resulting PAGs
# pag_ccd5psparse <- map2(ccd_5psparse, mat_5psparse,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_5psparse <- simdata_5psparse %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                     labels = colnames(.x)) %>% .@amat 
#   )

# save(fci_5psparse, file="data/fci_5psparse.RData")
load("data/fci_5psparse.RData")

# plot resulting PAGs
# pag_fci5psparse <- fci_5psparse %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_5psparse <- simdata_5psparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=alpha, 
#                     labels = colnames(.x), p = ncol(.x)) %>% .$maag  
#   )
# save(cci_5psparse, file="data/cci_5psparse.RData")
load("data/cci_5psparse.RData")

# plot resulting PAGs
# pag_cci5psparse <- cci_5psparse %>%
#   map_depth(2, ~plotAG(.x))



## ====================
## 5p - dense
## ====================
# specify B matrix
B5dense = matrix(c(0, 0, 0, 0, 0,
                   1.4, 0, 0.8, 0, 0,
                   0, 0, 0, 0.9, 0,
                   0, 0.7, 0, 0, 1,
                   1, 0, 0, 0, 0), 5, 5, byrow = T)
colnames(B5dense) <- c("X1", "X2", "X3", "X4", "X5")
par(mfrow=c(1,2))
true5pdense <- qgraph(t(B5dense), layout=layout5, labels = colnames(B5dense), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5dense)
# generate data (sample size as specified above)
simdata_5pdense <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B5dense, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


# True Ancestral Graph
dcg_5pdense <- matrix(c(0,1,0,0,1,
                        0,0,0,1,0,
                        0,1,0,0,0,
                        0,0,1,0,0,
                        0,0,0,1,0), 5,5,byrow=T)
trueag_5pdense <- true_ancestral(dcg_5pdense, gen_dat(B5dense), gaussCItest)
dimnames(trueag_5pdense) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
plotAG(trueag_5pdense)

## Run CCD algorithm
ccd_5pdense <- simdata_5pdense %>%
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
  )
mat_5pdense <- ccd_5pdense %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_5pdense, file="data/ccd_5pdense.RData")
# save(mat_5pdense, file="data/mat_5pdense.RData")
# load("data/ccd_5pdense.RData")
# par(mfrow=c(2,5))
# pag_ccd5pdense <- map2(ccd_5pdense, mat_5pdense,
#                         ~map2(.x, .y, plotPAG)
#                          )
## Run FCI algorithm
# fci_5pdense <- simdata_5pdense %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # extract amat
#   )
# save(fci_5pdense, file="data/fci_5pdense.RData")
load("data/fci_5pdense.RData")

# plot resulting PAGs
# pag_fci5pdense <- fci_5pdense %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_5pdense <- simdata_5pdense %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag
#   )
# save(cci_5pdense, file="data/cci_5pdense.RData")
load("data/cci_5pdense.RData")

# plot resulting PAGs
# pag_cci5pdense <- cci_5pdense %>%
#   map_depth(2, ~plotAG(.x))



## ====================
## 10p - sparse
## ====================

B10sparse = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0.8, 0, 0, 0, 0, 0, 0, 0, 
                     0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 
                     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0.8, 0, 0.5, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 0, 0, 0.8, 0, 
                     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 10, 10, byrow = T)
dimnames(B10sparse) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))

# specify layout
layout10 = matrix(c(0,1,
                    2,1,
                    1,0,
                    2,-1,
                    3,0,
                    4, -1,
                    5, 0,
                    6, -1,
                    4, 1,
                    7, 1),10,2,byrow = T)
par(mfrow=c(1,2))
true10psparse <- qgraph(t(B10sparse), layout = layout10, labels = colnames(B10sparse), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B10sparse)

# generate data (sample size as specified above)
simdata_10psparse <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B10sparse, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))

## True Ancestral Graph
dcg_10psparse <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                          0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 
                          0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                          0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                          0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                          0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                          0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                          0, 0, 0, 0, 0, 0, 0, 0, 1, 0), 10, 10, byrow = T)
trueag_10psparse <- true_ancestral(dcg_10psparse, gen_dat(B10sparse), gaussCItest)
dimnames(trueag_10psparse) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(trueag_10psparse)

## Run CCD algorithm
ccd_10psparse <- simdata_10psparse %>%
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
  )
mat_10psparse <- ccd_10psparse %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_10psparse, file="data/ccd_10psparse.RData")
# save(mat_10psparse, file="data/mat_10psparse.RData")
# load("data/ccd_10psparse.RData")

# plot resulting PAGs
# pag_ccd10psparse <- map2(ccd_10psparse, mat_10psparse,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_10psparse <- simdata_10psparse %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
#   )
# save(fci_10psparse, file="data/fci_10psparse.RData")
load("data/fci_10psparse.RData")

# plot resulting PAGs
# pag_fci10psparse <- fci_10psparse %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_10psparse <- simdata_10psparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  
#   )
# save(cci_10psparse, file="data/cci_10psparse.RData")
load("data/cci_10psparse.RData")

# plot resulting PAGs
# pag_cci10psparse <- cci_10psparse %>%
#   map_depth(2, ~plotAG(.x))



## ====================
## 10p - dense
## ====================
B10dense = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0.3, 0, 0.8, 0, 0, 0, 0, 0, 0, 0, 
                    0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 1, 0, 0, 0.5, 0, 0, 0, 0, 
                    0, 0, 0.6, 1, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0.9, 0, 0.5, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0.8, 1, 
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 10, 10, byrow = T)
colnames(B10dense) <- c(paste("X", 1:10, sep=""))
# specify layout
layout10 = matrix(c(0,1,
                    2,1,
                    1,0,
                    2,-1,
                    3,0,
                    4, -1,
                    5, 0,
                    6, -1,
                    4, 1,
                    7, 1),10,2,byrow = T)

par(mfrow=c(1,2))
true10pdense <- qgraph(t(B10dense), layout = layout10, labels = colnames(B10dense), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B10dense)
# generate data (sample size as specified above)
simdata_10pdense <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B10dense, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))

## True Ancestral Graph
dcg_10pdense <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 1, 0, 1, 0), 10, 10, byrow = T)
trueag_10pdense <- true_ancestral(dcg_10pdense, gen_dat(B10dense), gaussCItest)
dimnames(trueag_10pdense) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(trueag_10pdense)

## Run CCD algorithm
ccd_10pdense <- simdata_10pdense  %>%
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
  )
mat_10pdense  <- ccd_10pdense  %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))
# save(ccd_10pdense, file="data/ccd_10pdense.RData")
# save(mat_10pdense, file="data/mat_10pdense.RData")
# load("data/ccd_10pdense.RData")

# plot resulting PAGs
# pag_ccd10pdense <- map2(ccd_10pdense, mat_10pdense,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_10pdense <- simdata_10pdense  %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
#   )
# save(fci_10pdense, file="data/fci_10pdense.RData")
load("data/fci_10pdense.RData")

# plot resulting PAGs
# pag_fci10pdense <- fci_10pdense  %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_10pdense  <- simdata_10pdense %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )
# save(cci_10pdense, file="data/cci_10pdense.RData")
load("data/cci_10pdense.RData")

# plot resulting PAGs
# pag_cci10pdense  <- cci_10pdense  %>%
#   map_depth(2, ~plotAG(.x))



## ====================
## 5p with LV sparse
## ====================
# specify B matrix
B5_lvsparse = matrix(c(0, 0, 0, 0, 0, 1,
                       0, 0, 0.4, 0, 0, 1,
                       0, 0, 0, 0.5, 0,0,
                       0, 0.7, 0, 0, 1.5,0,
                       0, 0, 0, 0, 0,0,
                       0,0,0,0,0,0), 6, 6, byrow = T)
colnames(B5_lvsparse) <- c("X1", "X2", "X3", "X4", "X5", "L1")

# specify layout
layout5_lv = matrix(c(0,1,
                      0,0,
                      1,-1,
                      2,0,
                      2,1,
                      -1, 0.5),6,2,byrow = T)
par(mfrow=c(1,2))
true5p_lvsparse <- qgraph(t(B5_lvsparse), layout=layout5_lv, labels = colnames(B5_lvsparse), theme="colorblind")


## Data generating
# equilibrium check
equilibrium_check(B5_lvsparse)
# generate data (sample size as specified above)
simdata_5pLVsparse <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B5_lvsparse, N = z)[,-6],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## True Ancestral Graph
# [i,j] = [j,i] = 2: a LV exists between i and j
dcg_5psparseLV <- matrix(c(0, 2, 0, 0, 0, 
                           2, 0, 0, 1, 0, 
                           0, 1, 0, 0, 0,
                           0, 0, 1, 0, 0,
                           0, 0, 0, 1, 0), 5, 5, byrow = T)
trueag_5psparseLV <- true_ancestral(dcg_5psparseLV, gen_dat(B5_lvsparse), gaussCItest)
dimnames(trueag_5psparseLV) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
plotAG(trueag_5psparseLV)


## Run CCD algorithm
ccd_5pLVsparse <- simdata_5pLVsparse  %>%
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
  )
mat_5pLVsparse  <- ccd_5pLVsparse  %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))
# save(ccd_5pLVsparse, file="data/ccd_5pLVsparse.RData")
# save(mat_5pLVsparse, file="data/mat_5pLVsparse.RData")
# load("data/ccd_5pLV.RData")

# plot resulting PAGs
# pag_ccd5pLVsparse <- map2(ccd_5pLVsparse, mat_5pLVsparse,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_5pLVsparse <- simdata_5pLVsparse  %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # extract amat
#   )
# save(fci_5pLVsparse, file="data/fci_5pLVsparse.RData")
load("data/fci_5pLVsparse.RData")

# plot resulting PAGs
# pag_fci_5pLVsparse <- fci_5pLVsparse  %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_5pLVsparse  <- simdata_5pLVsparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )
# 
# save(cci_5pLVsparse, file="data/cci_5pLVsparse.RData")
load("data/cci_5pLVsparse.RData")

# plot resulting PAGs
# pag_cci_5pLVsparse <- cci_5pLVsparse  %>%
#   map_depth(2, ~plotAG(.x))



## ====================
## 5p with LV dense
## ====================
B5_lvdense = matrix(c(0, 0, 0, 0, 0, 1,
                      0, 0, 0.4, 0, 0, 1,
                      0, 0, 0, 0.5, 0,0,
                      0, 0.7, 0, 0, 0.5, 0,
                      0.6, 0, 0, 0, 0,0,
                      0,0,0,0,0,0), 6, 6, byrow = T)
colnames(B5_lvdense) <- c("X1", "X2", "X3", "X4", "X5", "L1")
# specify layout
layout5_lv = matrix(c(0,1,
                      0,0,
                      1,-1,
                      2,0,
                      2,1,
                      -1, 0.5),6,2,byrow = T)
par(mfrow=c(1,2))
true5p_lvdense <- qgraph(t(B5_lvdense), layout=layout5_lv, labels = colnames(B5_lvdense), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5_lvdense)
# generate data (sample size as specified above)
simdata_5pLVdense <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B5_lvdense, N = z)[,-6],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))
## True Ancestral Graph
# [i,j] = [j,i] = 2: a LV exists between i and j
dcg_5pdenseLV <- matrix(c(0, 2, 0, 0, 1, 
                          2, 0, 0, 1, 0, 
                          0, 1, 0, 0, 0,
                          0, 0, 1, 0, 0,
                          0, 0, 0, 1, 0), 5, 5, byrow = T)
trueag_5pdenseLV <- true_ancestral(dcg_5pdenseLV, gen_dat(B5_lvdense), gaussCItest)
dimnames(trueag_5pdenseLV) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
plotAG(trueag_5pdenseLV)

## Run CCD algorithm
ccd_5pLVdense <- simdata_5pLVdense  %>%
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
  )
mat_5pLVdense  <- ccd_5pLVdense  %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))
# save(ccd_5pLVdense, file="data/ccd_5pLVdense.RData")
# save(mat_5pLVdense, file="data/mat_5pLVdense.RData")
# load("data/ccd_5pLV.RData")
# par(mfrow=c(2,5))
# pag_ccd5pLVdense <- map2(ccd_5pLVdense, mat_5pLVdense,
#                         ~map2(.x, .y, plotPAG)
#                          )
## Run FCI algorithm
# fci_5pLVdense <- simdata_5pLVdense  %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat
#   )

# save(fci_5pLVdense, file="data/fci_5pLVdense.RData")
load("data/fci_5pLVdense.RData")

# plot resulting PAGs
# pag_fci5pLVdense <- fci_5pLVdense  %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_5pLVdense  <- simdata_5pLVdense %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )

# save(cci_5pLVdense, file="data/cci_5pLVdense.RData")
load("data/cci_5pLVdense.RData")

# plot resulting PAGs
# pag_cci5pLVdense <- cci_5pLVdense  %>%
#   map_depth(2, ~plotAG(.x))



## ====================
## 10p with LV sparse
## ====================
# specify B matrix

## with 2 LVs
B10_lvsparse = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0.7,
                        0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.7,
                        0, 0, 0, 0, 0.8, 0, 0.5, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0.8, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0.8, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 12, 12, byrow = T)
# colnames(B10_lvsparse) <- c(paste("X", 1:10, sep=""), "L1")
colnames(B10_lvsparse) <- c(paste("X", 1:10, sep=""), "L1", "L2")
# specify layout
layout10LV2 = matrix(c(0, 1,
                       2, 1,
                       1, 0,
                       2, -1,
                       3, 0,
                       4, -1,
                       5, 0,
                       6, -1,
                       4, 1,
                       7, 1,
                       8, 0,
                       3, 2), 12, 2, byrow = T)
par(mfrow=c(1,2))
true10pLVsparse <- qgraph(t(B10_lvsparse), layout = layout10LV2, labels = colnames(B10_lvsparse), theme="colorblind")
## Data generating
# equilibrium check
equilibrium_check(B10_lvsparse)
# generate data (sample size as specified above)
simdata_10pLVsparse <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B10_lvsparse, N = z)[,-c(11,12)],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))
## True Ancestral Graph
# [i,j] = [j,i] = 2: a LV exists between i and j
dcg_10psparseLV <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 
                            0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 
                            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                            0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 
                            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                            0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 
                            0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                            0, 0, 0, 0, 0, 0, 0, 2, 1, 0), 10, 10, byrow = T)
trueag_10psparseLV <- true_ancestral(dcg_10psparseLV, gen_dat(B10_lvsparse), gaussCItest)
dimnames(trueag_10psparseLV) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(trueag_10psparseLV)

## Run CCD algorithm
ccd_10pLVsparse  <- simdata_10pLVsparse   %>%
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
  )
mat_10pLVsparse   <- ccd_10pLVsparse %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))
# save(ccd_10pLVsparse, file="data/ccd_10pLVsparse.RData")
# save(mat_10pLVsparse, file="data/mat_10pLVsparse.RData")
# load("data/ccd_10pLVsparse.RData")
# par(mfrow=c(2,5))
# pag_ccd10pLVsparse <- map2(ccd_10pLVsparse, mat_10pLVsparse,
#                         ~map2(.x, .y, plotPAG)
#                          )
## Run FCI algorithm
# fci_10pLVsparse  <- simdata_10pLVsparse   %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                     labels = colnames(.x)) %>% .@amat  
#   )
# save(fci_10pLVsparse, file="data/fci_10pLVsparse.RData")
load("data/fci_10pLVsparse.RData")

# plot resulting PAGs
# pag_fci10pLV  <- fci_10pLVsparse %>%
#   map_depth(2, ~plotAG(.x))
## Run CCI algorithm
# cci_10pLVsparse  <- simdata_10pLVsparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )
# save(cci_10pLVsparse, file="data/cci_10pLVsparse.RData")
load("data/cci_10pLVsparse.RData")

# plot resulting PAGs
# pag_cci10pLVsparse   <- cci_10pLVsparse %>%
#   map_depth(2, ~plotAG(.x))



## ====================
## 10p with LV dense
## ====================
# specify B matrix

# B10_lvdense = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                        0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4,
#                        0.4, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                        0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 0, 0,
#                        0, 0, 0.6, 1, 0, 0, 0, 0, 0, 0, 0, 0.5,
#                        0, 0, 0, 0, 0.8, 0, 0.5, 0, 0, 0, 0, 0,
#                        0, 0, 0, 0, 0, 0, 0, 1, 0.8, 0.6, 0, 0,
#                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.8, 0,
#                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4, 0, 0,
#                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0,
#                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 12, 12, byrow = T)
B10_lvdense = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0.5, 0, 1.4, 0, 0, 0, 0, 0, 0, 0, 0, 1.1,
                       0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 1, 0, 0, 0.9, 0, 0, 0, 0, 0, 0,
                       0, 0, 1.2, 0.8, 0, 0, 0, 0, 0, 0, 0, 0.7,
                       0, 0, 0, 0, 0.8, 0, 0.5, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 1.4, 0.5, 0.4, 0,
                       0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0.6, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 12, 12, byrow = T)
colnames(B10_lvdense) <- c(paste("X", 1:10, sep=""), "L1", "L2")
# specify layout
layout10LV2 = matrix(c(0, 1,
                       2, 1,
                       1, 0,
                       2, -1,
                       3, 0,
                       4, -1,
                       5, 0,
                       6, -1,
                       4, 1,
                       7, 1,
                       8, 0,
                       3, 2), 12, 2, byrow = T)
par(mfrow=c(1,2))
true10pLVdense <- qgraph(t(B10_lvdense), layout = layout10LV2, labels = colnames(B10_lvdense), theme="colorblind")
## Data generating
# equilibrium check
equilibrium_check(B10_lvdense)
# generate data (sample size as specified above)
simdata_10pLVdense <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B10_lvdense, N = z)[,-c(11,12)],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))
## True Ancestral Graph
# [i,j] = [j,i] = 2: a LV exists between i and j
dcg_10pdenseLV <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 
                           0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                           0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 1, 2, 1, 0), 10, 10, byrow = T)
trueag_10pdenseLV <- true_ancestral(dcg_10pdenseLV, gen_dat(B10_lvdense), gaussCItest)
dimnames(trueag_10pdenseLV) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(trueag_10pdenseLV)

## Run CCD algorithm
ccd_10pLVdense  <- simdata_10pLVdense   %>%
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
  )
mat_10pLVdense   <- ccd_10pLVdense %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_10pLVdense, file="data/ccd_10pLVdense.RData")
# save(mat_10pLVdense, file="data/mat_10pLVdense.RData")
# load("data/ccd_10pLV.RData")
# plot resulting PAGs
# pag_ccd10pLVdense <- map2(ccd_10pLVdense, mat_10pLVdense,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_10pLVdense  <- simdata_10pLVdense   %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                     labels = colnames(.x)) %>% .@amat
#   )
# save(fci_10pLVdense, file="data/fci_10pLVdense.RData")
load("data/fci_10pLVdense.RData")

# plot resulting PAGs
# pag_fci10pLVdense <- fci_10pLVdense   %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_10pLVdense  <- simdata_10pLVdense  %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )
# save(cci_10pLVdense, file="data/cci_10pLVdense.RData")
load("data/cci_10pLVdense.RData")

# plot resulting PAGs
# pag_cci10pLVdense   <- cci_10pLVdense %>%
#   map_depth(2, ~plotAG(.x))



## ====================
## Running time 
## ====================

#remotes::install_github("joshuaulrich/microbenchmark")

library(microbenchmark)

times <- microbenchmark(
  ccd_5psparse = ccdKP(df=simdata_5psparse[[1]][[1]], dataType = "continuous", alpha = 0.05),
  fci_5psparse = fci(list(C = cor(simdata_5psparse[[1]][[1]]), n = 1e3),indepTest=gaussCItest, alpha = 0.05, selectionBias= FALSE, labels = colnames(simdata_5psparse[[1]][[1]])),
  cci_5psparse = cci(list(C = cor(simdata_5psparse[[1]][[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_5psparse[[1]][[1]])),

  ccd_5pdense = ccdKP(df=simdata_5pdense[[1]][[1]], dataType = "continuous", alpha = 0.05),
  fci_5pdense = fci(list(C = cor(simdata_5pdense[[1]][[1]]), n = 1e3),indepTest=gaussCItest, alpha = 0.05, selectionBias= FALSE, labels = colnames(simdata_5pdense[[1]][[1]])),
  cci_5pdense = cci(list(C = cor(simdata_5pdense[[1]][[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_5pdense[[1]][[1]])),

  ccd_10psparse = ccdKP(df=simdata_10psparse[[1]][[1]], dataType = "continuous", alpha = 0.05),
  fci_10psparse = fci(list(C = cor(simdata_10psparse[[1]][[1]]), n = 1e3),indepTest=gaussCItest, alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(simdata_10psparse[[1]][[1]])),
  cci_10psparse = cci(list(C = cor(simdata_10psparse[[1]][[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_10psparse[[1]][[1]])),

  ccd_10pdense = ccdKP(df=simdata_10pdense[[1]][[1]], dataType = "continuous", alpha = 0.05),
  fci_10pdense = fci(list(C = cor(simdata_10pdense[[1]][[1]]), n = 1e3),indepTest=gaussCItest, alpha = 0.05, selectionBias= FALSE, labels = colnames(simdata_10pdense[[1]][[1]])),
  cci_10pdense = cci(list(C = cor(simdata_10pdense[[1]][[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_10pdense[[1]][[1]])),

  ccd_5pLVsparse = ccdKP(df=simdata_5pLVsparse[[1]][[1]], dataType = "continuous", alpha = 0.05),
  fci_5pLVsparse = fci(list(C = cor(simdata_5pLVsparse[[1]][[1]]), n = 1e3),indepTest=gaussCItest,
                alpha = 0.05, selectionBias= FALSE, labels = colnames(simdata_5pLVsparse[[1]][[1]])),
  cci_5pLVsparse = cci(list(C = cor(simdata_5pLVsparse[[1]][[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_5pLVsparse[[1]][[1]])),

  ccd_5pLVdense = ccdKP(df=simdata_5pLVdense[[1]][[1]], dataType = "continuous", alpha = 0.05),
  fci_5pLVdense = fci(list(C = cor(simdata_5pLVdense[[1]][[1]]), n = 1e3),indepTest=gaussCItest,
                       alpha = 0.05, selectionBias= FALSE, labels = colnames(simdata_5pLVdense[[1]][[1]])),
  cci_5pLVdense = cci(list(C = cor(simdata_5pLVdense[[1]][[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_5pLVdense[[1]][[1]])),
  
  ccd_10pLVsparse = ccdKP(df=simdata_10pLVsparse[[1]][[1]], dataType = "continuous", alpha = 0.05),
  fci_10pLVsparse = fci(list(C = cor(simdata_10pLVsparse[[1]][[1]]), n = 1e3),indepTest=gaussCItest,
                alpha = 0.05, selectionBias= FALSE, labels = colnames(simdata_10pLVsparse[[1]][[1]])),
  cci_10pLVsparse = cci(list(C = cor(simdata_10pLVsparse[[1]][[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_10pLVsparse[[1]][[1]])),
  
  ccd_10pLVdense = ccdKP(df=simdata_10pLVdense[[1]][[1]], dataType = "continuous", alpha = 0.05),
  fci_10pLVdense = fci(list(C = cor(simdata_10pLVdense[[1]][[1]]), n = 1e3),indepTest=gaussCItest,
                        alpha = 0.05, selectionBias= FALSE, labels = colnames(simdata_10pLVdense[[1]][[1]])),
  cci_10pLVdense = cci(list(C = cor(simdata_10pLVdense[[1]][[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_10pLVdense[[1]][[1]]))
  
)

times <- times %>%
  mutate(algorithm = substr(expr, 1, 3),
         condition = stringr::str_split(expr, "_", simplify=T)[,2])

## plot the results
times %>%
  ggplot(aes(x=factor(condition, levels= c("5psparse", "5pdense", "10psparse", "10pdense", "5pLVsparse","5pLVdense", "10pLVsparse","10pLVdense")), y = log(time), col= factor(algorithm))) +
  geom_boxplot(position = "dodge", outlier.size = 0.8, outlier.alpha = 0.2) + theme_classic() +
  # scale_x_discrete(name ="Condition",
  #                  labels=c("", "5p-sparse", "", "","5p-dense","","", "10p-sparse","","","10p-dense","","","5p-LV","","","10p-LV","")) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(y = " log(ms)", x = "conditions", title = "Algorithm Running Time", subtitle = "Time in milliseconds (ms)") +
  theme(axis.text.x = element_text(face = "bold", angle=40, margin = margin(t = 13)))


