## =============================================================================
## Description
#
# This script generates data from different models that are used in 
# the main simulation study with fixed B matrices and apply three constraint-
# based algorithms (CCD, FCI, and CCI) to the simulated data.
#
# A total of eight different models are simulated: 5nodes-sparse, 5nodes-dense, 
# 10nodes-sparse, 10nodes-dense, 5nodes-sparse with latent variables (LV), 
# 5nodes-dense with LVs, 10nodes-sparse with LVs, and 10nodes-dense with LVs.
#
# For each model, we generate 500 datasets and 
# estimate partial ancestral graphs (PAGs) using CCD, FCI, and CCI algorithms.
#
# Due to the long processing time, the code that executes the algorithms 
# and generates the PAGs is currently commented out. 
# If interested, simply uncomment those lines and run the code.
## =============================================================================
## The script is organized as follows.
# 0. Preparation: Necessary functions and packages are sourced and loaded. 
#                 Also, the core simulation condition is specified (e.g., sample
#                 size, number of replication, and alpha level).
#
# 1 - 8. Simulation: Data is generated and algorithms are applied to simulated 
#                 data. For each condition, a directed cyclic graph (DCG) and 
#                 the corresponding true ancestral graph (AG) are plotted.
#
# 9. Analysis of algorithm running time (extra): The running time for each 
#                 algorithm is computed. This is an additional analysis included 
#                 as part of the supplementary material.
## =============================================================================


## =============================================================================
## 0. Preparation
## =============================================================================

## source all the necessary functions
source("  utils/CCD_fnc.R")             # function for running CCD algorithm
source("  utils/plot_fnc.R")            # function for plotting PAGs
source("  utils/searchAM_KP_fnc.R")     # function for extracting ancestral relations
source("  utils/data_generating_fnc.R") # function for generating data
source("  utils/eval_metrics_fnc.R")    # function for evaluating performance
source("  utils/true_ancestral_fnc.R")  # function for generation AGs

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
# slightly modified CCI package
# remotes::install_github("KyuriP/CCI_KP")
library(CCI.KP)
# remotes::install_github("bd2kccd/r-causal")
library(rcausal)
# bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("Rgraphviz")
# BiocManager::install("graph")
# BiocManager::install("RBGL")


## simulation specification
# specify the sample sizes
N <- c(50, 150, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000)
# specify replication number
n <- 500
# specify alpha level
alpha <- 0.01
# allow parallel processing
plan(multisession) 
# set the seed
set.seed(123)



## =============================================================================
## 1. 5p - sparse condition
## =============================================================================

## Model specification
# specify B matrix
B5sparse = matrix(c(0, 0, 0, 0, 0,
                    1, 0, 0.8, 0, 0,
                    0, 0, 0, 0.9, 0,
                    0, 0.7, 0, 0, 1.5,
                    0, 0, 0, 0, 0), 5, 5, byrow = T)
dimnames(B5sparse) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
# specify layout
layout5 = matrix(c(0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1),5,2,byrow = T)

par(mfrow=c(1,2))
# true 5p sparse DCG
true5psparse <- qgraph(t(B5sparse), layout = layout5, 
                       labels = colnames(B5sparse), theme = "colorblind")

# true ancestral graph (AG)
dcg_5psparse <- matrix(c(0,1,0,0,0,
                         0,0,0,1,0,
                         0,1,0,0,0,
                         0,0,1,0,0,
                         0,0,0,1,0), 5, 5, byrow=T)
trueag_5psparse <- true_ancestral(dcg_5psparse, gen_dat(B5sparse), gaussCItest)
dimnames(trueag_5psparse) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
# plot true AG
plotAG(trueag_5psparse)


## Data generating
# equilibrium check
equilibrium_check(B5sparse)
# generate data (sample size as specified above)
simdata_5psparse <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B5sparse, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## Run CCD algorithm
# ccd_5psparse <- simdata_5psparse %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
#   )
# mat_5psparse <- ccd_5psparse %>%
#   map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_5psparse, file="    simulation/output/fixedB_n500/ccd_5psparse.RData")
# save(mat_5psparse, file="    simulation/output/fixedB_n500/mat_5psparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/mat_5psparse.RData") 

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

# save(fci_5psparse, file="        simulation/output/fixedB_n500/fci_5psparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/fci_5psparse.RData")

# plot resulting PAGs
# pag_fci5psparse <- fci_5psparse %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_5psparse <- simdata_5psparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=alpha, 
#                     labels = colnames(.x), p = ncol(.x)) %>% .$maag  
#   )
# save(cci_5psparse, file="fixedB_n500/    simulation/output/fixedB_n500/cci_5psparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/cci_5psparse.RData")

# plot resulting PAGs
# pag_cci5psparse <- cci_5psparse %>%
#   map_depth(2, ~plotAG(.x))



## =============================================================================
## 2. 5p - dense condition
## =============================================================================

## Model specification
# specify B matrix
B5dense = matrix(c(0, 0, 0, 0, 0,
                   1, 0, 0.8, 0, 0,
                   0, 0, 0, 0.9, 0,
                   0, 0.7, 0, 0, 1.5,
                   1, 0, 0, 0, 0), 5, 5, byrow = T)
dimnames(B5dense) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
par(mfrow=c(1,2))
# true 5p dense DCG
true5pdense <- qgraph(t(B5dense), layout=layout5, 
                      labels = colnames(B5dense), theme="colorblind")
# true ancestral graph (AG)
dcg_5pdense <- matrix(c(0,1,0,0,1,
                        0,0,0,1,0,
                        0,1,0,0,0,
                        0,0,1,0,0,
                        0,0,0,1,0), 5,5,byrow=T)
trueag_5pdense <- true_ancestral(dcg_5pdense, gen_dat(B5dense), gaussCItest)
dimnames(trueag_5pdense) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
# plot true AG
plotAG(trueag_5pdense)


## Data generating
# equilibrium check
equilibrium_check(B5dense)
# generate data (sample size as specified above)
simdata_5pdense <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B5dense, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## Run CCD algorithm
# ccd_5pdense <- simdata_5pdense %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
#   )
# mat_5pdense <- ccd_5pdense %>% 
#   map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_5pdense, file="    simulation/output/fixedB_n500/ccd_5pdense.RData")
# save(mat_5pdense, file="    simulation/output/fixedB_n500/mat_5pdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/mat_5pdense.RData")

# pag_ccd5pdense <- map2(ccd_5pdense, mat_5pdense,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_5pdense <- simdata_5pdense %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                     labels = colnames(.x)) %>% .@amat # extract amat
#   )
# save(fci_5pdense, file="    simulation/output/fixedB_n500/fci_5pdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/fci_5pdense.RData")

# plot resulting PAGs
# pag_fci5pdense <- fci_5pdense %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_5pdense <- simdata_5pdense %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, 
#                     alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag
#   )
# save(cci_5pdense, file="    simulation/output/fixedB_n500/cci_5pdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/cci_5pdense.RData")

# plot resulting PAGs
# pag_cci5pdense <- cci_5pdense %>%
#   map_depth(2, ~plotAG(.x))



## =============================================================================
## 3. 10p - sparse condition
## =============================================================================

## Model specification
# specify B matrix
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
# true 10p sparse DCG
true10psparse <- qgraph(t(B10sparse), layout = layout10, 
                        labels = colnames(B10sparse), theme="colorblind")

# true ancestral graph (AG)
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
dimnames(trueag_10psparse) <- list(paste("X", 1:10, sep=""), 
                                   paste("X", 1:10, sep=""))
# plot true AG
plotAG(trueag_10psparse)


## Data generating
# equilibrium check
equilibrium_check(B10sparse)
# generate data (sample size as specified above)
simdata_10psparse <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B10sparse, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## Run CCD algorithm
# ccd_10psparse <- simdata_10psparse %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
#   )
# mat_10psparse <- ccd_10psparse %>% 
#   map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_10psparse, file="    simulation/output/fixedB_n500/ccd_10psparse.RData")
# save(mat_10psparse, file="    simulation/output/fixedB_n500/mat_10psparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/mat_10psparse.RData")

# plot resulting PAGs
# pag_ccd10psparse <- map2(ccd_10psparse, mat_10psparse,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_10psparse <- simdata_10psparse %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                     labels = colnames(.x)) %>% .@amat # exxtract amat
#   )
# save(fci_10psparse, file="    simulation/output/fixedB_n500/fci_10psparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/fci_10psparse.RData")

# plot resulting PAGs
# pag_fci10psparse <- fci_10psparse %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_10psparse <- simdata_10psparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, 
#                     alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag
#   )
# save(cci_10psparse, file="    simulation/output/fixedB_n500/cci_10psparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/cci_10psparse.RData")

# plot resulting PAGs
# pag_cci10psparse <- cci_10psparse %>%
#   map_depth(2, ~plotAG(.x))



## =============================================================================
## 4. 10p - dense condition
## =============================================================================

## Model specification
# specify B matrix
B10dense = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0.3, 0, 0.8, 0, 0, 0, 0, 0, 0, 0, 
                    0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 
                    0, 0.4, 0, 1, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0.8, 0, 0.5, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0.8, 1, 
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 0.4, 
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
# true 10p dense DCG
true10pdense <- qgraph(t(B10dense), layout = layout10, 
                       labels = colnames(B10dense), theme="colorblind")
# true ancestral graph (AG)
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
dimnames(trueag_10pdense) <- list(paste("X", 1:10, sep=""), 
                                  paste("X", 1:10, sep=""))
# plot AG
plotAG(trueag_10pdense)


## Data generating
# equilibrium check
equilibrium_check(B10dense)
# generate data (sample size as specified above)
simdata_10pdense <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B10dense, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## Run CCD algorithm
# ccd_10pdense <- simdata_10pdense  %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
#   )
# mat_10pdense  <- ccd_10pdense  %>% 
#   map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))
# save(ccd_10pdense, file="    simulation/output/fixedB_n500/ccd_10pdense.RData")
# save(mat_10pdense, file="    simulation/output/fixedB_n500/mat_10pdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/mat_10pdense.RData")

# plot resulting PAGs
# pag_ccd10pdense <- map2(ccd_10pdense, mat_10pdense,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_10pdense <- simdata_10pdense  %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                     labels = colnames(.x)) %>% .@amat # exxtract amat
#   )
# save(fci_10pdense, file="    simulation/output/fixedB_n500/fci_10pdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/fci_10pdense.RData")

# plot resulting PAGs
# pag_fci10pdense <- fci_10pdense  %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_10pdense  <- simdata_10pdense %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, 
#                     alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  
#   )
# save(cci_10pdense, file="    simulation/output/fixedB_n500/cci_10pdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/cci_10pdense.RData")

# plot resulting PAGs
# pag_cci10pdense  <- cci_10pdense  %>%
#   map_depth(2, ~plotAG(.x))



## =============================================================================
## 5. 5p with LV sparse condition
## =============================================================================

## Model specification
# specify B matrix
B5_lvsparse = matrix(c(0, 0, 0, 0, 0, 1,
                       0, 0, 0.8, 0, 0, 1,
                       0, 0, 0, 0.9, 0,0,
                       0, 0.7, 0, 0, 1.5,0,
                       0, 0, 0, 0, 0,0,
                       0,0,0,0,0,0), 6, 6, byrow = T)
dimnames(B5_lvsparse) <- list(c(paste("X", 1:5, sep=""), "L1"), 
                              c(paste("X", 1:5, sep=""), "L1"))
# specify layout
layout5_lv = matrix(c(0,1,
                      0,0,
                      1,-1,
                      2,0,
                      2,1,
                      -1, 0.5),6,2,byrow = T)
par(mfrow=c(1,2))
# true 5p sparse w/LV DCG
true5p_lvsparse <- qgraph(t(B5_lvsparse), layout=layout5_lv, 
                          labels = colnames(B5_lvsparse), theme="colorblind")
# true ancestral graph (AG)
# [i,j] = [j,i] = 2: a LV exists between i and j
dcg_5psparseLV <- matrix(c(0, 2, 0, 0, 0, 
                           2, 0, 0, 1, 0, 
                           0, 1, 0, 0, 0,
                           0, 0, 1, 0, 0,
                           0, 0, 0, 1, 0), 5, 5, byrow = T)
trueag_5psparseLV <- true_ancestral(dcg_5psparseLV, gen_dat(B5_lvsparse), gaussCItest)
dimnames(trueag_5psparseLV) <- list(paste("X", 1:5, sep=""), 
                                    paste("X", 1:5, sep=""))
# plot AG
plotAG(trueag_5psparseLV)


## Data generating
# equilibrium check
equilibrium_check(B5_lvsparse)
# generate data (sample size as specified above)
simdata_5pLVsparse <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B5_lvsparse, N = z)[,-6],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## Run CCD algorithm
# ccd_5pLVsparse <- simdata_5pLVsparse  %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
#   )
# mat_5pLVsparse  <- ccd_5pLVsparse  %>% 
#   map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))
# save(ccd_5pLVsparse, file="    simulation/output/fixedB_n500/ccd_5pLVsparse.RData")
# save(mat_5pLVsparse, file="    simulation/output/fixedB_n500/mat_5pLVsparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/mat_5pLVsparse.RData")

# plot resulting PAGs
# pag_ccd5pLVsparse <- map2(ccd_5pLVsparse, mat_5pLVsparse,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_5pLVsparse <- simdata_5pLVsparse  %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                     labels = colnames(.x)) %>% .@amat # extract amat
#   )
# save(fci_5pLVsparse, file="    simulation/output/fixedB_n500/fci_5pLVsparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/fci_5pLVsparse.RData")

# plot resulting PAGs
# pag_fci_5pLVsparse <- fci_5pLVsparse  %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_5pLVsparse  <- simdata_5pLVsparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, 
#                     alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  
#   )
# 
# save(cci_5pLVsparse, file="    simulation/output/fixedB_n500/cci_5pLVsparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/cci_5pLVsparse.RData")

# plot resulting PAGs
# pag_cci_5pLVsparse <- cci_5pLVsparse  %>%
#   map_depth(2, ~plotAG(.x))



## =============================================================================
## 6. 5p with LV dense condition
## =============================================================================

## Model specification
# specify B matrix
B5_lvdense = matrix(c(0, 0, 0, 0, 0, 1,
                      0, 0, 0.8, 0, 0, 1,
                      0, 0, 0, 0.9, 0, 0,
                      0, 0.7, 0, 0, 1.5, 0,
                      1, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0), 6, 6, byrow = T)
dimnames(B5_lvsparse) <- list(c(paste("X", 1:5, sep=""), "L1"), 
                              c(paste("X", 1:5, sep=""), "L1"))
# specify layout
layout5_lv = matrix(c(0,1,
                      0,0,
                      1,-1,
                      2,0,
                      2,1,
                      -1, 0.5),6,2,byrow = T)
par(mfrow=c(1,2))
# true 5p dense w/LV DCG
true5p_lvdense <- qgraph(t(B5_lvdense), layout=layout5_lv, 
                         labels = colnames(B5_lvdense), theme="colorblind")
# true ancestral graph (AG)
# [i,j] = [j,i] = 2: a LV exists between i and j
dcg_5pdenseLV <- matrix(c(0, 2, 0, 0, 1, 
                          2, 0, 0, 1, 0, 
                          0, 1, 0, 0, 0,
                          0, 0, 1, 0, 0,
                          0, 0, 0, 1, 0), 5, 5, byrow = T)
trueag_5pdenseLV <- true_ancestral(dcg_5pdenseLV, gen_dat(B5_lvdense), gaussCItest)
dimnames(trueag_5pdenseLV) <- list(paste("X", 1:5, sep=""), 
                                   paste("X", 1:5, sep=""))
# plot AG
plotAG(trueag_5pdenseLV)


## Data generating
# equilibrium check
equilibrium_check(B5_lvdense)
# generate data (sample size as specified above)
simdata_5pLVdense <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B5_lvdense, N = z)[,-6],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## Run CCD algorithm
# ccd_5pLVdense <- simdata_5pLVdense  %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
#   )
# mat_5pLVdense  <- ccd_5pLVdense  %>% 
#   map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))
# save(ccd_5pLVdense, file="    simulation/output/fixedB_n500/ccd_5pLVdense.RData")
# save(mat_5pLVdense, file="    simulation/output/fixedB_n500/mat_5pLVdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/mat_5pLVdense.RData")

# plot resulting PAGs
# pag_ccd5pLVdense <- map2(ccd_5pLVdense, mat_5pLVdense,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_5pLVdense <- simdata_5pLVdense  %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                     labels = colnames(.x)) %>% .@amat
#   )

# save(fci_5pLVdense, file="    simulation/output/fixedB_n500/fci_5pLVdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/fci_5pLVdense.RData")

# plot resulting PAGs
# pag_fci5pLVdense <- fci_5pLVdense  %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_5pLVdense  <- simdata_5pLVdense %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, 
#                     alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag 
#   )

# save(cci_5pLVdense, file="    simulation/output/fixedB_n500/cci_5pLVdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/cci_5pLVdense.RData")

# plot resulting PAGs
# pag_cci5pLVdense <- cci_5pLVdense  %>%
#   map_depth(2, ~plotAG(.x))



## =============================================================================
## 7. 10p with LV sparse condition
## =============================================================================

## Model specification
# specify B matrix
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
dimnames(B10_lvsparse) <- list(c(paste("X", 1:10, sep=""), "L1", "L2"), 
                               c(paste("X", 1:10, sep=""), "L1", "L2"))
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
# true 10p sparse w/LV DCG
true10pLVsparse <- qgraph(t(B10_lvsparse), layout = layout10LV2, 
                          labels = colnames(B10_lvsparse), theme="colorblind")
# true ancestral graph (AG)
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
dimnames(trueag_10psparseLV) <- list(paste("X", 1:10, sep=""), 
                                     paste("X", 1:10, sep=""))
# plot AG
plotAG(trueag_10psparseLV)


## Data generating
# equilibrium check
equilibrium_check(B10_lvsparse)
# generate data (sample size as specified above)
simdata_10pLVsparse <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B10_lvsparse, N = z)[,-c(11,12)],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## Run CCD algorithm
# ccd_10pLVsparse  <- simdata_10pLVsparse   %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
#   )
# mat_10pLVsparse   <- ccd_10pLVsparse %>% 
#   map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))
# save(ccd_10pLVsparse, file="    simulation/output/fixedB_n500/ccd_10pLVsparse.RData")
# save(mat_10pLVsparse, file="    simulation/output/fixedB_n500/mat_10pLVsparse.RData")
# load the saved result
load("    simulation/output/fixedB_n500/mat_10pLVsparse.RData")

# plot resulting PAGs
# pag_ccd10pLVsparse <- map2(ccd_10pLVsparse, mat_10pLVsparse,
#                         ~map2(.x, .y, plotPAG)
#                          )

## Run FCI algorithm
# fci_10pLVsparse  <- simdata_10pLVsparse   %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                     labels = colnames(.x)) %>% .@amat  
#   )
# save(fci_10pLVsparse, file="    simulation/output/fixedB_n500/fci_10pLVsparse.RData")
load("    simulation/output/fixedB_n500/fci_10pLVsparse.RData")

# plot resulting PAGs
# pag_fci10pLV  <- fci_10pLVsparse %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_10pLVsparse  <- simdata_10pLVsparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, 
#                     alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  
#   )
# save(cci_10pLVsparse, file="    simulation/output/fixedB_n500/cci_10pLVsparse.RData")
load("    simulation/output/fixedB_n500/cci_10pLVsparse.RData")

# plot resulting PAGs
# pag_cci10pLVsparse   <- cci_10pLVsparse %>%
#   map_depth(2, ~plotAG(.x))



## =============================================================================
## 8. 10p with LV dense condition
## =============================================================================

## Model specification
# specify B matrix
B10_lvdense = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0.5, 0, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0.7,
                       0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 0, 0,
                       0, 0, 1.2, 1, 0, 0, 0, 0, 0, 0, 0, 0.7,
                       0, 0, 0, 0, 0.8, 0, 0.5, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0.8, 0.5, 0.4, 0,
                       0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0.8, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 12, 12, byrow = T)
dimnames(B10_lvsparse) <- list(c(paste("X", 1:10, sep=""), "L1", "L2"), 
                               c(paste("X", 1:10, sep=""), "L1", "L2"))
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
# true 10p sparse w/LV DCG
true10pLVdense <- qgraph(t(B10_lvdense), layout = layout10LV2, 
                         labels = colnames(B10_lvdense), theme="colorblind")
# true ancestral graph (AG)
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
dimnames(trueag_10pdenseLV) <- list(paste("X", 1:10, sep=""), 
                                    paste("X", 1:10, sep=""))
# plot AG
plotAG(trueag_10pdenseLV)


## Data generating
# equilibrium check
equilibrium_check(B10_lvdense)
# generate data (sample size as specified above)
simdata_10pLVdense <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(B10_lvdense, N = z)[,-c(11,12)],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## Run CCD algorithm
# ccd_10pLVdense  <- simdata_10pLVdense   %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
#   )
# mat_10pLVdense   <- ccd_10pLVdense %>% 
#   map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_10pLVdense, file="    simulation/output/fixedB_n500/ccd_10pLVdense.RData")
# save(mat_10pLVdense, file="    simulation/output/fixedB_n500/mat_10pLVdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/mat_10pLVdense.RData")

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
# save(fci_10pLVdense, file="    simulation/output/fixedB_n500/fci_10pLVdense.RData")
# load the saved result
load("    simulation/output/fixedB_n500/fci_10pLVdense.RData")

# plot resulting PAGs
# pag_fci10pLVdense <- fci_10pLVdense   %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_10pLVdense  <- simdata_10pLVdense  %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, 
#                     alpha = alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag 
#   )
# save(cci_10pLVdense, file="    simulation/output/fixedB_n500/cci_10pLVdense.RData")
load("    simulation/output/fixedB_n500/cci_10pLVdense.RData")

# plot resulting PAGs
# pag_cci10pLVdense   <- cci_10pLVdense %>%
#   map_depth(2, ~plotAG(.x))



## =============================================================================
## 9. Algorithm running time (*extra*)
## =============================================================================

## load the package 
# remotes::install_github("joshuaulrich/microbenchmark")
library(microbenchmark)

## specify n, alpha, and dataset to use
n <- 1e3
alpha <- 0.05
dataset <- list(simdata_5psparse[[1]][[1]], simdata_5pdense[[1]][[1]],
                simdata_10psparse[[1]][[1]], simdata_10pdense[[1]][[1]],
                simdata_5pLVsparse[[1]][[1]], simdata_5pLVdense[[1]][[1]],
                simdata_10pLVsparse[[1]][[1]],simdata_10pLVsparse[[1]][[1]])


## compute the algorithm running time 
# times <- microbenchmark(
#   ccd_5psparse = ccdKP(df=dataset[[1]], dataType = "continuous", alpha = alpha),
#   fci_5psparse = fci(list(C = cor(dataset[[1]]), n = n),indepTest = gaussCItest, 
#                      alpha = alpha, selectionBias= FALSE, labels = colnames(dataset[[1]])),
#   cci_5psparse = cci(list(C = cor(dataset[[1]]), n = n), gaussCItest, 
#                      alpha = alpha, p=ncol(dataset[[1]])),
# 
#   ccd_5pdense = ccdKP(df=dataset[[2]], dataType = "continuous", alpha = alpha),
#   fci_5pdense = fci(list(C = cor(dataset[[2]]), n = n),indepTest = gaussCItest, 
#                     alpha = alpha, selectionBias= FALSE, labels = colnames(dataset[[2]])),
#   cci_5pdense = cci(list(C = cor(dataset[[2]]), n = n), gaussCItest, 
#                     alpha = alpha, p=ncol(dataset[[2]])),
# 
#   ccd_10psparse = ccdKP(df=dataset[[3]], dataType = "continuous", alpha = alpha),
#   fci_10psparse = fci(list(C = cor(dataset[[3]]), n = n),indepTest = gaussCItest, 
#                       alpha = alpha, selectionBias= FALSE, labels = colnames(dataset[[3]])),
#   cci_10psparse = cci(list(C = cor(dataset[[3]]), n = n), gaussCItest, 
#                       alpha = alpha, p=ncol(dataset[[3]])),
# 
#   ccd_10pdense = ccdKP(df=dataset[[4]], dataType = "continuous", alpha = alpha),
#   fci_10pdense = fci(list(C = cor(dataset[[4]]), n = 1e3),indepTest=gaussCItest, 
#                      alpha = alpha, selectionBias= FALSE, labels = colnames(dataset[[4]])),
#   cci_10pdense = cci(list(C = cor(dataset[[4]]), n = 1e3), gaussCItest, 
#                      alpha = alpha, p=ncol(dataset[[4]])),
# 
#   ccd_5pLVsparse = ccdKP(df=dataset[[5]], dataType = "continuous", alpha = alpha),
#   fci_5pLVsparse = fci(list(C = cor(dataset[[5]]), n = n),indepTest=gaussCItest,
#                 alpha = alpha, selectionBias= FALSE, labels = colnames(dataset[[5]])),
#   cci_5pLVsparse = cci(list(C = cor(dataset[[5]]), n = n), gaussCItest, 
#                        alpha = alpha, p=ncol(dataset[[5]])),
# 
#   ccd_5pLVdense = ccdKP(df=dataset[[6]], dataType = "continuous", alpha = alpha),
#   fci_5pLVdense = fci(list(C = cor(dataset[[6]]), n = n),indepTest=gaussCItest,
#                        alpha = alpha, selectionBias= FALSE, labels = colnames(dataset[[6]])),
#   cci_5pLVdense = cci(list(C = cor(dataset[[6]]), n = n), gaussCItest, 
#                       alpha=alpha, p=ncol(dataset[[6]])),
#   
#   ccd_10pLVsparse = ccdKP(df=dataset[[7]], dataType = "continuous", alpha = alpha),
#   fci_10pLVsparse = fci(list(C = cor(dataset[[7]]), n = n),indepTest=gaussCItest,
#                 alpha = alpha, selectionBias= FALSE, labels = colnames(dataset[[7]])),
#   cci_10pLVsparse = cci(list(C = cor(dataset[[7]]), n = n), gaussCItest, 
#                         alpha = alpha, p=ncol(dataset[[7]])),
#   
#   ccd_10pLVdense = ccdKP(df=dataset[[8]], dataType = "continuous", alpha = alpha),
#   fci_10pLVdense = fci(list(C = cor(dataset[[8]]), n = n),indepTest=gaussCItest,
#                         alpha = alpha, selectionBias= FALSE, labels = colnames(dataset[[8]])),
#   cci_10pLVdense = cci(list(C = cor(dataset[[8]]), n = n), gaussCItest, 
#                        alpha = alpha, p=ncol(dataset[[8]]))
# )

# save(times, file="    simulation/output/fixedB_n500/algo_runningtime.RData")
# load the saved result
load("    simulation/output/fixedB_n500/algo_runningtime.RData")

## neat up result object
times <- times %>%
  mutate(algorithm = toupper(substr(expr, 1, 3)),
         condition = stringr::str_split(expr, "_", simplify=T)[,2])


## create the figure of algorithm running time
timeplot <- times %>%
  # create ggplot object
  ggplot(aes(x=factor(condition, levels= c("5psparse", "5pdense", "10psparse", "10pdense", 
                                           "5pLVsparse","5pLVdense", "10pLVsparse","10pLVdense")), 
             y = log(time), col= factor(algorithm))) +
  # add boxplots
  geom_boxplot(position = "dodge", outlier.size = 0.8, outlier.alpha = 0.2) + 
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  # specify plot labels
  labs(y = " log(ms)", x = "conditions", title = "Algorithm Running Time", 
       subtitle = "Time in milliseconds (ms)") +
  scale_x_discrete(labels=c("5p sparse", "5p dense", "10p sparse", "10p dense", 
                            "5p sparse LC","5p dense LC", "10p sparse LC","10p dense LC")) +
  # apply the themes
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", margin = margin(t = 13)),
        plot.title = element_text(face = "bold", family = "Palatino", size = 15),
        plot.subtitle = element_text(face = "italic", family = "Palatino", size = 15),
        axis.text=element_text(face = "bold",family = "Palatino", size = 11),
        axis.title = element_text(face = "bold",family = "Palatino", size = 12),
        legend.text = element_text(face = "bold", family = "Palatino", size = 12),
        legend.position="bottom") 

## save the plot
# ggsave(timeplot, filename = "figures/algo_time.pdf", width = 25, height = 13, dpi = 300, units = "cm")
