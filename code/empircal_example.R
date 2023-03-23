
## =============================================================================
## Description
#
# This script concerns applying the CCD, FCI and CCI algorithms to an empirical data 
# from McNally et al. (2017) in order to test the practical applicability of the algorithms.
#
# The data is obtained from : https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900#supplementary-materials

# It contains the code to create "Figure 18. Estimated statistical network model and PAGs"
## =============================================================================



## =======================================
## 1. Preparation
## =======================================
## load necessary packages
library(qgraph)
library(pcalg)
library(ggplot2)
library(dplyr)

## source all the necessary functions
source("code/R/CCD_fnc.R")
source("code/R/plot_fnc.R")

## import the example empirical data
mcnally <- read.csv("data/McNally.csv")

# separate depression / OCD symptoms
# (original data contains both depression and OCD symptoms)
# (here we only use depression symptoms)
depression <- mcnally[,1:16]
ocd <- mcnally[,17:26]

# paranormal transformation using huge package
trans_dep <- huge::huge.npn(depression)
trans_ocd <- huge::huge.npn(ocd)


## =======================================
## 2. Estimate GGM with GLASSO 
## =======================================
## estimate GGM via graphical LASSO on depression symptoms
cordep <- cor(trans_dep)
# found the optimal sparsity by gamma = 1
glassoFitdep <- EBICglasso(cordep, n = nrow(trans_dep), gamma = 1)
qgraph(glassoFitdep, layout = "spring", theme="colorblind",
       nodeNames = colnames(depression), legend.cex = 0.4)


## =======================================
## 3. Estimate PAGs using CCD, FCI and CCI
## =======================================
## estimate the PAG on depression symptoms by running CCD
# run CCD
ccd_mcnally_dep <- ccdKP(df=trans_dep, dataType = "continuous", depth = -1, alpha=0.01)
# create an adjacency matrix for PAG
mat_mcnally_dep <- CreateAdjMat(ccd_mcnally_dep, p = ncol(trans_dep))
# plot the PAG
pag_mcnally_dep <- plotPAG(ccd_mcnally_dep, mat_mcnally_dep)

## estimate the PAG on depression symptoms by running FCI
fci(list(C = cor(trans_dep), n = nrow(trans_dep)), gaussCItest, alpha=0.01, 
    labels = colnames(trans_dep), verbose=TRUE) %>% .@amat %>% plotAG 

## estimate the PAG on depression symptoms by running CCI
cci(list(C = cor(trans_dep), n = nrow(trans_dep)), gaussCItest, alpha=0.01, 
    labels = colnames(trans_dep), p = ncol(trans_dep), verbose=TRUE) %>% .$maag %>% plotAG

