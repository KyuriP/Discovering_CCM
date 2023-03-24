
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
source("R/CCD_fnc.R")
source("R/plot_fnc.R")

## import the example empirical data
mcnally <- read.csv("../data/McNally.csv")

# separate depression / OCD symptoms
# (original data contains both depression and OCD symptoms)
# (here we only use depression symptoms)
depression <- mcnally[,1:16]
ocd <- mcnally[,17:26]

# paranormal transformation using huge package
# trans_dep <- huge::huge.npn(depression)
# trans_ocd <- huge::huge.npn(ocd)


## =======================================
## 2. Estimate GGM with GLASSO 
## =======================================
## estimate GGM via graphical LASSO on depression symptoms
cordep <- cor(depression)
# found the optimal sparsity by gamma = 1
glassoFitdep <- EBICglasso(cordep, n = nrow(depression), gamma = 1)
qgraph(glassoFitdep, layout = "spring", theme="colorblind",
       nodeNames = colnames(depression), legend.cex = 0.4)


## =======================================
## 3. Estimate PAGs using CCD, FCI and CCI
## =======================================
set.seed(123)
## estimate the PAG on depression symptoms by running CCD
# run CCD
ccd_mcnally_dep <- ccdKP(df=depression, dataType = "discrete", alpha=0.05)
# create an adjacency matrix for PAG
mat_mcnally_dep <- CreateAdjMat(ccd_mcnally_dep, p = ncol(depression))
# plot the PAG
pag_mcnally_dep <- plotPAG(ccd_mcnally_dep, mat_mcnally_dep)

## estimate the PAG on depression symptoms by running FCI
fci(list(C = cor(depression), n = nrow(depression)), gaussCItest, alpha=0.05, 
    labels = colnames(depression), verbose=TRUE) %>% .@amat %>% plotAG 

## estimate the PAG on depression symptoms by running CCI
cci(list(C = cor(depression), n = nrow(depression)), gaussCItest, alpha=0.05, 
    labels = colnames(depression), p = ncol(depression), verbose=TRUE) %>% .$maag %>% plotAG

1/sqrt(408)
