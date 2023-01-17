
## =============================================================================
## Description
#
# This script contains simulated models that are used in the original study.
# There are in total 6 models considered: 4nodes-sparse, 4nodes-dense, 5nodes-sparse,
# 5nodes-dense, 6nodes-sparse, 6nodes-dense.
#
# We estimate a GGM and PAG (using the CCD algorithm) for each model and also 
# retrieve a Markov equivalence class of DCGs based on the estimated PAG from each condition.
#
# They are commented out; simply un-comment them if one is interested in running the models.
# But note that computing an equivalence class can take several minutes.
## =============================================================================


## ========================
## Preparation
## ========================
## load necessary packages
library(qgraph)
library(pcalg)
library(dplyr)

## source all the necessary functions
source("code/R/CCD_fnc.R")
source("code/R/plot_fnc.R")
source("code/R/dsep_fnc.R")
source("code/R/searchAM_KP_fnc.R")
source("code/R/equivset_fnc.R")
source("code/R/data_generating_fnc.R")
source("code/R/eval_metric_fnc.R")


## ========================
## Model 1) 4 nodes - sparse
## ========================

# specify B matrix
p = 4
B4 = matrix(c(0, 0, 0, 0,
              1, 0, 0.5, 0,
              0, 0.5, 0, 0.9,
              0, 0, 0, 0), p, p, byrow = T)

colnames(B4) <- c("X1", "X2", "X3", "X4")
# specify layout
layout4 = matrix(c(-1,1,
                   -1,0,
                   1,0,
                   1,1),4,2,byrow = T)
# ## True graph
# true4p <- qgraph(t(B4), layout=layout4, labels = colnames(B4), theme="colorblind")
# 
# ## Data generating
# # equilibrium check
# equilibrium_check(B4)
# # generate data
# data4p <- gen_dat(B4, N =1e6, seed = 1)
# 
# ## Estimate GGM
# ggm4p <- qgraph(cor(data4p), layout=layout4, theme="colorblind")
# 
# ## Run CCD algorithm
# ccd_4p <- ccdKP(df=data4p, dataType = "continuous", alpha = 0.05)
# mat4p <- CreateAdjMat(ccd_4p, 4)
# 
# ## Estimate PAG
# pag4p <- plotPAG(ccd_4p, mat4p)

# ## Compute equivalence class of all DCGs given the PAG
# # (this takes relatively a long time, so we save the object)
# equiv4p <- semiequiv_dcg(ccd_4p, mat4p)
# save(equiv4p, file="data/equiv4p.RData")
# load("data/equiv4p.RData")


## ========================
## Model 2) 4 nodes - dense
## ========================
# specify B matrix
p = 4
B4_high = matrix(c(0, 0, 0, 0,
                   0.9, 0, 0.4, 0,
                   0, 0.5, 0, .5,
                   -0.8, 0, 0, 0), p, p, byrow = T)
colnames(B4_high) <- c("X1", "X2", "X3", "X4")

# ## True graph
# true4p_high <- qgraph(t(B4_high), layout=layout4, labels = colnames(B4_high), theme="colorblind")
# 
# ## Data generating
# # equilibrium check
# equilibrium_check(B4_high)
# # generate data
# data4p_high <- gen_dat(B4_high, N =1e6, seed = 1)
# 
# ## Estimate GGM
# ggm4p_high <- qgraph(t(cor(data4p_high)), layout=layout4, theme="colorblind")
# 
# ## Run CCD algorithm
# ccd_4p_high <- ccdKP(df=data4p_high, dataType = "continuous", alpha = 0.05)
# mat4p_high <- CreateAdjMat(ccd_4p_high, 4)
# 
# ## Estimate PAG
# pag4p <- plotPAG(ccd_4p_high, mat4p_high)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv4p_high <- semiequiv_dcg(ccd_4p_high, mat4p_high)
# save(equiv4p_high, file="data/equiv4p_high.RData")
# load("data/equiv4p_high.RData")


## ========================
## Model 3) 5 nodes - sparse
## ========================
# specify B matrix
p = 5

B5 = matrix(c(0, 1, 0, 0, 0,
              0, 0, 0, 0.7, 0,
              0, 0.4, 0, 0, 0,
              0, 0, .5, 0, 0,
              0, 0, 0, -1.5, 0), p, p, byrow = T)
colnames(B5) <- c("X1", "X2", "X3", "X4", "X5")
# specify layout
layout5 = matrix(c(0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1),5,2,byrow = T)
# ## True graph
# true5p <- qgraph(t(B5), layout=layout5, labels = colnames(B5), theme="colorblind")
# 
# ## Data generating
# # equilibrium check
# equilibrium_check(B5)
# # generate data
# data5p <- gen_dat(B5, N =1e6, seed = 123)
# 
# ## Estimate GGM
# ggm5p <- qgraph(cor(data5p), layout = layout5, theme="colorblind")
# 
# ## Run CCD algorithm
# ccd_5p <- ccdKP(df=data5p, dataType = "continuous", alpha = 0.05)
# mat5p <- CreateAdjMat(ccd_5p, 5)
# 
# ## Estimate PAG
# pag5p <- plotPAG(ccd_5p, mat5p)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv5p <- semiequiv_dcg(ccd_5p, mat5p)
# save(equiv5p, file="data/equiv5p.RData")
# load("data/equiv5p.RData")


## ========================
## Model 4) 5 nodes - dense
## ========================
# specify B matrix
p = 5

B5_high = matrix(c(0, 0.9, 0, 0, 0.6,
              0, 0, 0, 0.7, 0,
              0, 0.9, 0, 0, 0,
              0, 0, 0.5, 0, 0,
              0, 0, 0, 1, 0), p, p, byrow = T)
colnames(B5_high) <- c("X1", "X2", "X3", "X4", "X5")
# specify layout
layout5 = matrix(c(0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1),5,2,byrow = T)

# ## True graph
# true5p_high <- qgraph(t(B5_high), layout=layout5, labels = colnames(B5_high), theme="colorblind")
# 
# ## Data generating
# # equilibrium check
# equilibrium_check(B5_high)
# # generate data
# data5p_high <- gen_dat(B5_high, N =1e6, seed = 1)
# 
# ## Estimate GGM
# ggm5p_high <- qgraph(cor(data5p_high), layout = layout5, theme="colorblind")
# 
# ## Run CCD algorithm
# ccd_5p_high <- ccdKP(df=data5p_high, dataType = "continuous", alpha = 0.05)
# mat5p_high <- CreateAdjMat(ccd_5p_high, 5)
# 
# ## Estimate PAG
# pag5p_high <- plotPAG(ccd_5p_high, mat5p_high)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv5p_high <- semiequiv_dcg(ccd_5p_high, mat5p_high)
# save(equiv5p_high, file="data/equiv5p_high.RData")
# load("data/equiv5p_high.RData")


## ========================
## Model 5) 6 nodes - sparse
## ========================
# specify B matrix
p = 6
B6 = matrix(c(0, 0, 0, 0, 0, 0,
              0.3, 0, 0.4, 0, 0, 0,
              0, 0, 0, 0.9, 0, 0,
              0, 0, 0, 0, 0.4, 0,
              0, 0, 1, 0, 0, 0,
              1, 0, 0, 0, 0.5, 0), p, p, byrow = T)
colnames(B6) <- c("X1", "X2", "X3", "X4", "X5", "X6")
# specify layout
layout6 = matrix(c(1, 2,
                   0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1),6,2,byrow = T)

# ## True graph
# true6p <- qgraph(t(B6), layout=layout6, labels = colnames(B6), theme="colorblind")
# 
# ## Data generating
# # equilibrium check
# equilibrium_check(B6)
# # generate data
# data6p <- gen_dat(B6, N =1e6, seed = 123)
# 
# ## GGM
# ggm6p <- qgraph(cor(data6p), layout = layout6, theme="colorblind")
# 
# 
# ## Run CCD algorithm
# ccd_6p <- ccdKP(df=data6p, dataType = "continuous", alpha = 0.05)
# mat6p <- CreateAdjMat(ccd_6p, 6)
# 
# ## Estimate PAG
# pag6p <- plotPAG(ccd_6p, mat6p)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv6p <- semiequiv_dcg(ccd_6p, mat6p)
# save(equiv6p, file="data/equiv6p.RData")
# load("data/equiv6p.RData")


## ========================
## Model 6) 6 nodes - dense
## ========================
# specify B matrix
p = 6
B6_high = matrix(c(0, 0, 0, 0, 0, 0,
              0.7, 0, 0.4, 0, 0, 0.9,
              0, 0, 0, 0.9, 0, 0,
              0, 0, 0, 0, 0.4, 0,
              0, 0, 1, 0, 0, 0,
              1, 0, 0, 0, 0.5, 0), p, p, byrow = T)
# colnames for B matrix is necessary for running CCD
colnames(B6_high) <- c("X1", "X2", "X3", "X4", "X5", "X6")


# true6p_high <- qgraph(t(B6_high), layout=layout6, labels = colnames(B6), theme="colorblind")
# 
# ## Data generating
# # equilibrium check
# equilibrium_check(B6_high)
# # generate data
# data6p_high<- gen_dat(B6_high, N =1e6, seed = 123)
# 
# ## Estimate GGM
# ggm6p_high <- qgraph(cor(data6p_high), layout = layout6, theme="colorblind")
# 
# ## Run CCD algorithm
# ccd_6p_high <- ccdKP(df=data6p_high, dataType = "continuous", alpha = 0.05)
# mat6p_high <- CreateAdjMat(ccd_6p_high, 6)
# 
# ## Estimate PAG
# pag6p_high <- plotPAG(ccd_6p_high, mat6p_high)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv6p_high <- semiequiv_dcg(ccd_6p_high, mat6p_high)
# save(equiv6p_high, file="data/equiv6p_high.RData")
# load("data/equiv6p_high.RData")

