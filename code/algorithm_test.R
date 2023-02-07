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

## slightly modified CCI package
#install_github("KyuriP/CCI_KP")
library(CCI.KP)

# for reproducibility
#set.seed(12345)

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
## True graph
true4p <- qgraph(t(B4), layout=layout4, labels = colnames(B4), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B4)
# generate data
data4p <- gen_dat(B4, N =1e6, seed = 12345)
## Estimate GGM
ggm4p <- qgraph(cor(data4p), layout=layout4, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_4p <- ccdKP(df=data4p, dataType = "continuous", alpha = 0.05)
mat4p <- CreateAdjMat(ccd_4p, 4)
# Estimate PAG
pag_ccd4p <- plotPAG(ccd_4p, mat4p)

# ## Compute equivalence class of all DCGs given the PAG
# # (this takes relatively a long time, so we save the object)
# equiv4p <- semiequiv_dcg(ccd_4p, mat4p)
# save(equiv4p, file="data/equiv4p.RData")
# load("data/equiv4p.RData")

## Run FCI algorithm
suffStat_4p = list()
suffStat_4p$C = cor(data4p)
suffStat_4p$n = 1e6

res4p <- fci(suffStat_4p,indepTest=gaussCItest,
             alpha = 0.05, labels = colnames(data4p), selectionBias = FALSE, rules = rep(TRUE, 10))

pag_fci4p <- plotAG(res4p@amat)

## Run CCI algorithm
G4p = cci(suffStat_4p, gaussCItest, alpha=0.05, p=ncol(data4p)) 
dimnames(G4p$maag) <- list(colnames(data4p), colnames(data4p)) # give the labels
pag_cci4p <- plotAG(G4p$maag)

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

## True graph
true4p_high <- qgraph(t(B4_high), layout=layout4, labels = colnames(B4_high), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B4_high)
# generate data
data4p_high <- gen_dat(B4_high, N =1e6, seed = 12345)

## Estimate GGM
ggm4p_high <- qgraph(t(cor(data4p_high)), layout=layout4, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_4p_high <- ccdKP(df=data4p_high, dataType = "continuous", alpha = 0.05)
mat4p_high <- CreateAdjMat(ccd_4p_high, 4)
## Estimate PAG
pag_ccd4pH <- plotPAG(ccd_4p_high, mat4p_high)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv4p_high <- semiequiv_dcg(ccd_4p_high, mat4p_high)
# save(equiv4p_high, file="data/equiv4p_high.RData")
# load("data/equiv4p_high.RData")

## Run FCI algorithm
suffStat_4phigh = list()
suffStat_4phigh$C = cor(data4p_high)
suffStat_4phigh$n = 1e6

res4p_high <- fci(suffStat_4phigh,indepTest=gaussCItest,
             alpha = 0.05, doPdsep = FALSE, labels = colnames(data4p_high))

pag_fci4pH <- plotAG(res4p_high@amat)

## Run CCI algorithm
G4p_high = cci(suffStat_4phigh, gaussCItest, alpha=0.05, p=ncol(data4p_high)) 
dimnames(G4p_high$maag) <- list(colnames(data4p_high), colnames(data4p_high)) # give the labels
pag_cci4pH <- plotAG(G4p_high$maag)


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
## True graph
true5p <- qgraph(t(B5), layout=layout5, labels = colnames(B5), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5)
# generate data
data5p <- gen_dat(B5, N =1e6, seed = 1)

## Estimate GGM
ggm5p <- qgraph(cor(data5p), layout = layout5, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_5p <- ccdKP(df=data5p, dataType = "continuous", alpha = 0.05)
mat5p <- CreateAdjMat(ccd_5p, 5)

## Estimate PAG
pag_ccd5p <- plotPAG(ccd_5p, mat5p)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv5p <- semiequiv_dcg(ccd_5p, mat5p)
# save(equiv5p, file="data/equiv5p.RData")
# load("data/equiv5p.RData")

## Run FCI algorithm
suffStat_5p = list()
suffStat_5p$C = cor(data5p)
suffStat_5p$n = 1e6

res5p <- fci(suffStat_5p,indepTest=gaussCItest,
             alpha = 0.05, doPdsep = FALSE, labels = colnames(data5p))

pag_fci5p <- plotAG(res5p@amat)

## Run CCI algorithm
G5p = cci(suffStat_5p, gaussCItest, alpha=0.05, p=ncol(data5p)) 
dimnames(G5p$maag) <- list(colnames(data5p), colnames(data5p)) # give the labels
pag_cci5p <- plotAG(G5p$maag)

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

## True graph
true5p_high <- qgraph(t(B5_high), layout=layout5, labels = colnames(B5_high), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5_high)
# generate data
data5p_high <- gen_dat(B5_high, N =1e6, seed = 1)

## Estimate GGM
ggm5p_high <- qgraph(cor(data5p_high), layout = layout5, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_5p_high <- ccdKP(df=data5p_high, dataType = "continuous", alpha = 0.05)
mat5p_high <- CreateAdjMat(ccd_5p_high, 5)

## Estimate PAG
pag_ccd5pH <- plotPAG(ccd_5p_high, mat5p_high)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv5p_high <- semiequiv_dcg(ccd_5p_high, mat5p_high)
# save(equiv5p_high, file="data/equiv5p_high.RData")
# load("data/equiv5p_high.RData")

## Run FCI algorithm
suffStat_5phigh = list()
suffStat_5phigh$C = cor(data5p_high)
suffStat_5phigh$n = 1e6

res5p_high <- fci(suffStat_5phigh,indepTest=gaussCItest,
             alpha = 0.05, doPdsep = FALSE, labels = colnames(data5p_high))

pag_fci5pH <- plotAG(res5p_high@amat)

## Run CCI algorithm
G5p_high = cci(suffStat_5phigh, gaussCItest, alpha=0.05, p=ncol(data5p_high)) 
dimnames(G5p_high$maag) <- list(colnames(data5p_high), colnames(data5p_high)) # give the labels
pag_cci5pH <- plotAG(G5p_high$maag)


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

## True graph
true6p <- qgraph(t(B6), layout=layout6, labels = colnames(B6), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B6)
# generate data
data6p <- gen_dat(B6, N =1e6, seed = 123)

## GGM
ggm6p <- qgraph(cor(data6p), layout = layout6, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_6p <- ccdKP(df=data6p, dataType = "continuous", alpha = 0.05)
mat6p <- CreateAdjMat(ccd_6p, 6)

## Estimate PAG
pag_ccd6p <- plotPAG(ccd_6p, mat6p)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv6p <- semiequiv_dcg(ccd_6p, mat6p)
# save(equiv6p, file="data/equiv6p.RData")
# load("data/equiv6p.RData")

## Run FCI algorithm
suffStat_6p = list()
suffStat_6p$C = cor(data6p)
suffStat_6p$n = 1e6

res6p <- fci(suffStat_6p,indepTest=gaussCItest,
                  alpha = 0.05, labels = colnames(data6p), selectionBias = FALSE)

pag_fci6p <- plotAG(res6p@amat)

## Run CCI algorithm
G6p = cci(suffStat_6p, gaussCItest, alpha=0.05, p=ncol(data6p)) 
dimnames(G6p$maag) <- list(colnames(data6p), colnames(data6p)) # give the labels
pag_cci6p <- plotAG(G6p$maag)


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


true6p_high <- qgraph(t(B6_high), layout=layout6, labels = colnames(B6), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B6_high)
# generate data
data6p_high<- gen_dat(B6_high, N =1e6, seed = 123)

## Estimate GGM
ggm6p_high <- qgraph(cor(data6p_high), layout = layout6, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_6p_high <- ccdKP(df=data6p_high, dataType = "continuous", alpha = 0.05)
mat6p_high <- CreateAdjMat(ccd_6p_high, 6)

## Estimate PAG
pag_ccd6pH <- plotPAG(ccd_6p_high, mat6p_high)

## Compute equivalence class of all DCGs given the PAG
# (this takes relatively a long time, so we save the object)
# equiv6p_high <- semiequiv_dcg(ccd_6p_high, mat6p_high)
# save(equiv6p_high, file="data/equiv6p_high.RData")
# load("data/equiv6p_high.RData")

## Run FCI algorithm
suffStat_6phigh = list()
suffStat_6phigh$C = cor(data6p_high)
suffStat_6phigh$n = 1e6

# selection bias = FALSE (Mooij --> absence of selection bias assumed to use FCI for cyclci structure)
res6p_high <- fci(suffStat_6phigh,indepTest=gaussCItest,
             alpha = 0.05, labels = colnames(data6p_high), selectionBias = FALSE, rules = rep(TRUE, 10))

pag_fci6pH <- plotAG(res6p_high@amat)

## Run CCI algorithm
G6p_high = cci(suffStat_6phigh, gaussCItest, alpha=0.05, p=ncol(data6p_high)) 
dimnames(G6p_high$maag) <- list(colnames(data6p_high), colnames(data6p_high)) # give the labels
pag_cci6pH <- plotAG(G6p_high$maag)



## ========================
## Model 7-1) 4 nodes with a LV - sparse
## ========================
# specify B matrix
B4_LV = matrix(c(0, 0, 0, 0, -0.4,
              0, 0, 0.6, 0, 0,
              0, 0.6, 0, 0.9, 0,
              0, 0, 0, 0, -1,
              0, 0, 0, 0, 0), 5, 5, byrow = T)

colnames(B4_LV) <- c("X1", "X2", "X3", "X4", "L1")
# specify layout
layout4LV = matrix(c(-1,1,
                   -1,0,
                   1,0,
                   1,1,
                   0, 1.5),5,2,byrow = T)
## True graph
true4pLV <- qgraph(t(B4_LV), layout=layout4LV, labels = colnames(B4_LV), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B4_LV)
# generate data
data4pLV <- gen_dat(B4_LV, N =1e6, seed = 123)[,-5]
## Estimate GGM
ggm4pLV <- qgraph(cor(data4pLV), layout=layout4, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_4pLV <- ccdKP(df=data4pLV, dataType = "continuous", alpha = 0.05)
mat4pLV <- CreateAdjMat(ccd_4pLV, 4)
# Estimate PAG
pag_ccd4pLV <- plotPAG(ccd_4pLV, mat4pLV)

## Run FCI algorithm
suffStat_4pLV = list()
suffStat_4pLV$C = cor(data4pLV)
suffStat_4pLV$n = 1e6

res4pLV <- fci(suffStat_4pLV,indepTest=gaussCItest,
             alpha = 0.05, doPdsep = FALSE, labels = colnames(data4pLV))

pag_fci4pLV <- plotAG(res4pLV@amat)

## Run CCI algorithm
G4pLV = cci(suffStat_4pLV, gaussCItest, alpha=0.05, p=ncol(data4pLV)) 
dimnames(G4pLV$maag) <- list(colnames(data4pLV), colnames(data4pLV)) # give the labels
pag_cci4pLV <- plotAG(G4pLV$maag)


## ========================
## Model 7-2) 4 nodes with a LV - dense
## ========================
# specify B matrix
B4_LVdense = matrix(c(0, 0, 0, 0, 0.4,
                 0.8, 0, 0.6, 0, 0.6,
                 0, 0.6, 0, 0.9, 0,
                 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0), 5, 5, byrow = T)

colnames(B4_LVdense) <- c("X1", "X2", "X3", "X4", "L1")
# specify layout
layout4LVdense = matrix(c(-1,1,
                     -1,0,
                     1,0,
                     1,1,
                     -2, 0.5),5,2,byrow = T)
## True graph
true4pdenseLV <- qgraph(t(B4_LVdense), layout=layout4LVdense, labels = colnames(B4_LVdense), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B4_LVdense)
# generate data
data4pLVdense <- gen_dat(B4_LVdense, N =1e6, seed = 123)[,-5]
## Estimate GGM
ggm4pLVdense <- qgraph(cor(data4pLVdense), layout=layout4, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_4pLVdense <- ccdKP(df=data4pLVdense, dataType = "continuous", alpha = 0.05)
mat4pLVdense <- CreateAdjMat(ccd_4pLVdense, 4)
# Estimate PAG
pag_ccd4pLVdense <- plotPAG(ccd_4pLVdense, mat4pLVdense)

## Run FCI algorithm
suffStat_4pLVdense = list()
suffStat_4pLVdense$C = cor(data4pLVdense)
suffStat_4pLVdense$n = 1e6

res4pLVdense <- fci(suffStat_4pLVdense,indepTest=gaussCItest,
               alpha = 0.05, labels = colnames(data4pLV), selectionBias = FALSE)

pag_fci4pLVdense <- plotAG(res4pLVdense@amat)

## Run CCI algorithm
G4pLVdense = cci(suffStat_4pLVdense, gaussCItest, alpha=0.05, p=ncol(data4pLVdense)) 
dimnames(G4pLVdense$maag) <- list(colnames(data4pLVdense), colnames(data4pLVdense)) # give the labels
pag_cci4pLVdense <- plotAG(G4pLVdense$maag)

## ========================
## Model 7-3) 4 nodes with a LV - sparse
## ========================
# specify B matrix
B4_LV2 = matrix(c(0, 0, 0, 0, 1,
                 0, 0, 0.6, 0, 1,
                 0, 0.6, 0, 0.9, 0,
                 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0), 5, 5, byrow = T)

colnames(B4_LV2) <- c("X1", "X2", "X3", "X4", "L1")
# specify layout
layout4LV2 = matrix(c(-1,1,
                     -1,0,
                     1,0,
                     1,1,
                     -2, 0.5),5,2,byrow = T)
## True graph
true4pLV2 <- qgraph(t(B4_LV2), layout=layout4LV2, labels = colnames(B4_LV2), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B4_LV2)
# generate data
data4pLV2 <- gen_dat(B4_LV2, N =1e6, seed = 123)[,-5]
## Estimate GGM
ggm4pLV2 <- qgraph(cor(data4pLV2), layout=layout4, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_4pLV2 <- ccdKP(df=data4pLV2, dataType = "continuous", alpha = 0.05)
mat4pLV2 <- CreateAdjMat(ccd_4pLV2, 4)
# Estimate PAG
pag_ccd4pLV2 <- plotPAG(ccd_4pLV2, mat4pLV2)

## Run FCI algorithm
suffStat_4pLV2 = list()
suffStat_4pLV2$C = cor(data4pLV2)
suffStat_4pLV2$n = 1e6

res4pLV2 <- fci(suffStat_4pLV2,indepTest=gaussCItest,
               alpha = 0.05, selectionBias = FALSE, labels = colnames(data4pLV2))

pag_fci4pLV2 <- plotAG(res4pLV2@amat)

## Run CCI algorithm
G4pLV2 = cci(suffStat_4pLV2, gaussCItest, alpha=0.05, p=ncol(data4pLV2)) 
dimnames(G4pLV2$maag) <- list(colnames(data4pLV2), colnames(data4pLV2)) # give the labels
pag_cci4pLV2 <- plotAG(G4pLV2$maag)


## ========================
## Model 8-1) 5 nodes with a LV - sparse
## ========================
# specify B matrix

B5_lv = matrix(c(0, 0, 0, 0, 0, 1,
              0, 0, 0.4, 0, 0, 1,
              0, 0, 0, 0.5, 0,0,
              0, 0.7, 0, 0, 1.5,0,
              0, 0, 0, 0, 0,0,
              0,0,0,0,0,0), 6, 6, byrow = T)

colnames(B5_lv) <- c("X1", "X2", "X3", "X4", "X5", "L1")
# specify layout
layout5_lv = matrix(c(0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1,
                   -1, 0.5),6,2,byrow = T)

true5p_lv <- qgraph(t(B5_lv), layout=layout5_lv, labels = colnames(B5_lv), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5_lv)
# generate data and exclude the LV
data5pLV <- gen_dat(B5_lv, N =1e6, seed = 123)[,-6]


## Estimate GGM
ggm5pLV <- qgraph(cor(data5pLV), layout = layout5, theme="colorblind",graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_5pLV <- ccdKP(df=data5pLV, dataType = "continuous", alpha = 0.05)
mat5pLV <- CreateAdjMat(ccd_5pLV, 5)

## Estimate PAG
pag_ccd5pLV <- plotPAG(ccd_5pLV, mat5pLV)

## Run FCI algorithm
suffStat_5pLV = list()
suffStat_5pLV$C = cor(data5pLV)
suffStat_5pLV$n = 1e6

res5pLV <- fci(suffStat_5pLV,indepTest=gaussCItest,
                  alpha = 0.05, doPdsep = FALSE, labels = colnames(data5pLV))

pag_fci5pLV <- plotAG(res5pLV@amat)

## Run CCI algorithm
G5pLV = cci(suffStat_5pLV, gaussCItest, alpha=0.05, p=ncol(data5pLV)) 
dimnames(G5pLV$maag) <- list(colnames(data5pLV), colnames(data5pLV)) # give the labels
pag_cci5pLV <- plotAG(G5pLV$maag)


## ========================
## Model 8-2) 5 nodes with a LV - dense
## ========================
# specify B matrix

B5_lv_dense = matrix(c(0, 0, 0, 0, 0, 1,
                 0, 0, 0.4, 0, 0, 1,
                 0, 0, 0, 0.5, 0,0,
                 0, 0.7, 0, 0, 1.5,0,
                 0.6, 0, 0, 0, 0,0,
                 0,0,0,0,0,0), 6, 6, byrow = T)

colnames(B5_lv_dense) <- c("X1", "X2", "X3", "X4", "X5", "L1")
# specify layout
layout5_lv = matrix(c(0,1,
                      0,0,
                      1,-1,
                      2,0,
                      2,1,
                      -1, 0.5),6,2,byrow = T)

true5p_lvdense <- qgraph(t(B5_lv_dense), layout=layout5_lv, labels = colnames(B5_lv_dense), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5_lv_dense)
# generate data and exclude the LV
data5pLVdense <- gen_dat(B5_lv_dense, N =1e6, seed = 123)[,-6]


## Estimate GGM
ggm5pLVdense <- qgraph(cor(data5pLVdense), layout = layout5, theme="colorblind",graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_5pLVdense <- ccdKP(df=data5pLVdense, dataType = "continuous", alpha = 0.05)
mat5pLVdense <- CreateAdjMat(ccd_5pLVdense, 5)

## Estimate PAG
pag_ccd5pLVdense <- plotPAG(ccd_5pLVdense, mat5pLVdense)

## Run FCI algorithm
suffStat_5pLVdense = list()
suffStat_5pLVdense$C = cor(data5pLVdense)
suffStat_5pLVdense$n = 1e6

res5pLVdense <- fci(suffStat_5pLVdense,indepTest=gaussCItest,
               alpha = 0.05, doPdsep = FALSE, labels = colnames(data5pLVdense))

pag_fci5pLVdense <- plotAG(res5pLVdense@amat)

## Run CCI algorithm
G5pLVdense = cci(suffStat_5pLVdense, gaussCItest, alpha=0.05, p=ncol(data5pLVdense)) 
dimnames(G5pLVdense$maag) <- list(colnames(data5pLVdense), colnames(data5pLVdense)) # give the labels
pag_cci5pLVdense <- plotAG(G5pLVdense$maag)




## ========================
## Model 8-3) 5 nodes with a LV - sparse
## ========================
# specify B matrix

B5_lv2 = matrix(c(0, 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 0, 1,
                 0, 0.5, 0, 0, 0.6, 0,
                 1, 0, 0, 0, 0.5, 0,
                 0, 0, 0, 0.5, 0,0,
                 0,0,0,0,0,0), 6, 6, byrow = T)

colnames(B5_lv2) <- c("X1", "X2", "X3", "X4", "X5", "L1")
# specify layout
layout5_lv2 = matrix(c(0,1,
                      0,0,
                      1,0.5,
                      2,1,
                      2,0,
                      -1, 0.5),6,2,byrow = T)

true5p_lv2 <- qgraph(t(B5_lv2), layout=layout5_lv2, labels = colnames(B5_lv2), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5_lv2)
# generate data and exclude the LV
data5pLV2 <- gen_dat(B5_lv2, N =1e6, seed = 12345)[,-6]


## Estimate GGM
ggm5pLV2 <- qgraph(cor(data5pLV2), layout = layout5, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_5pLV2 <- ccdKP(df=data5pLV2, dataType = "continuous", alpha = 0.05)
mat5pLV2 <- CreateAdjMat(ccd_5pLV2, 5)
## Estimate PAG
pag_ccd5pLV2 <- plotPAG(ccd_5pLV2, mat5pLV2)

## Run FCI algorithm
suffStat_5pLV2 = list()
suffStat_5pLV2$C = cor(data5pLV2)
suffStat_5pLV2$n = 1e6

res5pLV2 <- fci(suffStat_5pLV2,indepTest=gaussCItest,
               alpha = 0.05, doPdsep = FALSE, labels = colnames(data5pLV2))

pag_fci5pLV2 <- plotAG(res5pLV2@amat)

## Run CCI algorithm
G5pLV2 = cci(suffStat_5pLV2, gaussCItest, alpha=0.05, p=ncol(data5pLV2)) 
dimnames(G5pLV2$maag) <- list(colnames(data5pLV2), colnames(data5pLV2)) # give the labels
pag_cci5pLV2 <- plotAG(G5pLV2$maag)



## ========================
## Model 9-1) 10 nodes with a LV - sparse
## ========================
# specify B matrix

B10_lv = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0.4, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0.2, 0, 0.5, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1, 0.8, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.8,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 11, 11, byrow = T)

colnames(B10_lv) <- c(paste("X", 1:10, sep=""), "L1")

# specify layout
layout10LV = matrix(c(0,1,
                      2,1,
                      1,0,
                      2,-1,
                      3,0,
                      4, -1,
                      5, 0,
                      6, -1,
                      4, 1,
                      7, 1,
                      8, 0),11,2,byrow = T)

true10pLV <- qgraph(t(B10_lv), layout = layout10LV, labels = colnames(B10_lv), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B10_lv)
# generate data and exclude the LV
data10pLV <- gen_dat(B10_lv, N =1e6, seed = 123)[,-11]


## Estimate GGM
ggm10pLV <- qgraph(cor(data10pLV), layout = layout10LV, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_10pLV <- ccdKP(df=data10pLV, dataType = "continuous", alpha = 0.05)
mat10pLV <- CreateAdjMat(ccd_10pLV, 10)
## Estimate PAG
pag_ccd10pLV <- plotPAG(ccd_10pLV, mat10pLV)

## Run FCI algorithm
suffStat_10pLV = list()
suffStat_10pLV$C = cor(data10pLV)
suffStat_10pLV$n = 1e6

res10pLV <- fci(suffStat_10pLV,indepTest=gaussCItest,
                alpha = 0.05, doPdsep = FALSE, labels = colnames(data10pLV))

pag_fci10pLV <- plotAG(res10pLV@amat)

## Run CCI algorithm
G10pLV = cci(suffStat_10pLV, gaussCItest, alpha=0.05, p=ncol(data10pLV)) 
dimnames(G10pLV$maag) <- list(colnames(data10pLV), colnames(data10pLV)) # give the labels
pag_cci10pLV <- plotAG(G10pLV$maag)


## ========================
## Model 9-2) 10 nodes with a LV - dense
## ========================
# specify B matrix

B10_lvdense = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0.4, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 0,
                  0, 0, 0.5, 1, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0.2, 0, 0.5, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1, 0.8, 0.5, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.8,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 11, 11, byrow = T)

colnames(B10_lvdense) <- c(paste("X", 1:10, sep=""), "L1")

# specify layout
layout10LV = matrix(c(0,1,
                      2,1,
                      1,0,
                      2,-1,
                      3,0,
                      4, -1,
                      5, 0,
                      6, -1,
                      4, 1,
                      7, 1,
                      8, 0),11,2,byrow = T)

true10pLVdense <- qgraph(t(B10_lvdense), layout = layout10LV, labels = colnames(B10_lv), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B10_lvdense)
# generate data and exclude the LV
data10pLVdense <- gen_dat(B10_lvdense, N =1e6, seed = 123)[,-11]


## Estimate GGM
ggm10pLVdense <- qgraph(cor(data10pLVdense), layout = layout10LV, theme="colorblind", graph = "pcor")

layout(t(1:3))
## Run CCD algorithm
ccd_10pLVdense <- ccdKP(df=data10pLVdense, dataType = "continuous", alpha = 0.05)
mat10pLVdense <- CreateAdjMat(ccd_10pLVdense, 10)
## Estimate PAG
pag_ccd10pLVdense <- plotPAG(ccd_10pLVdense, mat10pLVdense)

## Run FCI algorithm
suffStat_10pLVdense = list()
suffStat_10pLVdense$C = cor(data10pLVdense)
suffStat_10pLVdense$n = 1e6

res10pLVdense <- fci(suffStat_10pLVdense,indepTest=gaussCItest,
                alpha = 0.05, doPdsep = FALSE, labels = colnames(data10pLVdense))

pag_fci10pLVdense <- plotAG(res10pLVdense@amat)

## Run CCI algorithm
G10pLVdense = cci(suffStat_10pLVdense, gaussCItest, alpha=0.05, p=ncol(data10pLVdense)) 
dimnames(G10pLVdense$maag) <- list(colnames(data10pLVdense), colnames(data10pLVdense)) # give the labels
pag_cci10pLVdense <- plotAG(G10pLVdense$maag)



#####################
## McNally example ##
#####################
# empirical data 408 rows by 26 columns (p = 26)
mcnally <- read.csv("data/McNally.csv") 
# check the data
dplyr::glimpse(mcnally)
skimr::skim(mcnally)

# run ccd algorithm (discrete)
ccd_mcnally <- ccdKP(df=mcnally, dataType = "discrete", depth = -1) 

ccd_mcnally$graph 
ccd_mcnally$edges
graph_mcnally <- tetradrunner.tetradGraphToDot(ccd_mcnally$graph)
dot(graph_mcnally)

mat_mcnally <- CreateAdjMat(ccd_mcnally, p = 26)
plotPAG(ccd_mcnally, mat_mcnally) # not pretty

# separate dep / ocd symptoms
depression <- mcnally[,1:16]
ocd <- mcnally[,17:26]

layout(1)
## run ccd on depression symptoms
ccd_mcnally_dep <- ccdKP(df=depression, dataType = "discrete", depth = -1) 
mat_mcnally_dep <- CreateAdjMat(ccd_mcnally_dep, p = ncol(depression))
plotPAG(ccd_mcnally_dep, mat_mcnally_dep) 

## run fci on depression symptoms
fci_mcdep <- tetradrunner(algoId = 'fci', df = depression,
                          dataType = 'discrete', verbose=TRUE)
fci_mcdep$graph
fci_mcdep$edges
fci_mcdep$graph$getUnderLines()

matfci_mcdep <- CreateAdjMat(fci_mcdep, ncol(depression))
plotPAG(fci_mcdep, matfci_mcdep)


# cant use discrete test with pcalg
# # disCItest not working (n is too small)
# suffStat <- list(dm = depression, nlev = rep(4, length(depression)), adaptDF = FALSE)
# Gdepression = cci(suffStat, disCItest, alpha=0.05, p=ncol(depression)) 
# 
# # not working with gaussCItest well
# ## Run CCI on depression symptoms
# suffStat_dep = list()
# suffStat_dep$C = cor(depression)
# suffStat_dep$n = 1e5 #nrow(depression)
# 
# Gdepression = cci(suffStat_dep, gaussCItest, alpha=0.05, p=ncol(depression)) 
# dimnames(Gdepression$maag) <- list(colnames(depression), colnames(depression)) # give the labels
# pag_ccidep <- plotAG(Gdepression$maag) # plot looks dirty
# 

## run ccd on ocd symptoms
ccd_mcnally_ocd <- ccdKP(df=ocd, dataType = "discrete", depth = -1) 
mat_mcnally_ocd <- CreateAdjMat(ccd_mcnally_ocd, p = ncol(ocd))
plotPAG(ccd_mcnally_ocd, mat_mcnally_ocd) 


## run fci on ocd symptoms
fci_mcocd <- tetradrunner(algoId = 'fci', df = ocd,
                          dataType = 'discrete', verbose=TRUE)
fci_mcocd$graph
fci_mcocd$edges
fci_mcocd$graph$getUnderLines()

matfci_mcocd <- CreateAdjMat(fci_mcocd, ncol(ocd))
plotPAG(fci_mcocd, matfci_mcocd)
