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
ggm4p <- qgraph(cor(data4p), layout=layout4, theme="colorblind")

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
suffStat_4p$n = 1000

res4p <- fci(suffStat_4p,indepTest=gaussCItest,
             alpha = 0.05, doPdsep = FALSE, labels = colnames(data4p))

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
ggm4p_high <- qgraph(t(cor(data4p_high)), layout=layout4, theme="colorblind")

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
suffStat_4phigh$n = 1000

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
ggm5p <- qgraph(cor(data5p), layout = layout5, theme="colorblind")

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
suffStat_5p$n = 1000

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
ggm5p_high <- qgraph(cor(data5p_high), layout = layout5, theme="colorblind")

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
suffStat_5phigh$n = 1000

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
ggm6p <- qgraph(cor(data6p), layout = layout6, theme="colorblind")

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
suffStat_6p$n = 1000

res6p <- fci(suffStat_6p,indepTest=gaussCItest,
                  alpha = 0.05, doPdsep = FALSE, labels = colnames(data6p))

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
ggm6p_high <- qgraph(cor(data6p_high), layout = layout6, theme="colorblind")

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
suffStat_6phigh$n = 1000

res6p_high <- fci(suffStat_6phigh,indepTest=gaussCItest,
             alpha = 0.05, doPdsep = FALSE, labels = colnames(data6p_high))

pag_fci6pH <- plotAG(res6p_high@amat)

## Run CCI algorithm
G6p_high = cci(suffStat_6phigh, gaussCItest, alpha=0.05, p=ncol(data6p_high)) 
dimnames(G6p_high$maag) <- list(colnames(data6p_high), colnames(data6p_high)) # give the labels
pag_cci6pH <- plotAG(G6p_high$maag)


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
