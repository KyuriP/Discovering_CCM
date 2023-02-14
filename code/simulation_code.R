## load functions & packages
# source("code/R/CCD_fnc.R")
# source("code/R/plot_fnc.R")
# source("code/R/data_generating_fnc.R")
# source("code/eval_metrics.R")
# source("../code/true_ancestral.R")


library(magrittr)
library(purrr)
library(furrr)
library(dplyr)
library(qgraph)
library(ggplot2)

## slightly modified CCI package
#install_github("KyuriP/CCI_KP")
library(CCI.KP)

## simulate data 1000 times ## 
set.seed(123)
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

true5psparse <- qgraph(t(B5sparse), layout=layout5, labels = colnames(B5sparse), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5sparse)

# generate data 
# specify the sample sizes
N <- c(50, 150, 500, 1000, 5000) 

simdata_5psparse <- N %>% future_map(function(z) {
  replicate(n=1000,
            expr = gen_dat(B5sparse, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


# data5psparse <- gen_dat(B5sparse, N =1e6, seed = 123)

# ## Estimate GGM
# ggm5p <- qgraph(cor(data5psparse), layout = layout5, theme="colorblind", graph = "pcor")
# 
# layout(t(1:3))

## True Ancestral Graph
dcg_5psparse <- matrix(c(0,1,0,0,0,
                         0,0,0,1,0,
                         0,1,0,0,0,
                         0,0,1,0,0,
                         0,0,0,1,0), 5,5,byrow=T)
trueag_5psparse <- true_ancestral(dcg_5psparse, gen_dat(B5sparse), gaussCItest)
dimnames(trueag_5psparse) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
plotAG(trueag_5psparse)

## Run CCD algorithm
# ccd_5psparse <- simdata_5psparse %>% 
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
#       ) 
mat_5psparse <- ccd_5psparse %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_5psparse, file="data/ccd_5psparse.RData")
# load("data/ccd_5psparse.RData")

# pag_ccd5psparse <- map2(ccd_5psparse, mat_5psparse, plotPAG)


# ## Run FCI algorithm
# fci_5psparse <- simdata_5psparse %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#            alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # extract amat
#       )

# save(fci_5psparse, file="data/fci_5psparse.RData")
# load("data/fci_5psparse.RData")

# pag_fci5psparse <- fci_5psparse %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
# cci_5psparse <- simdata_5psparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#       )

# save(cci_5psparse, file="data/cci_5psparse.RData")
# load("data/cci_5psparse.RData")

# pag_cci5psparse <- cci_5psparse %>% 
#   map(~plotAG(.x))



## evaluation 
# CCD
res_ccd5psparse <- mat_5psparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_ccd5psparse <- mat_5psparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% do.call("cbind", .) %>% apply(., 2, unlist) %>%  as.data.frame %>% 
  rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd5psparse, na.rm=T)

# SHD
SHD_ccd5psparse <- mat_5psparse %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% do.call("cbind", .) %>% apply(., 2, unlist) %>%  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd5psparse)

# FCI
res_fci5psparse <- fci_5psparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci5psparse <- fci_5psparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_fci5psparse, na.rm=T)

# SHD
SHD_fci5psparse <- fci_5psparse %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci5psparse)

# CCI
res_cci5psparse <- cci_5psparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_cci5psparse <- cci_5psparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_cci5psparse)

# SHD
SHD_cci5psparse <- cci_5psparse %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_cci5psparse)


## ====================
## 5p - dense
## ====================
# specify B matrix

B5dense = matrix(c(0, 0, 0, 0, 0,
                    1, 0, 0.8, 0, 0,
                    0, 0, 0, 0.9, 0,
                    0, 0.7, 0, 0, 1.5,
                    1, 0, 0, 0, 0), 5, 5, byrow = T)

colnames(B5dense) <- c("X1", "X2", "X3", "X4", "X5")

true5pdense <- qgraph(t(B5dense), layout=layout5, labels = colnames(B5dense), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5dense)

# generate data (sample size as specified above)
simdata_5pdense <- N %>% future_map(function(z) {
  replicate(n=1000,
            expr = gen_dat(B5dense, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))

# # generate data and exclude the LV
# data5pdense <- gen_dat(B5dense, N =1e6, seed = 12345)
# 
# ## Estimate GGM
# ggm5pdense <- qgraph(cor(data5pdense), layout = layout5, theme="colorblind", graph = "pcor")
# 
# layout(t(1:3))


## True Ancestral Graph
dcg_5pdense <- matrix(c(0,1,0,0,1,
                        0,0,0,1,0,
                        0,1,0,0,0,
                        0,0,1,0,0,
                        0,0,0,1,0), 5,5,byrow=T)

trueag_5pdense <- true_ancestral(dcg_5pdense, gen_dat(B5dense), gaussCItest)
dimnames(trueag_5pdense) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
plotAG(trueag_5pdense)


## Run CCD algorithm
# ccd_5pdense <- simdata_5pdense %>% 
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
#   ) 
mat_5pdense <- ccd_5pdense %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_5pdense, file="data/ccd_5pdense.RData")
# load("data/ccd_5pdense.RData")

# pag_ccd5pdense <- map2(ccd_5pdense, mat_5pdense, plotPAG)


## Run FCI algorithm
# fci_5pdense <- simdata_5pdense %>% 
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#            alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # extract amat
#   )

# save(fci_5pdense, file="data/fci_5pdense.RData")
# load("data/fci_5pdense.RData")


# pag_fci5pdense <- fci_5pdense %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
# cci_5pdense <- simdata_5pdense %>% 
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )

# save(cci_5pdense, file="data/cci_5pdense.RData")
# load("data/cci_5pdense.RData")

# pag_cci5pdense <- cci_5pdense %>% 
#   map(~plotAG(.x))


## evaluation
# CCD
res_ccd5pdense <- mat_5pdense %>% 
  map_depth(2, ~precision_recall(trueag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame()

# UNCERTAINTY
uncer_ccd5pdense <- mat_5pdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_ccd5pdense, na.rm=T)

# SHD
SHD_ccd5pdense <- mat_5pdense %>% 
  map_depth(2, ~SHD(trueag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(SHD_ccd5pdense)

# FCI
res_fci5pdense <- fci_5pdense %>% 
  map_depth(2, ~precision_recall(trueag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci5pdense <- fci_5pdense%>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_fci5pdense, na.rm=T)

# SHD
SHD_fci5pdense <- fci_5pdense %>% 
  map_depth(2, ~SHD(trueag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(SHD_fci5pdense)

# CCI
res_cci5pdense <- cci_5pdense %>% 
  map_depth(2, ~precision_recall(trueag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY 
uncer_cci5pdense <- cci_5pdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_cci5pdense)

# SHD
SHD_cci5pdense <- cci_5pdense %>% 
  map_depth(2, ~SHD(trueag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(SHD_cci5pdense)


## ====================
## 10p - sparse
## ====================

B10sparse = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                  0.4, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 
                  0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 
                  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                  0, 0, 0, 0, 0.2, 0, 0.5, 0, 0, 0, 
                  0, 0, 0, 0, 0, 0, 0, 1, 0.8, 0, 
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4, 
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
 
true10psparse <- qgraph(t(B10sparse), layout = layout10, labels = colnames(B10sparse), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B10sparse)

# generate data (sample size as specified above)
simdata_10psparse <- N %>% future_map(function(z) {
  replicate(n=1000,
            expr = gen_dat(B10sparse, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


# # generate data and exclude the LV
# data10psparse <- gen_dat(B10sparse, N =1e6, seed = 123)
# 
# ## Estimate GGM
# ggm10psparse <- qgraph(cor(data10psparse), layout = layout10, theme="colorblind", graph = "pcor")
# layout(t(1:3))

## True Ancestral Graph
dcg_10psparse <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                          0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                          0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                          0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                          0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                          0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                          0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                          0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                          0, 0, 0, 0, 0, 0, 0, 0, 1, 0), 10, 10, byrow = T)

trueag_10psparse <- true_ancestral(dcg_10psparse, gen_dat(B10sparse), gaussCItest)
dimnames(trueag_10psparse) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(trueag_10psparse)


## Run CCD algorithm
# ccd_10psparse <- simdata_10psparse %>% 
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
#   ) 
mat_10psparse <- ccd_10psparse %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_10psparse, file="data/ccd_10psparse.RData")
# load("data/ccd_10psparse.RData")

# pag_ccd10psparse <- map2(ccd_10psparse, mat_10psparse, plotPAG)


## Run FCI algorithm
# fci_10psparse <- simdata_10psparse %>% 
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#            alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
#   )

# save(fci_10psparse, file="data/fci_10psparse.RData")
# load("data/fci_10psparse.RData")

# pag_fci10psparse <- fci_10psparse %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
# cci_10psparse <- simdata_10psparse %>% 
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )

# save(cci_10psparse, file="data/cci_10psparse.RData")
# load("data/cci_10psparse.RData")

# pag_cci10psparse <- cci_10psparse %>% 
#   map(~plotAG(.x))



## evaluation
# CCD
res_ccd10psparse <- mat_10psparse %>% 
  map_depth(2, 
  ~precision_recall(trueag_10psparse, .x)) %>%
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame()
  

# UNCERTAINTY
uncer_ccd10psparse <- mat_10psparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd10psparse, na.rm=T)

# SHD
SHD_ccd10psparse <- mat_10psparse %>% 
  map_depth(2, ~SHD(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd10psparse)

# FCI
res_fci10psparse <- fci_10psparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci10psparse <- fci_10psparse%>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_fci10psparse, na.rm=T)

# SHD
SHD_fci10psparse <- fci_10psparse %>% 
  map_depth(2, ~SHD(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(SHD_fci10psparse)

# CCI
res_cci10psparse <- cci_10psparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_cci10psparse <- cci_10psparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_cci10psparse, na.rm=T)

# SHD
SHD_cci10psparse <- cci_10psparse %>% 
  map_depth(2, ~SHD(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_cci10psparse)



## ====================
## 10p - dense
## ====================

B10dense = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                     0.4, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 
                     0, 0.4, 0, 1, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0.9, 0, 0.5, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 0, 1, 0.8, 1, 
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4, 
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
layout(1)
true10pdense <- qgraph(t(B10dense), layout = layout10, labels = colnames(B10dense), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B10dense)

# generate data (sample size as specified above)
simdata_10pdense <- N %>% future_map(function(z) {
  replicate(n=1000,
            expr = gen_dat(B10dense, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))

# ## Estimate GGM
# ggm10pdense <- qgraph(cor(data10pdense), layout = layout10, theme="colorblind", graph = "pcor")
# layout(t(1:3))

## True Ancestral Graph
dcg_10pdense <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 
                         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 1, 1, 1, 0), 10, 10, byrow = T)

trueag_10pdense <- true_ancestral(dcg_10pdense, gen_dat(B10dense), gaussCItest)
dimnames(trueag_10pdense) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(trueag_10pdense)


## Run CCD algorithm
# ccd_10pdense <- simdata_10pdense  %>% 
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
#   ) 
mat_10pdense  <- ccd_10pdense  %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_10pdense, file="data/ccd_10pdense.RData")
# load("data/ccd_10pdense.RData")

# pag_ccd10pdense  <- map2(ccd_10pdense , mat_10pdense , plotPAG)


## Run FCI algorithm
# fci_10pdense <- simdata_10pdense  %>% 
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#            alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
#   )

# save(fci_10pdense, file="data/fci_10pdense.RData")
# load("data/fci_10pdense.RData")

# pag_fci10pdense <- fci_10pdense  %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
# cci_10pdense  <- simdata_10pdense %>% 
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )

# save(cci_10pdense, file="data/cci_10pdense.RData")
# load("data/cci_10pdense.RData")

# pag_cci10pdense  <- cci_10pdense  %>% 
#   map(~plotAG(.x))


## evaluation
# CCD
res_ccd10pdense  <- mat_10pdense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_ccd10pdense  <- mat_10pdense  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd10pdense , na.rm=T)

# SHD
SHD_ccd10pdense  <- mat_10pdense %>% 
  map_depth(2, ~SHD(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd10pdense)

# FCI
res_fci10pdense  <- fci_10pdense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci10pdense <- fci_10pdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_fci10pdense, na.rm=T)

# SHD
SHD_fci10pdense  <- fci_10pdense %>% 
  map_depth(2, ~SHD(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci10pdense)

# CCI
res_cci10pdense <- cci_10pdense %>% 
  map_depth(2, ~precision_recall(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_cci10pdense <- cci_10pdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_cci10pdense, na.rm=T)

# SHD
SHD_cci10pdense <- cci_10pdense %>% 
  map_depth(2, ~SHD(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# averae SHD
colMeans(SHD_cci10pdense)


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

true5p_lvsparse <- qgraph(t(B5_lvsparse), layout=layout5_lv, labels = colnames(B5_lvsparse), theme="colorblind")
 
# B5_lv2 = matrix(c(0, 0, 0, 0, 0, 1,
#                   0, 0, 0, 0, 0, 1,
#                   0, 0.5, 0, 0, 0.6, 0,
#                   1, 0, 0, 0, 0.5, 0,
#                   0, 0, 0, 0.5, 0,0,
#                   0,0,0,0,0,0), 6, 6, byrow = T)
# 
# colnames(B5_lv2) <- c("X1", "X2", "X3", "X4", "X5", "L1")
# # specify layout
# layout5_lv2 = matrix(c(0,1,
#                        0,0,
#                        1,0.5,
#                        2,1,
#                        2,0,
#                        -1, 0.5),6,2,byrow = T)
# 
# true5p_lv2 <- qgraph(t(B5_lv2), layout=layout5_lv2, labels = colnames(B5_lv2), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5_lvsparse)

# generate data (sample size as specified above)
simdata_5pLVsparse <- N %>% future_map(function(z) {
  replicate(n=1000,
            expr = gen_dat(B5_lvsparse, N = z)[,-6],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


# # generate data and exclude the LV
# data5pLV2 <- gen_dat(B5_lv2, N =1e6, seed = 12345)[,-6]


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
# ccd_5pLVsparse <- simdata_5pLVsparse  %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
#   )
mat_5pLVsparse  <- ccd_5pLVsparse  %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_5pLVsparse, file="data/ccd_5pLVsparse.RData")
# load("data/ccd_5pLV.RData")

# pag_ccd5pLV2  <- map2(ccd_5pLV2 , mat_5pLV2 , plotPAG)


## Run FCI algorithm
# fci_5pLVsparse <- simdata_5pLVsparse  %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#            alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
#   )

# save(fci_5pLVsparse, file="data/fci_5pLVsparse.RData")
# load("data/fci_5pLV2.RData")

# pag_fci5pLV2 <- fci_5pLV2  %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
# cci_5pLVsparse  <- simdata_5pLVsparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )
# 
# save(cci_5pLVsparse, file="data/cci_5pLVsparse.RData")
# load("data/cci_5pLV2.RData")

# pag_cci5pLV2  <- cci_5pLV2  %>% 
#   map(~plotAG(.x))



## evaluation
# CCD
res_ccd5pLVsparse  <- mat_5pLVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_ccd5pLVsparse  <- mat_5pLVsparse  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd5pLVsparse , na.rm=T)

# SHD
SHD_ccd5pLVsparse <- mat_5pLVsparse %>% 
  map_depth(2, ~SHD(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd5pLVsparse)

# FCI
res_fci5pLVsparse  <- fci_5pLVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci5pLVsparse <- fci_5pLVsparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_fci5pLVsparse, na.rm=T)

# SHD
SHD_fci5pLVsparse <- fci_5pLVsparse %>% 
  map_depth(2, ~SHD(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci5pLVsparse)

# CCI
res_cci5pLVsparse <- cci_5pLVsparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_cci5pLVsparse <- cci_5pLVsparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average unceratinty
colMeans(uncer_cci5pLVsparse, na.rm=T)

# SHD
SHD_cci5pLVsparse <- cci_5pLVsparse %>% 
  map_depth(2, ~SHD(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_cci5pLVsparse)


## ====================
## 5p with LV dense
## ====================
# specify B matrix

B5_lvdense = matrix(c(0, 0, 0, 0, 0, 1,
                       0, 0, 0.4, 0, 0, 1,
                       0, 0, 0, 0.5, 0,0,
                       0, 0.7, 0, 0, 1.5,0,
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

true5p_lvdense <- qgraph(t(B5_lvdense), layout=layout5_lv, labels = colnames(B5_lvdense), theme="colorblind")


## Data generating
# equilibrium check
equilibrium_check(B5_lvdense)

# generate data (sample size as specified above)
simdata_5pLVdense <- N %>% future_map(function(z) {
  replicate(n=1000,
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
# ccd_5pLVdense <- simdata_5pLVdense  %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
#   )
mat_5pLVdense  <- ccd_5pLVdense  %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_5pLVdense, file="data/ccd_5pLVdense.RData")
# load("data/ccd_5pLV.RData")

# pag_ccd5pLV2  <- map2(ccd_5pLV2 , mat_5pLV2 , plotPAG)


## Run FCI algorithm
# fci_5pLVdense <- simdata_5pLVdense  %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                     alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
#   )
# 
# save(fci_5pLVdense, file="data/fci_5pLVdense.RData")
# load("data/fci_5pLV2.RData")

# pag_fci5pLVdense <- fci_5pLVdense  %>%
#   map_depth(2, ~plotAG(.x))

## Run CCI algorithm
# cci_5pLVdense  <- simdata_5pLVdense %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )
# 
# save(cci_5pLVdense, file="data/cci_5pLVdense.RData")
# load("data/cci_5pLV2.RData")

# pag_cci5pLV2  <- cci_5pLV2  %>% 
#   map(~plotAG(.x))



## evaluation
# CCD
res_ccd5pLVdense  <- mat_5pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_5pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_ccd5pLVdense  <- mat_5pLVdense  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd5pLVdense , na.rm=T)

# SHD
SHD_ccd5pLVdense <- mat_5pLVdense %>% 
  map_depth(2, ~SHD(trueag_5pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd5pLVdense)

# FCI
res_fci5pLVdense  <- fci_5pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_5pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci5pLVdense <- fci_5pLVdense  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_fci5pLVdense , na.rm=T)

# SHD
SHD_fci5pLVdense  <- fci_5pLVdense  %>% 
  map_depth(2, ~SHD(trueag_5pdenseLV , .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci5pLVdense )

# CCI
res_cci5pLVdense  <- cci_5pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_5pdenseLV , .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_cci5pLVdense  <- cci_5pLVdense  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average unceratinty
colMeans(uncer_cci5pLVdense , na.rm=T)

# SHD
SHD_cci5pLVdense  <- cci_5pLVdense  %>% 
  map_depth(2, ~SHD(trueag_5pdenseLV , .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_cci5pLVdense )



## ====================
## 10p with LV sparse
## ====================
# specify B matrix

B10_lvsparse = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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

colnames(B10_lvsparse) <- c(paste("X", 1:10, sep=""), "L1")

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

true10pLVsparse <- qgraph(t(B10_lvsparse), layout = layout10LV, labels = colnames(B10_lvsparse), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B10_lvsparse)

# generate data (sample size as specified above)
simdata_10pLVsparse <- N %>% future_map(function(z) {
  replicate(n=1000,
            expr = gen_dat(B10_lvsparse, N = z)[,-11],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))

# # generate data and exclude the LV
# data10pLV <- gen_dat(B10_lv, N =1e6, seed = 123)[,-11]
# 
# ## Estimate GGM
# ggm10pLV <- qgraph(cor(data10pLV), layout = layout10LV, theme="colorblind", graph = "pcor")
# 
# layout(t(1:3))

## True Ancestral Graph
# [i,j] = [j,i] = 2: a LV exists between i and j
dcg_10psparseLV <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                            0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 
                            0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                            0, 0, 0, 0, 0, 0, 0, 2, 1, 0), 10, 10, byrow = T)

trueag_10psparseLV <- true_ancestral(dcg_10psparseLV, gen_dat(B10_lvsparse), gaussCItest)
dimnames(trueag_10psparseLV) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(trueag_10psparseLV)


## Run CCD algorithm
# ccd_10pLVsparse  <- simdata_10pLVsparse   %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
#   )
mat_10pLVsparse   <- ccd_10pLVsparse %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_10pLVsparse, file="data/ccd_10pLVsparse.RData")
# load("data/ccd_10pLV.RData")

# pag_ccd10pLV <- map2(ccd_10pLV, mat_10pLV  , plotPAG)


## Run FCI algorithm
# fci_10pLVsparse  <- simdata_10pLVsparse   %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#            alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
#   )

# save(fci_10pLVsparse, file="data/fci_10pLVsparse.RData")
# load("data/fci_10pLV.RData")

# pag_fci10pLV  <- fci_10pLV   %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
# cci_10pLVsparse  <- simdata_10pLVsparse %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )

# save(cci_10pLVsparse, file="data/cci_10pLVsparse.RData")
# load("data/cci_10pLV.RData")

# pag_cci10pLV   <- cci_10pLV %>% 
#   map(~plotAG(.x))



## evaluation
# CCD
res_ccd10pLVsparse   <- mat_10pLVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_ccd10pLVsparse  <- mat_10pLVsparse  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd10pLVsparse , na.rm=T)

# SHD
SHD_ccd10pLVsparse <- mat_10pLVsparse %>% 
  map_depth(2, ~SHD(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd10pLVsparse)

# FCI
res_fci10pLVsparse <- fci_10pLVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci10pLVsparse <- fci_10pLVsparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_fci10pLVsparse, na.rm=T)


# SHD
SHD_fci10pLVsparse <- fci_10pLVsparse %>% 
  map_depth(2, ~SHD(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci10pLVsparse)


# CCI 
res_cci10pLVsparse <- cci_10pLVsparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 


# UNCERTAINTY
uncer_cci10pLVsparse <- cci_10pLVsparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_cci10pLVsparse, na.rm=T)

# SHD
SHD_cci10pLVsparse <- cci_10pLVsparse %>% 
  map_depth(2, ~SHD(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_cci10pLVsparse)




## ====================
## 10p with LV dense
## ====================
# specify B matrix

B10_lvdense = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0.4, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0.7, 0, 0, 0.9, 0, 0, 0, 0, 0,
                  0, 0, 0.6, 1, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0.2, 0, 0.5, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1, 0.8, 0.6, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.8,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 11, 11, byrow = T)

colnames(B10_lvdense) <- c(paste("X", 1:10, sep=""), "L1")

# specify layout
layout10LV = matrix(c(0, 1,
                      2, 1,
                      1, 0,
                      2, -1,
                      3, 0,
                      4, -1,
                      5, 0,
                      6, -1,
                      4, 1,
                      7, 1,
                      8, 0), 11, 2, byrow = T)

true10pLVdense <- qgraph(t(B10_lvdense), layout = layout10LV, labels = colnames(B10_lvdense), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B10_lvdense)

# generate data (sample size as specified above)
simdata_10pLVdense <- N %>% future_map(function(z) {
  replicate(n=1000,
            expr = gen_dat(B10_lvdense, N = z)[,-11],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))



## True Ancestral Graph
# [i,j] = [j,i] = 2: a LV exists between i and j
dcg_10pdenseLV <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 
                           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 1, 2, 1, 0), 10, 10, byrow = T)

trueag_10pdenseLV <- true_ancestral(dcg_10pdenseLV, gen_dat(B10_lvdense), gaussCItest)
dimnames(trueag_10pdenseLV) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(trueag_10pdenseLV)


## Run CCD algorithm
# ccd_10pLVdense  <- simdata_10pLVdense   %>%
#   map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
#   )
mat_10pLVdense   <- ccd_10pLVdense %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# save(ccd_10pLVdense, file="data/ccd_10pLVdense.RData")
# load("data/ccd_10pLV.RData")

# pag_ccd10pLV <- map2(ccd_10pLV, mat_10pLV  , plotPAG)


## Run FCI algorithm
# fci_10pLVdense  <- simdata_10pLVdense   %>%
#   map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#            alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
#   )

# save(fci_10pLVdense, file="data/fci_10pLVdense.RData")
# load("data/fci_10pLV.RData")

# pag_fci10pLV  <- fci_10pLV   %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
# cci_10pLVdense  <- simdata_10pLVdense  %>%
#   map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
#   )

# save(cci_10pLVdense, file="data/cci_10pLVdense.RData")
# load("data/cci_10pLV.RData")

# pag_cci10pLV   <- cci_10pLV %>% 
#   map(~plotAG(.x))



## evaluation
# CCD
res_ccd10pLVdense   <- mat_10pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_ccd10pLVdense  <- mat_10pLVdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd10pLVdense , na.rm=T)

# SHD
SHD_ccd10pLVdense <- mat_10pLVdense %>% 
  map_depth(2, ~SHD(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N)) 
# average SHD
colMeans(SHD_ccd10pLVdense)

# FCI
res_fci10pLVdense <- fci_10pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci10pLVdense <- fci_10pLVdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_fci10pLVdense, na.rm=T)

# SHD
SHD_fci10pLVdense <- fci_10pLVdense %>% 
  map_depth(2, ~SHD(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci10pLVdense)

# CCI
res_cci10pLVdense <- cci_10pLVdense %>% 
  map_depth(2, ~precision_recall(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_cci10pLVdense <- cci_10pLVdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_cci10pLVdense, na.rm=T)


# SHD
SHD_cci10pLVdense <- cci_10pLVdense %>% 
  map_depth(2, ~SHD(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_cci10pLVdense)


##################
## Running time ##
##################
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
  geom_boxplot(position = "dodge",   outlier.size = 0.8, outlier.alpha = 0.2) + theme_classic() +
  # scale_x_discrete(name ="Condition",
  #                  labels=c("", "5p-sparse", "", "","5p-dense","","", "10p-sparse","","","10p-dense","","","5p-LV","","","10p-LV","")) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(y = " log(ms)", x = "conditions", title = "Algorithm Running Time", subtitle = "Time in milliseconds (ms)") +
  theme(axis.text.x = element_text(face = "bold", angle=40, margin = margin(t = 13)))


