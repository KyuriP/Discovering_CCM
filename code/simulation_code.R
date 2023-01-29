## load functions & packages
source("code/R/CCD_fnc.R")
source("code/R/plot_fnc.R")
source("code/R/data_generating_fnc.R")
source("code/eval_metrics.R")

library(magrittr)
library(purrr)
library(furrr)
library(dplyr)
library(qgraph)

## slightly modified CCI package
install_github("KyuriP/CCI_KP")
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

## True PAG
truepag_5psparse <- matrix(c(
       0, 2, 2, 0, 0,
       3, 0, 3, 3, 3,
       3, 3, 0, 3, 0,
       0, 3, 3, 0, 3,
       0, 2, 0, 2, 0), 5, 5, byrow = TRUE)
dimnames(truepag_5psparse) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
plotAG(truepag_5psparse)

## Run CCD algorithm
ccd_5psparse <- simdata_5psparse %>% 
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
      ) 
mat_5psparse <- ccd_5psparse %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# pag_ccd5psparse <- map2(ccd_5psparse, mat_5psparse, plotPAG)


## Run FCI algorithm
fci_5psparse <- simdata_5psparse %>% 
  map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
           alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
      )

# pag_fci5psparse <- fci_5psparse %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
cci_5psparse <- simdata_5psparse %>% 
  map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
      )

# pag_cci5psparse <- cci_5psparse %>% 
#   map(~plotAG(.x))



## evaluation 
# CCD
res_ccd5psparse <- mat_5psparse %>% 
  map_depth(2, ~precision_recall(truepag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_ccd5psparse <- mat_5psparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% do.call("cbind", .) %>% apply(., 2, unlist) %>%  as.data.frame %>% 
  rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd5psparse, na.rm=T)

# SHD
SHD_ccd5psparse <- mat_5psparse %>% 
  map_depth(2, ~SHD(truepag_5psparse, .x)) %>% do.call("cbind", .) %>% apply(., 2, unlist) %>%  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd5psparse)

# FCI
res_fci5psparse <- fci_5psparse %>% 
  map_depth(2, ~precision_recall(truepag_5psparse, .x)) %>% 
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
  map_depth(2, ~SHD(truepag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci5psparse)

# CCI
res_cci5psparse <- cci_5psparse %>% 
  map_depth(2, ~precision_recall(truepag_5psparse, .x)) %>% 
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
  map_depth(2, ~SHD(truepag_5psparse, .x)) %>% 
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


## True PAG
truepag_5pdense <- matrix(c(
  0, 2, 2, 0, 2,
  3, 0, 3, 3, 3,
  3, 3, 0, 3, 0,
  0, 3, 3, 0, 3,
  3, 2, 0, 2, 0), 5, 5, byrow = TRUE)
dimnames(truepag_5pdense) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
plotAG(truepag_5pdense)


## Run CCD algorithm
ccd_5pdense <- simdata_5pdense %>% 
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
  ) 
mat_5pdense <- ccd_5pdense %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# pag_ccd5pdense <- map2(ccd_5pdense, mat_5pdense, plotPAG)


## Run FCI algorithm
fci_5pdense <- simdata_5pdense %>% 
  map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
           alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # extract amat
  )

# pag_fci5pdense <- fci_5pdense %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
cci_5pdense <- simdata_5pdense %>% 
  map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
  )

# pag_cci5pdense <- cci_5pdense %>% 
#   map(~plotAG(.x))


## evaluation
# CCD
res_ccd5pdense <- mat_5pdense %>% 
  map_depth(2, ~precision_recall(truepag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame()

# UNCERTAINTY
uncer_ccd5pdense <- mat_5pdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_ccd5pdense, na.rm=T)

# SHD
SHD_ccd5pdense <- mat_5pdense %>% 
  map_depth(2, ~SHD(truepag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(SHD_ccd5pdense)

# FCI
res_fci5pdense <- fci_5pdense %>% 
  map_depth(2, ~precision_recall(truepag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci5pdense <- fci_5pdense%>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_fci5pdense, na.rm=T)

# SHD
SHD_fci5pdense <- fci_5pdense %>% 
  map_depth(2, ~SHD(truepag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(SHD_fci5pdense)

# CCI
res_cci5pdense <- cci_5pdense %>% 
  map_depth(2, ~precision_recall(truepag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY 
uncer_cci5pdense <- cci_5pdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_cci5pdense)

# SHD
SHD_cci5pdense <- cci_5pdense %>% 
  map_depth(2, ~SHD(truepag_5pdense, .x)) %>% 
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

## True PAG
truepag_10psparse <- matrix(c(
         0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 
         3, 3, 0, 2, 0, 0, 0, 0, 0, 0, 
         0, 0, 3, 0, 3, 3, 0, 0, 0, 0, 
         0, 0, 0, 3, 0, 3, 0, 0, 0, 0, 
         0, 0, 0, 3, 3, 0, 3, 0, 0, 0, 
         0, 0, 0, 0, 0, 2, 0, 3, 3, 0, 
         0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 
         0, 0, 0, 0, 0, 0, 2, 0, 0, 3, 
         0, 0, 0, 0, 0, 0, 0, 0, 2, 0), 10, 10, byrow = T)

dimnames(truepag_10psparse) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(truepag_10psparse)


## Run CCD algorithm
ccd_10psparse <- simdata_10psparse %>% 
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
  ) 
mat_10psparse <- ccd_10psparse %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# pag_ccd10psparse <- map2(ccd_10psparse, mat_10psparse, plotPAG)


## Run FCI algorithm
fci_10psparse <- simdata_10psparse %>% 
  map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
           alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
  )

# pag_fci10psparse <- fci_10psparse %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
cci_10psparse <- simdata_10psparse %>% 
  map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
  )

# pag_cci10psparse <- cci_10psparse %>% 
#   map(~plotAG(.x))



## evaluation
# CCD
res_ccd10psparse <- mat_10psparse %>% 
  map_depth(2, 
  ~precision_recall(truepag_10psparse, .x)) %>%
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
  map_depth(2, ~SHD(truepag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd10psparse)

# FCI
res_fci10psparse <- fci_10psparse %>% 
  map_depth(2, ~precision_recall(truepag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci10psparse <- fci_10psparse%>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_fci10psparse, na.rm=T)

# SHD
SHD_fci10psparse <- fci_10psparse %>% 
  map_depth(2, ~SHD(truepag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(SHD_fci10psparse)

# CCI
res_cci10psparse <- cci_10psparse %>% 
  map_depth(2, ~precision_recall(truepag_10psparse, .x)) %>% 
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
  map_depth(2, ~SHD(truepag_10psparse, .x)) %>% 
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

## True PAG
truepag_10pdense <- matrix(c(
  0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 
  3, 3, 0, 2, 0, 2, 0, 0, 0, 0, 
  0, 3, 3, 0, 3, 3, 0, 0, 0, 0, 
  0, 3, 0, 3, 0, 3, 3, 0, 0, 0, 
  0, 0, 0, 3, 3, 0, 3, 0, 0, 0, 
  0, 0, 0, 0, 0, 2, 0, 3, 3, 3, 
  0, 0, 0, 0, 0, 0, 2, 0, 0, 3, 
  0, 0, 0, 0, 0, 0, 2, 0, 0, 3, 
  0, 0, 0, 0, 0, 0, 2, 2, 2, 0), 10, 10, byrow = T)

dimnames(truepag_10pdense) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(truepag_10pdense)


## Run CCD algorithm
ccd_10pdense <- simdata_10pdense  %>% 
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
  ) 
mat_10pdense  <- ccd_10pdense  %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# pag_ccd10pdense  <- map2(ccd_10pdense , mat_10pdense , plotPAG)


## Run FCI algorithm
fci_10pdense <- simdata_10pdense  %>% 
  map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
           alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
  )

# pag_fci10pdense <- fci_10pdense  %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
cci_10pdense  <- simdata_10pdense %>% 
  map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
  )

# pag_cci10pdense  <- cci_10pdense  %>% 
#   map(~plotAG(.x))



## evaluation
# CCD
res_ccd10pdense  <- mat_10pdense  %>% 
  map_depth(2, ~precision_recall(truepag_10pdense, .x)) %>% 
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
  map_depth(2, ~SHD(truepag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
mean(SHD_ccd10pdense)

# FCI
res_fci10pdense  <- fci_10pdense  %>% 
  map_depth(2, ~precision_recall(truepag_10pdense, .x)) %>% 
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
  map_depth(2, ~SHD(truepag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci10pdense)

# CCI
res_cci10pdense <- cci_10pdense %>% 
  map_depth(2, ~precision_recall(truepag_10pdense, .x)) %>% 
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
  map_depth(2, ~SHD(truepag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# averae SHD
colMeans(SHD_cci10pdense)


## ====================
## 5p with LV 
## ====================
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

# generate data (sample size as specified above)
simdata_5pLV2 <- N %>% future_map(function(z) {
  replicate(n=1000,
            expr = gen_dat(B5_lv2, N = z)[,-6],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


# # generate data and exclude the LV
# data5pLV2 <- gen_dat(B5_lv2, N =1e6, seed = 12345)[,-6]


## True PAG
truepag_5pLV2 <- matrix(c(
  0, 2, 2, 0, 0,
  2, 0, 3, 3, 3,
  2, 3, 0, 3, 0,
  0, 3, 3, 0, 3,
  0, 2, 0, 2, 0), 5, 5, byrow = TRUE)
dimnames(truepag_5pLV2) <- list(paste("X", 1:5, sep=""), paste("X", 1:5, sep=""))
plotAG(truepag_5pLV2)


## Run CCD algorithm
ccd_5pLV2 <- simdata_5pLV2  %>% 
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
  ) 
mat_5pLV2  <- ccd_5pLV2  %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# pag_ccd5pLV2  <- map2(ccd_5pLV2 , mat_5pLV2 , plotPAG)


## Run FCI algorithm
fci_5pLV2 <- simdata_5pLV2  %>% 
  map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
           alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
  )

# pag_fci5pLV2 <- fci_5pLV2  %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
cci_5pLV2  <- simdata_5pLV2 %>% 
  map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
  )

# pag_cci5pLV2  <- cci_5pLV2  %>% 
#   map(~plotAG(.x))



## evaluation
# CCD
res_ccd5pLV2  <- mat_5pLV2  %>% 
  map_depth(2, ~precision_recall(truepag_5pLV2, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_ccd5pLV2  <- mat_5pLV2  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd5pLV2 , na.rm=T)

# SHD
SHD_ccd5pLV2 <- mat_5pLV2 %>% 
  map_depth(2, ~SHD(truepag_5pLV2, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd5pLV2)

# FCI
res_fci5pLV2  <- fci_5pLV2  %>% 
  map_depth(2, ~precision_recall(truepag_5pLV2, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci5pLV2 <- fci_5pLV2 %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_fci5pLV2, na.rm=T)

# SHD
SHD_fci5pLV2 <- fci_5pLV2 %>% 
  map_depth(2, ~SHD(truepag_5pLV2, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci5pLV2)

# CCI
res_cci5pLV2 <- cci_5pLV2 %>% 
  map_depth(2, ~precision_recall(truepag_5pLV2, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_cci5pLV2 <- cci_5pLV2 %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average unceratinty
colMeans(uncer_cci5pLV2, na.rm=T)

# SHD
SHD_cci5pLV2 <- cci_5pLV2 %>% 
  map_depth(2, ~SHD(truepag_5pLV2, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_cci5pLV2)



## ====================
## 10p with LV
## ====================
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

# generate data (sample size as specified above)
simdata_10pLV <- N %>% future_map(function(z) {
  replicate(n=1000,
            expr = gen_dat(B10_lv, N = z)[,-11],  
            simplify = FALSE)
}, .options = furrr_options(seed=123))

# # generate data and exclude the LV
# data10pLV <- gen_dat(B10_lv, N =1e6, seed = 123)[,-11]
# 
# ## Estimate GGM
# ggm10pLV <- qgraph(cor(data10pLV), layout = layout10LV, theme="colorblind", graph = "pcor")
# 
# layout(t(1:3))

## True PAG
truepag_10pLV <- matrix(c(
  0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 
  3, 3, 0, 2, 0, 2, 0, 0, 0, 0, 
  0, 0, 3, 0, 3, 3, 0, 0, 0, 0, 
  0, 0, 0, 3, 0, 3, 3, 0, 0, 0, 
  0, 0, 3, 3, 3, 0, 3, 0, 0, 0, 
  0, 0, 0, 0, 2, 2, 0, 3, 3, 0, 
  0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 
  0, 0, 0, 0, 0, 0, 2, 0, 0, 3, 
  0, 0, 0, 0, 0, 0, 0, 2, 2, 0), 10, 10, byrow = T)

dimnames(truepag_10pLV) <- list(paste("X", 1:10, sep=""), paste("X", 1:10, sep=""))
plotAG(truepag_10pLV)


## Run CCD algorithm
ccd_10pLV  <- simdata_10pLV   %>% 
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = 0.05)
  ) 
mat_10pLV   <- ccd_10pLV %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))

# pag_ccd10pLV <- map2(ccd_10pLV, mat_10pLV  , plotPAG)


## Run FCI algorithm
fci_10pLV  <- simdata_10pLV   %>% 
  map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
           alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # exxtract amat
  )

# pag_fci10pLV  <- fci_10pLV   %>% 
#   map(~plotAG(.x))

## Run CCI algorithm
cci_10pLV  <- simdata_10pLV  %>% 
  map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=0.05, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
  )

# pag_cci10pLV   <- cci_10pLV %>% 
#   map(~plotAG(.x))



## evaluation
# CCD
res_ccd10pLV   <- mat_10pLV  %>% 
  map_depth(2, ~precision_recall(truepag_10pLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_ccd10pLV  <- mat_10pLV  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd10pLV , na.rm=T)

# SHD
SHD_ccd10pLV <- mat_10pLV %>% 
  map_depth(2, ~SHD(truepag_10pLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd10pLV)

# FCI
res_fci10pLV <- fci_10pLV  %>% 
  map_depth(2, ~precision_recall(truepag_10pLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci10pLV <- fci_10pLV %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_fci10pLV, na.rm=T)

# SHD
SHD_fci10pLV <- fci_10pLV %>% 
  map_depth(2, ~SHD(truepag_10pLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_fci10pLV)

# CCI
res_cci10pLV <- cci_10pLV %>% 
  map_depth(2, ~precision_recall(truepag_10pLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_cci10pLV <- cci_10pLV %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_cci10pLV, na.rm=T)

# SHD
SHD_cci10pLV<- cci_10pLV %>% 
  map_depth(2, ~SHD(truepag_10pLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_cci10pLV)


##################
## Running time ##
##################
#remotes::install_github("joshuaulrich/microbenchmark")
library(microbenchmark)

times <- microbenchmark(
  ccd_5psparse = ccdKP(df=simdata_5psparse[[1]], dataType = "continuous", alpha = 0.05),
  fci_5psparse = fci(list(C = cor(simdata_5psparse[[1]]), n = 1e3),indepTest=gaussCItest, alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(simdata_5psparse[[1]])),
  cci_5psparse = cci(list(C = cor(simdata_5psparse[[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_5psparse[[1]])),
  
  ccd_5pdense = ccdKP(df=simdata_5pdense[[1]], dataType = "continuous", alpha = 0.05),
  fci_5pdense = fci(list(C = cor(simdata_5pdense[[1]]), n = 1e3),indepTest=gaussCItest, alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(simdata_5pdense[[1]])),
  cci_5pdense = cci(list(C = cor(simdata_5pdense[[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_5pdense[[1]])),
  
  ccd_10psparse = ccdKP(df=simdata_10psparse[[1]], dataType = "continuous", alpha = 0.05),
  fci_10psparse = fci(list(C = cor(simdata_10psparse[[1]]), n = 1e3),indepTest=gaussCItest, alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(simdata_10psparse[[1]])),
  cci_10psparse = cci(list(C = cor(simdata_10psparse[[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_10psparse[[1]])),
  
  ccd_10pdense = ccdKP(df=simdata_10pdense[[1]], dataType = "continuous", alpha = 0.05),
  fci_10pdense = fci(list(C = cor(simdata_10pdense[[1]]), n = 1e3),indepTest=gaussCItest, alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(simdata_10pdense[[1]])),
  cci_10pdense = cci(list(C = cor(simdata_10pdense[[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_10pdense[[1]])),
  
  ccd_5pLV = ccdKP(df=simdata_5pLV2[[1]], dataType = "continuous", alpha = 0.05),
  fci_5pLV = fci(list(C = cor(simdata_5pLV2[[1]]), n = 1e3),indepTest=gaussCItest,
                  alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(simdata_5pLV2[[1]])),
  cci_5pLV = cci(list(C = cor(simdata_5pLV2[[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_5pLV2[[1]])),
  
  ccd_10pLV = ccdKP(df=simdata_10pLV[[1]], dataType = "continuous", alpha = 0.05),
  fci_10pLV = fci(list(C = cor(simdata_10pLV[[1]]), n = 1e3),indepTest=gaussCItest,
                alpha = 0.05, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(simdata_10pLV[[1]])),
  cci_10pLV = cci(list(C = cor(simdata_10pLV[[1]]), n = 1e3), gaussCItest, alpha=0.05, p=ncol(simdata_10pLV[[1]]))
)


times <- times %>% 
  mutate(algorithm = substr(expr, 1, 3),
         condition = stringr::str_split(expr, "_", simplify=T)[,2])

## plot the results
times %>%  
  ggplot(aes(x=factor(condition, levels= c("5psparse", "5pdense", "10psparse", "10pdense", "5pLV", "10pLV")), y = log(time), col= factor(algorithm))) +
  geom_boxplot(position = "dodge",   outlier.size = 0.8, outlier.alpha = 0.2) + theme_classic() +
  # scale_x_discrete(name ="Condition", 
  #                  labels=c("", "5p-sparse", "", "","5p-dense","","", "10p-sparse","","","10p-dense","","","5p-LV","","","10p-LV","")) + 
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(y = " log(ms)", x = "conditions", title = "Algorithm Running Time", subtitle = "Time in milliseconds (ms)") + 
  theme(axis.text.x = element_text(face = "bold", angle=40, margin = margin(t = 13)))
  
  
  