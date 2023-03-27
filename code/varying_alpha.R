## =============================================================================
## Description
#
# This script contains all the code for the secondary analysis 
# with varying alpha level.
# As is the case with the main simulation study, there are in total 8 models and
# we generate 500 datasets from each model.
#
# The content is as follows:
# 0. Preparation: we source and load necessary functions & packages and generate data.
# 1. Run algorithms: we again run three algorithms CCD, FCI, and CCI then estimate PAGs.
# 2. Evaluate performance: we compute structural Hamming distance, precision, recall,
# and uncertainty rate for each condition.
# 3. Organize results: we make neat data frames of resulting values of evaluation
# metrics from each algorithm.
# 4. Create figures: we create figures for each evaluation metric comparing the 
# performance of each algorithm per condition.
## =============================================================================

## ========================
## 0. Preparation
## ========================
# source the simulation study results
source("code/simulation_code.R")
source("code/R/eval_metrics.R")
source("code/R/true_ancestral.R")

# load packages
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(magrittr) # for assigning pipes %<>%

## Data generating
# specify the sample sizes
N <- c(50, 150, 500, 1000, 1500, 2000, 3000, 4000, 5000, 10000)
# specify replication number
n <- 500
# vary alpha depending on N
alpha <- 1/sqrt(N)

# allow parallel processing
plan(multisession) 

generatesimdat <- function(B, N, LV = NULL, seed=123){
  # generate data n times for each N
  simdat <- N %>% future_map(function(z) {
    replicate(n = n,
              expr = if(is.null(LV)) gen_dat(B, N = z) else gen_dat(B, N = z)[,-LV],
              simplify = FALSE)
  }, .options = furrr_options(seed=seed) 
  ) %>% 
    rlang::set_names(.,  N)
  return(simdat)
}


# generate data
# simdat_alpha <- list(B5sparse = B5sparse, B5dense = B5dense, B10sparse =  B10sparse, B10dense = B10dense, B10_lvsparse = B10_lvsparse, B10_lvdense = B10_lvdense) %>% 
#   map(~generatesimdat(.x, N)
#         )


## create all simulated data using random B 
simdat_alphawoLV <- list(B5sparse = B5sparse, B5dense = B5dense, B10sparse =  B10sparse, B10dense = B10dense) %>% 
  map(~generatesimdat(.x, N)
  )

simdata_alpha5pwLV <- list(B5_lvsparse = B5_lvsparse, B5_lvdense = B5_lvdense) %>% 
  map(~
        generatesimdat(.x, LV = 6, N)
  )

simdata_alpha10pwLV <- list(B10_lvsparse = B10_lvsparse, B10_lvdense = B10_lvdense) %>% 
  map(~
        generatesimdat(.x, LV = c(11, 12), N)
  )

simdat_alpha <- append(simdat_alphawoLV, append(simdata_alpha5pwLV, simdata_alpha10pwLV))

simdat_alpha2 <- simdat_alpha %>% bind_rows(.id="id")

## ============================
## 1. Running algorithms
## ============================
## =============
## B5 sparse 
## =============
CCDB5sparse <- list()
for(i in 1:length(N)){
  CCDB5sparse[[i]] <- simdat_alpha2 %>% filter(id == "B5sparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  # can I use sth else not paste? lol 
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
FCIB5sparse <- list()
for(i in 1:length(N)){
  FCIB5sparse[[i]] <- simdat_alpha2 %>% filter(id == "B5sparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest, alpha = alpha[i], 
        doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
    ) %>% 
    rlang::set_names(., N[i])
}
CCIB5sparse <- list()
for(i in 1:length(N)){
  CCIB5sparse[[i]] <- simdat_alpha2 %>% filter(id == "B5sparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha[i], 
               labels = colnames(.x), p = ncol(.x)) %>% .$maag
    ) %>% 
    rlang::set_names(., N[i])
}

## =============
## B5 dense 
## =============
CCDB5dense <- list()
for(i in 1:length(N)){
  CCDB5dense[[i]] <- simdat_alpha2 %>% filter(id == "B5dense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  # can I use sth else not paste? lol 
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
FCIB5dense <- list()
for(i in 1:length(N)){
  FCIB5dense[[i]] <- simdat_alpha2 %>% filter(id == "B5dense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest, alpha = alpha[i], 
              doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
    ) %>% 
    rlang::set_names(., N[i])
}
CCIB5dense <- list()
for(i in 1:length(N)){
  CCIB5dense[[i]] <- simdat_alpha2 %>% filter(id == "B5dense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha[i], 
             labels = colnames(.x), p = ncol(.x)) %>% .$maag
    ) %>% 
    rlang::set_names(., N[i])
}

## =============
## B10 sparse 
## =============
CCDB10sparse <- list()
for(i in 1:length(N)){
  CCDB10sparse[[i]] <- simdat_alpha2 %>% filter(id == "B10sparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  # can I use sth else not paste? lol 
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
FCIB10sparse <- list()
for(i in 1:length(N)){
  FCIB10sparse[[i]] <- simdat_alpha2 %>% filter(id == "B10sparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest, alpha = alpha[i], 
              doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
    ) %>% 
    rlang::set_names(., N[i])
}
CCIB10sparse <- list()
for(i in 1:length(N)){
  CCIB10sparse[[i]] <- simdat_alpha2 %>% filter(id == "B10sparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha[i], 
             labels = colnames(.x), p = ncol(.x)) %>% .$maag
    ) %>% 
    rlang::set_names(., N[i])
}

## =============
## B10 dense 
## =============
CCDB10dense <- list()
for(i in 1:length(N)){
  CCDB10dense[[i]] <- simdat_alpha2 %>% filter(id == "B10dense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  # can I use sth else not paste? lol 
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
FCIB10dense <- list()
for(i in 1:length(N)){
  FCIB10dense[[i]] <- simdat_alpha2 %>% filter(id == "B10dense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest, alpha = alpha[i], 
              doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
    ) %>% 
    rlang::set_names(., N[i])
}
CCIB10dense <- list()
for(i in 1:length(N)){
  CCIB10dense[[i]] <- simdat_alpha2 %>% filter(id == "B10dense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha[i], 
             labels = colnames(.x), p = ncol(.x)) %>% .$maag
    ) %>% 
    rlang::set_names(., N[i])
}


## =============
## B5 sparse LV
## =============
CCDB5_LVsparse <- list()
for(i in 1:length(N)){
  CCDB5_LVsparse[[i]] <- simdat_alpha2 %>% filter(id == "B5_lvsparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  # can I use sth else not paste? lol 
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
FCIB5_LVsparse <- list()
for(i in 1:length(N)){
  FCIB5_LVsparse[[i]] <- simdat_alpha2 %>% filter(id == "B5_lvsparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest, alpha = alpha[i], 
              doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
    ) %>% 
    rlang::set_names(., N[i])
}
CCIB5_LVsparse <- list()
for(i in 1:length(N)){
  CCIB5_LVsparse[[i]] <- simdat_alpha2 %>% filter(id == "B5_lvsparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha[i], 
             labels = colnames(.x), p = ncol(.x)) %>% .$maag
    ) %>% 
    rlang::set_names(., N[i])
}

## =============
## B5 dense LV
## =============
CCDB5_LVdense <- list()
for(i in 1:length(N)){
  CCDB5_LVdense[[i]] <- simdat_alpha2 %>% filter(id == "B5_lvdense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  # can I use sth else not paste? lol 
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
FCIB5_LVdense <- list()
for(i in 1:length(N)){
  FCIB5_LVdense[[i]] <- simdat_alpha2 %>% filter(id == "B5_lvdense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest, alpha = alpha[i], 
              doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
    ) %>% 
    rlang::set_names(., N[i])
}
CCIB5_LVdense <- list()
for(i in 1:length(N)){
  CCIB5_LVdense[[i]] <- simdat_alpha2 %>% filter(id == "B5_lvdense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha[i], 
             labels = colnames(.x), p = ncol(.x)) %>% .$maag
    ) %>% 
    rlang::set_names(., N[i])
}

## =============
## B10 sparse LV
## =============
CCDB10_LVsparse <- list()
for(i in 1:length(N)){
  CCDB10_LVsparse[[i]] <- simdat_alpha2 %>% filter(id == "B10_lvsparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  # can I use sth else not paste? lol 
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
FCIB10_LVsparse <- list()
for(i in 1:length(N)){
  FCIB10_LVsparse[[i]] <- simdat_alpha2 %>% filter(id == "B10_lvsparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest, alpha = alpha[i], 
              doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
    ) %>% 
    rlang::set_names(., N[i])
}
CCIB10_LVsparse <- list()
for(i in 1:length(N)){
  CCIB10_LVsparse[[i]] <- simdat_alpha2 %>% filter(id == "B10_lvsparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha[i], 
             labels = colnames(.x), p = ncol(.x)) %>% .$maag
    ) %>% 
    rlang::set_names(., N[i])
}

## =============
## B10 dense LV
## =============
CCDB10_LVdense <- list()
for(i in 1:length(N)){
  CCDB10_LVdense[[i]] <- simdat_alpha2 %>% filter(id == "B10_lvdense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  # can I use sth else not paste? lol 
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
FCIB10_LVdense <- list()
for(i in 1:length(N)){
  FCIB10_LVdense[[i]] <- simdat_alpha2 %>% filter(id == "B10_lvdense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest, alpha = alpha[i], 
              doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
    ) %>% 
    rlang::set_names(., N[i])
}
CCIB10_LVdense <- list()
for(i in 1:length(N)){
  CCIB10_LVdense[[i]] <- simdat_alpha2 %>% filter(id == "B10_lvdense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>% 
    map(~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha = alpha[i], 
             labels = colnames(.x), p = ncol(.x)) %>% .$maag
    ) %>% 
    rlang::set_names(., N[i])
}


## ============================
## 2. Evaluating performance 
## ============================
## =============
## B5 sparse
## =============
## CCD
res_ccd5psparse2 <- CCDB5sparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_ccd5psparse2 <- CCDB5sparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% do.call("cbind", .) %>% 
  apply(., 2, unlist) %>%  as.data.frame %>% 
  rename_with(~ paste0("N = ", N))
# SHD
SHD_ccd5psparse2 <- CCDB5sparse %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% do.call("cbind", .) %>% 
  apply(., 2, unlist) %>%  as.data.frame %>% rename_with(~ paste0("N = ", N))

## FCI
res_fci5psparse2 <- FCIB5sparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_fci5psparse2 <- FCIB5sparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_fci5psparse2 <- FCIB5sparse %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## CCI
res_cci5psparse2 <- CCIB5sparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_cci5psparse2 <- CCIB5sparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_cci5psparse2 <- CCIB5sparse %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))


## =============
## B5 dense
## =============
## CCD
res_ccd5pdense2 <- CCDB5dense %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_ccd5pdense2 <- CCDB5dense %>% 
  map_depth(2, ~uncertainty(.x)) %>% do.call("cbind", .) %>% 
  apply(., 2, unlist) %>%  as.data.frame %>% 
  rename_with(~ paste0("N = ", N))
# SHD
SHD_ccd5pdense2 <- CCDB5dense %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% do.call("cbind", .) %>% 
  apply(., 2, unlist) %>%  as.data.frame %>% rename_with(~ paste0("N = ", N))

## FCI
res_fci5pdense2 <- FCIB5dense %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_fci5pdense2 <- FCIB5dense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_fci5pdense2 <- FCIB5dense %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## CCI
res_cci5pdense2 <- CCIB5dense %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_cci5pdense2 <- CCIB5dense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_cci5pdense2 <- CCIB5dense %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))


## =============
## B10 sparse
## =============
## CCD
res_ccd10psparse2 <- CCDB10sparse %>% 
  map_depth(2, 
            ~precision_recall(trueag_10psparse, .x)) %>%
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame()
# UNCERTAINTY
uncer_ccd10psparse2 <- CCDB10sparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_ccd10psparse2 <- CCDB10sparse %>% 
  map_depth(2, ~SHD(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## FCI
res_fci10psparse2 <- FCIB10sparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_fci10psparse2 <- FCIB10sparse%>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_fci10psparse2 <- FCIB10sparse %>% 
  map_depth(2, ~SHD(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## CCI
res_cci10psparse2 <- CCIB10sparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_cci10psparse2 <- CCIB10sparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_cci10psparse2 <- CCIB10sparse %>% 
  map_depth(2, ~SHD(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))


## =============
## B10 dense
## =============
## CCD
res_ccd10pdense2  <- CCDB10dense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_ccd10pdense2  <- CCDB10dense  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_ccd10pdense2  <- CCDB10dense %>% 
  map_depth(2, ~SHD(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## FCI
res_fci10pdense2  <- FCIB10dense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_fci10pdense2 <- FCIB10dense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_fci10pdense2  <- FCIB10dense %>% 
  map_depth(2, ~SHD(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## CCI
res_cci10pdense2 <- CCIB10dense %>% 
  map_depth(2, ~precision_recall(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_cci10pdense2 <- CCIB10dense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_cci10pdense2 <- CCIB10dense %>% 
  map_depth(2, ~SHD(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## =============
## B5 sparse LV
## =============
## CCD
res_ccd5pLVsparse2  <- CCDB5_LVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_ccd5pLVsparse2  <- CCDB5_LVsparse  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_ccd5pLVsparse2 <- CCDB5_LVsparse %>% 
  map_depth(2, ~SHD(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## FCI
res_fci5pLVsparse2  <- FCIB5_LVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_fci5pLVsparse2 <- FCIB5_LVsparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_fci5pLVsparse2 <- FCIB5_LVsparse %>% 
  map_depth(2, ~SHD(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## CCI
res_cci5pLVsparse2 <- CCIB5_LVsparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_cci5pLVsparse2 <- CCIB5_LVsparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_cci5pLVsparse2 <- CCIB5_LVsparse %>% 
  map_depth(2, ~SHD(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## =============
## B5 dense LV
## =============
## CCD
res_ccd5pLVdense2  <- CCDB5_LVdense  %>% 
  map_depth(2, ~precision_recall(trueag_5pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_ccd5pLVdense2  <- CCDB5_LVdense  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_ccd5pLVdense2 <- CCDB5_LVdense %>% 
  map_depth(2, ~SHD(trueag_5pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## FCI
res_fci5pLVdense2  <- FCIB5_LVdense  %>% 
  map_depth(2, ~precision_recall(trueag_5pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_fci5pLVdense2 <- FCIB5_LVdense  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_fci5pLVdense2  <- FCIB5_LVdense  %>% 
  map_depth(2, ~SHD(trueag_5pdenseLV , .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## CCI
res_cci5pLVdense2  <- CCIB5_LVdense  %>% 
  map_depth(2, ~precision_recall(trueag_5pdenseLV , .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_cci5pLVdense2  <- CCIB5_LVdense  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_cci5pLVdense2  <- CCIB5_LVdense  %>% 
  map_depth(2, ~SHD(trueag_5pdenseLV , .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## =============
## B10 sparse LV
## =============
## CCD
res_ccd10pLVsparse2   <- CCDB10_LVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_ccd10pLVsparse2  <- CCDB10_LVsparse  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_ccd10pLVsparse2 <- CCDB10_LVsparse %>% 
  map_depth(2, ~SHD(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## FCI
res_fci10pLVsparse2 <- FCIB10_LVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_fci10pLVsparse2 <- FCIB10_LVsparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_fci10pLVsparse2 <- FCIB10_LVsparse %>% 
  map_depth(2, ~SHD(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## CCI 
res_cci10pLVsparse2 <- CCIB10_LVsparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_cci10pLVsparse2 <- CCIB10_LVsparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_cci10pLVsparse2 <- CCIB10_LVsparse %>% 
  map_depth(2, ~SHD(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## =============
## B10 dense LV
## =============
## CCD
res_ccd10pLVdense2   <- CCDB10_LVdense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_ccd10pLVdense2  <- CCDB10_LVdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_ccd10pLVdense2 <- CCDB10_LVdense %>% 
  map_depth(2, ~SHD(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N)) 

## FCI
res_fci10pLVdense2 <- FCIB10_LVdense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_fci10pLVdense2 <- FCIB10_LVdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_fci10pLVdense2 <- FCIB10_LVdense %>% 
  map_depth(2, ~SHD(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

## CCI
res_cci10pLVdense2 <- CCIB10_LVdense %>% 
  map_depth(2, ~precision_recall(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_cci10pLVdense2 <- CCIB10_LVdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# SHD
SHD_cci10pLVdense2 <- CCIB10_LVdense %>% 
  map_depth(2, ~SHD(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))



## ============================
## 3. Create neat dataframes
## ============================

## Compute average precision & recall and corresponding sd for each condition
pre_rec2 <- list(
  # put all the results together in a list
  res_ccd5psparse2, res_fci5psparse2, res_cci5psparse2, res_ccd10psparse2, res_fci10psparse2, res_cci10psparse2, res_ccd5pdense2, res_fci5pdense2, res_cci5pdense2, res_ccd10pdense2, res_fci10pdense2, res_cci10pdense2, res_ccd5pLVsparse2, res_fci5pLVsparse2, res_cci5pLVsparse2, res_ccd5pLVdense2, res_fci5pLVdense2, res_cci5pLVdense2, res_ccd10pLVsparse2, res_fci10pLVsparse2, res_cci10pLVsparse2, res_ccd10pdense2,  res_fci10pLVdense2, res_cci10pLVdense2
) %>% 
  # transpose df
  map(~ sjmisc::rotate_df(.x) %>%
        # add sample size (N) info
        rename_with(~paste0(.x, "N = ", rep(N, each=8)))  %>%
        # think about how to deal with NAs or do I want to define sth. else instead of NAs.
        # na.omit(.x) %>% 
        # get the average and sd
        dplyr::summarise(across(everything(.), list(mean = ~mean(., na.rm=T), sd = ~sd(., na.rm=T))))) %>% 
  bind_rows() %>% 
  mutate(algorithm = rep(c("CCD", "FCI", "CCI"), 8),
         condition = rep(c("5p_sparse", "10p_sparse", "5p_dense", "10p_dense", "5p_LVsparse", "5p_LVdense", "10p_LVsparse", "10p_LVdense"), each=3),
         netsize = stringr::str_split(condition, "_", simplify=T)[,1],
         latentvar = ifelse(stringr::str_detect(condition, "LV")==TRUE, "with LV", "without LV"),
         densities = stringr::str_remove(stringr::str_split(condition, "_", simplify=T)[,2], "LV")
  ) %>%
  # brings the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric)) %>% 
  # convert it to a long format
  tidyr::pivot_longer(!c(algorithm, condition, netsize, latentvar, densities), names_to = "metric", values_to = "value") %>% 
  # Add sample size column (N) & clean up the column name 
  mutate(N = stringr::str_extract(metric, "(?<=[N =])\\d+"),
         metric = stringr::str_replace_all(metric, "[0-9.]+|[N =]", "")) 



## Compute average uncertainty rate and corresponding sd for each condition
uncertainties2 <- bind_rows(
  # bind all results from each condition
  "CCD_5p-sparse" = uncer_ccd5psparse2, "FCI_5p-sparse" = uncer_fci5psparse2, "CCI_5p-sparse"=uncer_cci5psparse2, "CCD_10p-sparse"=uncer_ccd10psparse2, "FCI_10p-sparse" = uncer_fci10psparse2, "CCI_10p-sparse" = uncer_cci10psparse2, "CCD_5p-dense"=uncer_ccd5pdense2, "FCI_5p-dense"=uncer_fci5pdense2, "CCI_5p-dense"=uncer_cci5pdense2, "CCD_10p-dense"=uncer_ccd10pdense2, "FCI_10p-dense"=uncer_fci10pdense2, "CCI_10p-dense"=uncer_cci10pdense2, "CCD_5p-LVsparse"=uncer_ccd5pLVsparse2, "FCI_5p-LVsparse"=uncer_fci5pLVsparse2, "CCI_5p-LVsparse"=uncer_cci5pLVsparse2, "CCD_10p-LVsparse"=uncer_ccd10pLVsparse2, "FCI_10p-LVsparse"=uncer_fci10pLVsparse2, "CCI_10p-LVsparse"=uncer_cci10pLVsparse2,
  "CCD_5p-LVdense"=uncer_ccd5pLVdense2, "FCI_5p-LVdense"=uncer_fci5pLVdense2, "CCI_5p-LVdense"=uncer_cci5pLVdense2, "CCD_10p-LVdense"=uncer_ccd10pLVdense2, "FCI_10p-LVdense"=uncer_fci10pLVdense2, "CCI_10p-LVdense"=uncer_cci10pLVdense2, .id="id"
) %>% 
  group_by(id) %>% 
  # get the average and sd
  summarise_all(list(means = mean, sds = sd)) %>%  
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2],
         netsize = stringr::str_split(condition, "-", simplify=T)[,1],
         latentvar = ifelse(stringr::str_detect(condition, "LV")==TRUE, "with LV", "without LV"),
         densities = stringr::str_remove(stringr::str_split(condition, "-", simplify=T)[,2], "LV")
  ) %>% 
  tidyr::pivot_longer(!c(algorithm, condition, id, netsize, latentvar, densities), names_to = "name", values_to = "value") %>% 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  relocate(where(is.character), .before = where(is.numeric))



## Compute average SHD values and corresponding sd for each condition
SHDs2 <- bind_rows(
  # bind all results from each condition
  "CCD_5p-sparse" = SHD_ccd5psparse2, "FCI_5p-sparse" = SHD_fci5psparse2, "CCI_5p-sparse"=SHD_cci5psparse2, "CCD_10p-sparse"= SHD_ccd10psparse2, "FCI_10p-sparse" = SHD_fci10psparse2, "CCI_10p-sparse" = SHD_cci10psparse2, "CCD_5p-dense"= SHD_ccd5pdense2, "FCI_5p-dense"=SHD_fci5pdense2, "CCI_5p-dense"=SHD_cci5pdense2, "CCD_10p-dense"= SHD_ccd10pdense2, "FCI_10p-dense"=SHD_fci10pdense2, "CCI_10p-dense"=SHD_cci10pdense2, "CCD_5p-LVsparse"=SHD_ccd5pLVsparse2, "FCI_5p-LVsparse"=SHD_fci5pLVsparse2, "CCI_5p-LVsparse"=SHD_cci5pLVsparse2, "CCD_10p-LVsparse"=SHD_ccd10pLVsparse2, "FCI_10p-LVsparse"=SHD_fci10pLVsparse2, "CCI_10p-LVsparse"=SHD_cci10pLVsparse2, "CCD_5p-LVdense"=SHD_ccd5pLVdense2, "FCI_5p-LVdense"=SHD_fci5pLVdense2, "CCI_5p-LVdense"=SHD_cci5pLVdense2, "CCD_10p-LVdense"=SHD_ccd10pLVdense2, "FCI_10p-LVdense"=SHD_fci10pLVdense2, "CCI_10p-LVdense"=SHD_cci10pLVdense2, .id="id"
) %>% 
  group_by(id) %>% 
  # get the average and sd
  summarise_all(list(means = mean, sds = sd)) %>%  
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2],
         netsize = stringr::str_split(condition, "-", simplify=T)[,1],
         latentvar = ifelse(stringr::str_detect(condition, "LV")==TRUE, "with LV", "without LV"),
         densities = stringr::str_remove(stringr::str_split(condition, "-", simplify=T)[,2], "LV")
  ) %>% 
  tidyr::pivot_longer(!c(algorithm, condition, id, netsize, latentvar, densities), names_to = "name", values_to = "value") %>% 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  relocate(where(is.character), .before = where(is.numeric)) 



## ============================
## 4. Create figures
## ============================

## specify the common figure theme
MyTheme <-  theme(plot.title = element_text(face = "bold", family = "Palatino", size = 15, hjust=0.5),
                  plot.subtitle = element_text(face = "italic", family = "Palatino", size = 15, hjust=0.5),
                  axis.text=element_text(face = "bold",family = "Palatino", size = 11),
                  axis.text.x = element_text(angle = 45, hjust = 1.2, vjust =1.2),
                  axis.title = element_text(face = "bold",family = "Palatino", size = 12),
                  legend.text = element_text(face = "bold", family = "Palatino", size = 12),
                  legend.position="bottom",
                  strip.text = element_text(face="bold", size=13, family = "Palatino"),
                  strip.background = element_rect(fill="#f0f0f0", linetype = "solid", color="gray"),
                  strip.placement = "outside",
                  panel.border = element_rect(color = "#DCDCDC", fill = NA),
                  panel.spacing.y = unit(4, "mm")
)


## SHD figure
shdplot <- SHDs2 %>%
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000")), y=means, group = algorithm, colour = algorithm, fill = algorithm)) +
  # add line plots
  geom_line(aes(group = algorithm)) +
  # add scattered points
  geom_point(size=1) + 
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by 3)
  geom_ribbon(aes(ymin=means-qnorm(0.975)*sds/sqrt(as.numeric(N))*3, ymax=means+qnorm(0.975)*sds/sqrt(as.numeric(N))*3), alpha=0.2, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="", y="", title = "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create a facet
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  scales = "free_y", switch="y") +
  ggtitle("(a) SHD") +
  # remove the x-axis texts & ticks
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


## Precision figure
precision_plot <- pre_rec2 %>% 
  filter(grepl("average_precision", metric)) %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000")), y=average_precision_mean, group = algorithm, colour = algorithm, fill=algorithm)) +
  # add line plots
  geom_line(aes(group = algorithm)) +
  # add scattered points
  geom_point(size=1) +
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by 2)
  geom_ribbon(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))*2, ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))*2), alpha=0.2, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create a facet
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  switch="y") +
  #labs(title = "Precision", x = "N", y = "")
  labs(title = "(b) Precision", x = "", y = "") +
  # remove the x-axis texts & ticks
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

## Recall figure
recall_plot <- pre_rec2 %>% 
  filter(grepl("average_recall", metric)) %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000")), y=average_recall_mean, group = algorithm, colour = algorithm, fill= algorithm)) +
  # add line plots
  geom_line(aes(group = algorithm)) +
  # add scattered points
  geom_point(size=1) +
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by 2)
  geom_ribbon(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))*2, ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))*2), alpha=0.2, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create a facet
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  switch="y") +
  # labs(title = "Recall", x = "N", y = "")
  labs(title = "(c) Recall", x = "", y = "") +
  # remove the x-axis texts & ticks
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# combine the plots together
#ggpubr::ggarrange(precision_plot, recall_plot, nrow=2, common.legend = TRUE, legend = "bottom")

## Uncertainty figure
uncertainty_plot <- uncertainties2 %>%
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000")), y=means, group = algorithm, colour = algorithm, fill = algorithm)) +
  # add line plots
  geom_line(aes(group = algorithm)) +
  # add scattered points
  geom_point(size=1) + 
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by )
  geom_ribbon(aes(ymin=means-qnorm(0.975)*sds/sqrt(as.numeric(N))*2, ymax=means+qnorm(0.975)*sds/sqrt(as.numeric(N))*2), alpha=0.2, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="N", y="", title = "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create a facet
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  scales = "free_y", switch="y") +
  ggtitle("(d) Uncertainty")


ggpubr::ggarrange(shdplot, precision_plot, nrow=2, common.legend = TRUE, legend = "bottom")
# ggsave(filename = "results/varyingalpha_result1.pdf", width = 25, height = 20, dpi = 300, units = "cm")

ggpubr::ggarrange(recall_plot, uncertainty_plot, nrow=2, common.legend = TRUE, legend = "bottom")

# ggsave(filename = "results/varyingalpha_result2.pdf", width = 25, height = 20, dpi = 300, units = "cm")

ggpubr::ggarrange(shdplot, precision_plot, recall_plot, uncertainty_plot, nrow=4, common.legend = TRUE, legend = "bottom")
# ggsave(filename = "results/varyingalpha_result.pdf", width = 25, height = 35, dpi = 300, units = "cm")
