## =============================================================================
## Description 
# 
# This script contains the code for the secondary analysis 
# with varying alpha levels based on the sample size (N).
#
# As is the case with the main simulation study, there are in total 8 models and
# we generate 500 datasets from each model.
## =============================================================================
# The content is as follows.
# 0. Preparation: Source and load necessary functions & packages and 
#    generate data with varying alpha levels based on the sample size (N).
#
# 1. Run algorithms: Apply three algorithms CCD, FCI, and CCI to the simulated 
#    data then estimate PAGs.
#
# 2. Evaluate performance: Compute structural Hamming distance, precision, 
#    recall, and uncertainty rate for each condition.
#
# 3. Organize results: Make neat data frames of resulting values of evaluation 
#    metrics from each algorithm.
#
# 4. Create figures: Plot figures for each evaluation metric comparing the 
#    performance of each algorithm per condition.
## =============================================================================


## =============================================================================
## 0. Preparation
## =============================================================================
# source the simulation study results
source("simulation/01_simulation.R")

# load packages
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(magrittr) # for assigning pipes %<>%

## Specify conditions
# specify the sample sizes
N <- c(50, 150, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000)
# specify replication number
n <- 500
# vary alpha depending on N
alpha <- 1/sqrt(N)
# allow parallel processing
plan(multisession) 


## Generate data
# models without LV
simdat_alphawoLV <- list(B5sparse = B5sparse, B5dense = B5dense, 
                         B10sparse =  B10sparse, B10dense = B10dense) %>% 
  map(~generatesimdat(.x, N)
  )
# 5p models with LV
simdata_alpha5pwLV <- list(B5_lvsparse = B5_lvsparse, B5_lvdense = B5_lvdense) %>% 
  map(~
        generatesimdat(.x, LV = 6, N)
  )
# 10p models with LV
simdata_alpha10pwLV <- list(B10_lvsparse = B10_lvsparse, B10_lvdense = B10_lvdense) %>% 
  map(~
        generatesimdat(.x, LV = c(11, 12), N)
  )
# append all data sets together
simdat_alpha <- append(simdat_alphawoLV, append(simdata_alpha5pwLV, simdata_alpha10pwLV))
simdat_alpha2 <- simdat_alpha %>% bind_rows(.id="id")


## =============================================================================
## 1. Running algorithms
## =============================================================================
## =============
## B5 sparse 
## =============
# run CCD on B5 sparse data
CCDB5sparse <- list()
for(i in 1:length(N)){
  CCDB5sparse[[i]] <- simdat_alpha2 %>% filter(id == "B5sparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
# run FCI on B5 sparse data
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
# run CCI on B5 sparse data
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
# run CCD on B5 dense data
CCDB5dense <- list()
for(i in 1:length(N)){
  CCDB5dense[[i]] <- simdat_alpha2 %>% filter(id == "B5dense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
# run FCI on B5 dense data
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
# run CCI on B5 dense data
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
# run CCD on B10 sparse data
CCDB10sparse <- list()
for(i in 1:length(N)){
  CCDB10sparse[[i]] <- simdat_alpha2 %>% filter(id == "B10sparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
# run FCI on B10 sparse data
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
# run CCI on B10 sparse data
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
# run CCD on B10 dense data
CCDB10dense <- list()
for(i in 1:length(N)){
  CCDB10dense[[i]] <- simdat_alpha2 %>% filter(id == "B10dense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
# run FCI on B10 dense data
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
# run CCI on B10 dense data
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
# run CCD on B5 sparse LV data
CCDB5_LVsparse <- list()
for(i in 1:length(N)){
  CCDB5_LVsparse[[i]] <- simdat_alpha2 %>% filter(id == "B5_lvsparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
# run FCI on B5 sparse LV data
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
# run CCI on B5 sparse LV data
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
# run CCD on B5 dense LV data
CCDB5_LVdense <- list()
for(i in 1:length(N)){
  CCDB5_LVdense[[i]] <- simdat_alpha2 %>% filter(id == "B5_lvdense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
# run FCI on B5 dense LV data
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
# run CCI on B5 dense LV data
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
# run CCD on B10 sparse LV data
CCDB10_LVsparse <- list()
for(i in 1:length(N)){
  CCDB10_LVsparse[[i]] <- simdat_alpha2 %>% filter(id == "B10_lvsparse") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
# run FCI on B10 sparse LV data
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
# run CCI on B10 sparse LV data
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
# run CCD on B10 dense LV data
CCDB10_LVdense <- list()
for(i in 1:length(N)){
  CCDB10_LVdense[[i]] <- simdat_alpha2 %>% filter(id == "B10_lvdense") %>% 
    dplyr::select(paste(N[i])) %>% 
    .[[paste(N[i])]] %>%  
    map(~ccdKP(df=.x, dataType = "continuous", alpha = alpha[i]) %>% 
          CreateAdjMat(., length(.$nodes))
    ) %>% 
    rlang::set_names(., N[i])
}
# run FCI on B10 dense LV data
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
# run CCI on B10 dense LV data
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


## =============================================================================
## 2. Evaluating performance 
## =============================================================================
## =============
## B5 sparse
## =============
## CCD
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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
# Precision & Recall
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



## =============================================================================
## 3. Create neat dataframes
## =============================================================================

## Compute average precision & recall and corresponding sd for each condition
pre_rec2 <- list(
  # put all the results together in a list
  res_ccd5psparse2, res_fci5psparse2, res_cci5psparse2, res_ccd10psparse2, res_fci10psparse2, 
  res_cci10psparse2, res_ccd5pdense2, res_fci5pdense2, res_cci5pdense2, res_ccd10pdense2, 
  res_fci10pdense2, res_cci10pdense2, res_ccd5pLVsparse2, res_fci5pLVsparse2, res_cci5pLVsparse2,
  res_ccd5pLVdense2, res_fci5pLVdense2, res_cci5pLVdense2, res_ccd10pLVsparse2, res_fci10pLVsparse2, res_cci10pLVsparse2, res_ccd10pdense2,  res_fci10pLVdense2, res_cci10pLVdense2
) %>% 
  # transpose df
  map(~ sjmisc::rotate_df(.x) %>%
        # add sample size (N) info
        rename_with(~paste0(.x, "N = ", rep(N, each=8)))  %>%
        # get the mean and sd
        dplyr::summarise(across(everything(.), list(mean = ~mean(., na.rm=T), 
                                                    sd = ~sd(., na.rm=T))))) %>%
  bind_rows() %>% 
  # create columns (conditions)
  mutate(algorithm = rep(c("CCD", "FCI", "CCI"), 8),
         condition = rep(c("5p_sparse", "10p_sparse", "5p_dense", "10p_dense", 
                           "5p_LVsparse", "5p_LVdense", "10p_LVsparse", "10p_LVdense"), each=3),
         netsize = stringr::str_split(condition, "_", simplify=T)[,1],
         latentvar = ifelse(stringr::str_detect(condition, "LV")==TRUE, "with LC", "without LC"),
         densities = stringr::str_remove(stringr::str_split(condition, "_", simplify=T)[,2], "LV")
  ) %>%
  # brings the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric)) %>% 
  # convert it to a long format
  tidyr::pivot_longer(!c(algorithm, condition, netsize, latentvar, densities), 
                      names_to = "metric", values_to = "value") %>% 
  # add sample size column (N) & clean up the column name 
  mutate(N = stringr::str_extract(metric, "(?<=[N =])\\d+"),
         metric = stringr::str_replace_all(metric, "[0-9.]+|[N =]", "")) 


## Compute average uncertainty rate and corresponding sd for each condition
uncertainties2 <- bind_rows(
  # bind all results from each condition
  "CCD_5p-sparse" = uncer_ccd5psparse2, "FCI_5p-sparse" = uncer_fci5psparse2, 
  "CCI_5p-sparse"= uncer_cci5psparse2, "CCD_10p-sparse"= uncer_ccd10psparse2, 
  "FCI_10p-sparse" = uncer_fci10psparse2, "CCI_10p-sparse" = uncer_cci10psparse2, 
  "CCD_5p-dense" = uncer_ccd5pdense2, "FCI_5p-dense"= uncer_fci5pdense2, 
  "CCI_5p-dense" = uncer_cci5pdense2, "CCD_10p-dense"= uncer_ccd10pdense2, 
  "FCI_10p-dense" = uncer_fci10pdense2, "CCI_10p-dense"= uncer_cci10pdense2, 
  "CCD_5p-LVsparse" = uncer_ccd5pLVsparse2, "FCI_5p-LVsparse" = uncer_fci5pLVsparse2, 
  "CCI_5p-LVsparse" = uncer_cci5pLVsparse2, "CCD_10p-LVsparse" = uncer_ccd10pLVsparse2, 
  "FCI_10p-LVsparse" = uncer_fci10pLVsparse2, "CCI_10p-LVsparse" = uncer_cci10pLVsparse2,
  "CCD_5p-LVdense" = uncer_ccd5pLVdense2, "FCI_5p-LVdense" = uncer_fci5pLVdense2, 
  "CCI_5p-LVdense" = uncer_cci5pLVdense2, "CCD_10p-LVdense" = uncer_ccd10pLVdense2, 
  "FCI_10p-LVdense" = uncer_fci10pLVdense2, "CCI_10p-LVdense" = uncer_cci10pLVdense2, .id="id"
) %>% 
  group_by(id) %>% 
  # get the mean and sd
  summarise_all(list(means = mean, sds = sd)) %>%  
  # create columns (conditions)
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2],
         netsize = stringr::str_split(condition, "-", simplify=T)[,1],
         latentvar = ifelse(stringr::str_detect(condition, "LV")==TRUE, "with LC", "without LC"),
         densities = stringr::str_remove(stringr::str_split(condition, "-", simplify=T)[,2], "LV")
  ) %>% 
  # convert it to a long format
  tidyr::pivot_longer(!c(algorithm, condition, id, netsize, latentvar, densities), 
                      names_to = "name", values_to = "value") %>% 
  # add sample size column (N) & clean up the column name 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  
  # bring the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric))


## Compute average SHD values and corresponding sd for each condition
SHDs2 <- bind_rows(
  # bind all results from each condition
  "CCD_5p-sparse" = SHD_ccd5psparse2, "FCI_5p-sparse" = SHD_fci5psparse2, 
  "CCI_5p-sparse" = SHD_cci5psparse2, "CCD_10p-sparse"= SHD_ccd10psparse2, 
  "FCI_10p-sparse" = SHD_fci10psparse2, "CCI_10p-sparse" = SHD_cci10psparse2, 
  "CCD_5p-dense" = SHD_ccd5pdense2, "FCI_5p-dense" = SHD_fci5pdense2, 
  "CCI_5p-dense" =SHD_cci5pdense2, "CCD_10p-dense" = SHD_ccd10pdense2, 
  "FCI_10p-dense" = SHD_fci10pdense2, "CCI_10p-dense" = SHD_cci10pdense2, 
  "CCD_5p-LVsparse" = SHD_ccd5pLVsparse2, "FCI_5p-LVsparse" = SHD_fci5pLVsparse2, 
  "CCI_5p-LVsparse" = SHD_cci5pLVsparse2, "CCD_10p-LVsparse" = SHD_ccd10pLVsparse2, 
  "FCI_10p-LVsparse" = SHD_fci10pLVsparse2, "CCI_10p-LVsparse" = SHD_cci10pLVsparse2, 
  "CCD_5p-LVdense" = SHD_ccd5pLVdense2, "FCI_5p-LVdense" = SHD_fci5pLVdense2, 
  "CCI_5p-LVdense" = SHD_cci5pLVdense2, "CCD_10p-LVdense" = SHD_ccd10pLVdense2, 
  "FCI_10p-LVdense" = SHD_fci10pLVdense2, "CCI_10p-LVdense" = SHD_cci10pLVdense2, .id="id"
) %>% 
  group_by(id) %>% 
  # get the mean and sd
  summarise_all(list(means = mean, sds = sd)) %>%  
  # create columns (conditions)
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2],
         netsize = stringr::str_split(condition, "-", simplify=T)[,1],
         latentvar = ifelse(stringr::str_detect(condition, "LV")==TRUE, "with LC", "without LC"),
         densities = stringr::str_remove(stringr::str_split(condition, "-", simplify=T)[,2], "LV")
  ) %>% 
  # convert it to a long format
  tidyr::pivot_longer(!c(algorithm, condition, id, netsize, latentvar, densities), 
                      names_to = "name", values_to = "value") %>% 
  # add sample size column (N) & clean up the column name 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  
  # bring the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric)) 



## =============================================================================
## 4. Create figures
## =============================================================================

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
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  # create ggplot object
  ggplot(aes(x= as.numeric(N), y=means, group = algorithm, 
             colour = algorithm, fill = algorithm)) +
  # add line plots
  geom_line(aes(group = algorithm)) +
  # add scattered points
  geom_point(size=1) + 
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, 
                  ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="", y="", title = "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create a facet
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels = c("p = 5", "p = 10")) ~ 
                        factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")),  
                      scales = "free_y", switch="y") +
  # manually specify the x-axis break
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  ggtitle("(a) SHD") +
  guides(color = "none", fill = "none")

# save the plot
# ggsave(filename = "results/varyingalpha_shd.pdf", width = 25, height = 10, dpi = 300, units = "cm")



## Precision figure
precision_plot <- pre_rec2 %>% 
  # select only the average precision
  filter(grepl("average_precision", metric)) %>% 
  # convert it to a wide format
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  # create ggplot object
  ggplot(aes(x= as.numeric(N), y=average_precision_mean, group = algorithm, 
             colour = algorithm, fill=algorithm)) +
  # add line plots
  geom_line(aes(group = algorithm)) +
  # add scattered points
  geom_point(size=1) +
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=average_precision_mean+qnorm(0.25)*average_precision_sd, 
                  ymax=average_precision_mean+qnorm(0.75)*average_precision_sd), 
                  alpha=0.15, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create a facet
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ 
                        factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")),  switch="y") +
  # manually specify the x-axis break
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  labs(title = "(b) Precision", x = "", y = "") +
  guides(color = "none", fill = "none")

# save the plot
# ggsave(filename = "results/varyingalpha_prec.pdf", width = 25, height = 10, dpi = 300, units = "cm")



## Recall figure
recall_plot <- pre_rec2 %>% 
  # select only the average recall
  filter(grepl("average_recall", metric)) %>% 
  # convert it to a wide format
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  # create ggplot object
  ggplot(aes(x= as.numeric(N), y=average_recall_mean, group = algorithm, 
             colour = algorithm, fill= algorithm)) +
  # add line plots
  geom_line(aes(group = algorithm)) +
  # add scattered points
  geom_point(size=1) +
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=average_recall_mean+qnorm(0.25)*average_recall_sd, 
                  ymax=average_recall_mean+qnorm(0.75)*average_recall_sd), 
                  alpha=0.15, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create a facet
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ 
                        factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")),  switch="y") +
  # manually specify the x-axis break
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  labs(title = "(c) Recall", x = "", y = "") + 
  guides(color = "none", fill = "none")

# save the plot
# ggsave(filename = "results/varyingalpha_rec.pdf", width = 25, height = 10, dpi = 300, units = "cm")

# combine the plots together
#ggpubr::ggarrange(precision_plot, recall_plot, nrow=2, common.legend = TRUE, legend = "bottom")



## Uncertainty figure
uncertainty_plot <- uncertainties2 %>%
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  # create ggplot object
  ggplot(aes(x= as.numeric(N), y=means, group = algorithm, colour = algorithm, fill = algorithm)) +
  # add line plots
  geom_line(aes(group = algorithm)) +
  # add scattered points
  geom_point(size=1) + 
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="N", y="", title = "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create a facet
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ 
                        factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")),  scales = "free_y", switch="y") +
  # manually specify the x-axis break
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  ggtitle("(d) Uncertainty")

# save the plot
# ggsave(filename = "results/varyingalpha_unc.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")


# combine the plots
ggpubr::ggarrange(shdplot, precision_plot, nrow=2, common.legend = TRUE, legend = "bottom")
# ggsave(filename = "results/varyingalpha_result1.pdf", width = 25, height = 22, dpi = 300, units = "cm")

ggpubr::ggarrange(recall_plot, uncertainty_plot, nrow=2, common.legend = TRUE, legend = "bottom")
# ggsave(filename = "results/varyingalpha_result2.pdf", width = 25, height = 22, dpi = 300, units = "cm")

ggpubr::ggarrange(shdplot, precision_plot, recall_plot, uncertainty_plot, nrow=4, common.legend = TRUE, legend = "bottom")
# ggsave(filename = "results/varyingalpha_result.pdf", width = 25, height = 35, dpi = 300, units = "cm")

