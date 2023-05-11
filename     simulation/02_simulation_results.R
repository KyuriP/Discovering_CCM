## ============================================================================
## Description
#
# This script evaluates the performance of algorithms using structural Hamming 
# distance (SHD), precision and recall.
# Then, it creates the overall result figures (Figure 14 & Figure 15 in paper).
#
# (i) First part of this script concerns creating a neat dataframe of 
# each evaluation metrics (precision, recall, and uncertainty rate)
# including their mean and standard deviation (SD) values.
#
# (ii) Second part of this script concerns creating the figures.
## =============================================================================
## The script is organized as follows.
# 0. Preparation: Source and load necessary functions & packages.
#
# 1. Evaluate performance: Compute structural Hamming distance (SHD), precision, 
#           recall, and uncertainty rate for each algorithm under each condition 
#           (8 conditions in total).
#
# 2. Organize results: Make neat data frames of resulting values of 
#           evaluation metrics from each algorithm.
#
# 3. Create figures: Plot the figures for each evaluation metric comparing the 
#           performance of each algorithm per condition.
## ============================================================================



## =============================================================================
## 0. Preparation
## =============================================================================

# source the simulation study results
source("    simulation/01_simulation.R")

# load packages
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggh4x)



## =============================================================================
## 1. Evaluating performance
## =============================================================================
## ================
## 5P sparse 
## ================
## CCD
res_ccd5psparse <- mat_5psparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>% as.data.frame() 

# UNCERTAINTY
uncer_ccd5psparse <- mat_5psparse %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% 
  apply(., 2, unlist) %>%  
  as.data.frame %>% 
  rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_ccd5psparse, na.rm=T)

# SHD
SHD_ccd5psparse <- mat_5psparse %>% 
  map_depth(2, ~SHD(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% 
  apply(., 2, unlist) %>%  
  as.data.frame %>% 
  rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_ccd5psparse)

## FCI
res_fci5psparse <- fci_5psparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## CCI
res_cci5psparse <- cci_5psparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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


## ================
## 5P dense 
## ================
## CCD
res_ccd5pdense <- mat_5pdense %>% 
  map_depth(2, ~precision_recall(trueag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame()
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

## FCI
res_fci5pdense <- fci_5pdense %>% 
  map_depth(2, ~precision_recall(trueag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## CCI
res_cci5pdense <- cci_5pdense %>% 
  map_depth(2, ~precision_recall(trueag_5pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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


## ================
## 10P sparse 
## ================
## CCD
res_ccd10psparse <- mat_10psparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparse, .x)) %>%
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

## FCI
res_fci10psparse <- fci_10psparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## CCI
res_cci10psparse <- cci_10psparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparse, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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


## ================
## 5P dense 
## ================
## CCD
res_ccd10pdense  <- mat_10pdense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## FCI
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

## CCI
res_cci10pdense <- cci_10pdense %>% 
  map_depth(2, ~precision_recall(trueag_10pdense, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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
# average SHD
colMeans(SHD_cci10pdense)


## ================
## 5P sparse w/ LV 
## ================
## CCD
res_ccd5pLVsparse  <- mat_5pLVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## FCI
res_fci5pLVsparse  <- fci_5pLVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## CCI
res_cci5pLVsparse <- cci_5pLVsparse %>% 
  map_depth(2, ~precision_recall(trueag_5psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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


## ================
## 5P dense w/ LV 
## ================
## CCD
res_ccd5pLVdense  <- mat_5pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_5pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## FCI
res_fci5pLVdense  <- fci_5pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_5pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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
colMeans(SHD_fci5pLVdense)

## CCI
res_cci5pLVdense  <- cci_5pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_5pdenseLV , .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
# UNCERTAINTY
uncer_cci5pLVdense  <- cci_5pLVdense  %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_cci5pLVdense , na.rm=T)
# SHD
SHD_cci5pLVdense  <- cci_5pLVdense  %>% 
  map_depth(2, ~SHD(trueag_5pdenseLV , .x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average SHD
colMeans(SHD_cci5pLVdense )


## ================
## 10P sparse w/ LV 
## ================
## CCD
res_ccd10pLVsparse   <- mat_10pLVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## FCI
res_fci10pLVsparse <- fci_10pLVsparse  %>% 
  map_depth(2, ~precision_recall(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## CCI 
res_cci10pLVsparse <- cci_10pLVsparse %>% 
  map_depth(2, ~precision_recall(trueag_10psparseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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


## ================
## 10P dense w/ LV 
## ================
## CCD
res_ccd10pLVdense   <- mat_10pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## FCI
res_fci10pLVdense <- fci_10pLVdense  %>% 
  map_depth(2, ~precision_recall(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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

## CCI
res_cci10pLVdense <- cci_10pLVdense %>% 
  map_depth(2, ~precision_recall(trueag_10pdenseLV, .x)) %>% 
  do.call("cbind", .) %>% t() %>%  
  apply(., 2, unlist) %>%  as.data.frame() 
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



## =============================================================================
## 2. Create neat dataframe
## =============================================================================

## Compute average precision & recall and corresponding sd for each condition
pre_rec <- list(
  # put all the results together in a list
  res_ccd5psparse, res_fci5psparse, res_cci5psparse, res_ccd10psparse, 
  res_fci10psparse, res_cci10psparse, res_ccd5pdense, res_fci5pdense, 
  res_cci5pdense, res_ccd10pdense, res_fci10pdense, res_cci10pdense, 
  res_ccd5pLVsparse, res_fci5pLVsparse, res_cci5pLVsparse, res_ccd5pLVdense, 
  res_fci5pLVdense, res_cci5pLVdense, res_ccd10pLVsparse, res_fci10pLVsparse,
  res_cci10pLVsparse, res_ccd10pdense,  res_fci10pLVdense, res_cci10pLVdense
) %>% 
  # transpose df
  map(~ sjmisc::rotate_df(.x) %>%
        # add sample size (N) info
        rename_with(~paste0(.x, "N = ", rep(N, each=8)))  %>%
        # get the average and sd
        dplyr::summarise(across(everything(.), list(mean = ~mean(., na.rm=T), 
                                                    sd = ~sd(., na.rm=T))))) %>% 
  bind_rows() %>% 
  # create columns 
  mutate(algorithm = rep(c("CCD", "FCI", "CCI"), 8),
         condition = rep(c("5p_sparse", "10p_sparse", "5p_dense", "10p_dense", 
                           "5p_LVsparse", "5p_LVdense", "10p_LVsparse", 
                           "10p_LVdense"), each=3),
         netsize = stringr::str_split(condition, "_", simplify=T)[,1],
         latentvar = ifelse(stringr::str_detect(condition, "LV")==TRUE, 
                            "with LC", "without LC"),
         densities = stringr::str_remove(
           stringr::str_split(condition, "_", simplify=T)[,2], "LV")
  ) %>%
  # bring the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric)) %>% 
  # convert it to a long format
  tidyr::pivot_longer(!c(algorithm, condition, netsize, latentvar, densities), 
                      names_to = "metric", values_to = "value") %>% 
  # add sample size column (N) & clean up the column name 
  mutate(N = stringr::str_extract(metric, "(?<=[N =])\\d+"),
         metric = stringr::str_replace_all(metric, "[0-9.]+|[N =]", "")) 



## Compute average uncertainty rate and corresponding sd for each condition
uncertainties <- bind_rows(
  # bind all results from each condition
  "CCD_5p-sparse" = uncer_ccd5psparse, "FCI_5p-sparse" = uncer_fci5psparse, 
  "CCI_5p-sparse" = uncer_cci5psparse, "CCD_10p-sparse"= uncer_ccd10psparse, 
  "FCI_10p-sparse" = uncer_fci10psparse, "CCI_10p-sparse" = uncer_cci10psparse, 
  "CCD_5p-dense" = uncer_ccd5pdense, "FCI_5p-dense" = uncer_fci5pdense, 
  "CCI_5p-dense" = uncer_cci5pdense, "CCD_10p-dense" = uncer_ccd10pdense, 
  "FCI_10p-dense" = uncer_fci10pdense, "CCI_10p-dense" = uncer_cci10pdense, 
  "CCD_5p-LVsparse" = uncer_ccd5pLVsparse, "FCI_5p-LVsparse" = uncer_fci5pLVsparse, 
  "CCI_5p-LVsparse" = uncer_cci5pLVsparse, "CCD_10p-LVsparse" = uncer_ccd10pLVsparse, 
  "FCI_10p-LVsparse" = uncer_fci10pLVsparse, "CCI_10p-LVsparse" = uncer_cci10pLVsparse,
  "CCD_5p-LVdense" = uncer_ccd5pLVdense, "FCI_5p-LVdense" = uncer_fci5pLVdense, 
  "CCI_5p-LVdense" = uncer_cci5pLVdense, "CCD_10p-LVdense" = uncer_ccd10pLVdense, 
  "FCI_10p-LVdense" = uncer_fci10pLVdense, "CCI_10p-LVdense" = uncer_cci10pLVdense, 
  .id = "id"
) %>% 
  group_by(id) %>% 
  # get the average and sd
  summarise_all(list(means = mean, sds = sd)) %>%  
  # create columns
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2],
         netsize = stringr::str_split(condition, "-", simplify=T)[,1],
         latentvar = ifelse(stringr::str_detect(condition, "LV")==TRUE, 
                            "with LC", "without LC"),
         densities = stringr::str_remove(
           stringr::str_split(condition, "-", simplify=T)[,2], "LV")
  ) %>% 
  # convert it to a long format
  tidyr::pivot_longer(!c(algorithm, condition, id, netsize, latentvar, densities), 
                      names_to = "name", values_to = "value") %>% 
  # add sample size column (N) & clean up the column name 
  mutate(N = stringr::str_extract(
    stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  
  # bring the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric))



## Compute average SHD values and corresponding sd for each condition
SHDs <- bind_rows(
  # bind all results from each condition
  "CCD_5p-sparse" = SHD_ccd5psparse, "FCI_5p-sparse" = SHD_fci5psparse, 
  "CCI_5p-sparse" = SHD_cci5psparse, "CCD_10p-sparse"= SHD_ccd10psparse, 
  "FCI_10p-sparse" = SHD_fci10psparse, "CCI_10p-sparse" = SHD_cci10psparse, 
  "CCD_5p-dense" = SHD_ccd5pdense, "FCI_5p-dense" = SHD_fci5pdense, 
  "CCI_5p-dense" = SHD_cci5pdense, "CCD_10p-dense" = SHD_ccd10pdense, 
  "FCI_10p-dense" = SHD_fci10pdense, "CCI_10p-dense" = SHD_cci10pdense, 
  "CCD_5p-LVsparse" = SHD_ccd5pLVsparse, "FCI_5p-LVsparse" = SHD_fci5pLVsparse, 
  "CCI_5p-LVsparse" = SHD_cci5pLVsparse, "CCD_10p-LVsparse" = SHD_ccd10pLVsparse, 
  "FCI_10p-LVsparse" = SHD_fci10pLVsparse, "CCI_10p-LVsparse" = SHD_cci10pLVsparse, 
  "CCD_5p-LVdense" = SHD_ccd5pLVdense, "FCI_5p-LVdense" = SHD_fci5pLVdense, 
  "CCI_5p-LVdense" = SHD_cci5pLVdense, "CCD_10p-LVdense" = SHD_ccd10pLVdense, 
  "FCI_10p-LVdense" = SHD_fci10pLVdense, "CCI_10p-LVdense" = SHD_cci10pLVdense, 
  .id = "id"
) %>% 
  group_by(id) %>% 
  # get the average and sd
  summarise_all(list(means = mean, sds = sd)) %>%  
  # create columns
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2],
         netsize = stringr::str_split(condition, "-", simplify=T)[,1],
         latentvar = ifelse(stringr::str_detect(condition, "LV")==TRUE, 
                            "with LC", "without LC"),
         densities = stringr::str_remove(stringr::str_split
                                         (condition, "-", simplify=T)[,2], "LV")
  ) %>% 
  # convert it to a long format
  tidyr::pivot_longer(!c(algorithm, condition, id, netsize, latentvar, densities), 
                      names_to = "name", values_to = "value") %>% 
  # add sample size column (N) & clean up the column name 
  mutate(N = stringr::str_extract(
    stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  
  # bring the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric)) 



## =============================================================================
## 3. Create figures
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
SHDs %>%
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from = value) %>% 
  # create a ggplot object
  ggplot(aes(x = as.numeric(N), y = means, group = algorithm, 
             colour = algorithm, fill = algorithm)) +
  # add line plots
  geom_line(aes(group = algorithm)) +
  # add scattered points
  geom_point(size = 1) + 
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin = means+qnorm(0.25)*sds, ymax = means+qnorm(0.75)*sds), 
              alpha=0.15, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="N", y="", title = "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create facets
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels = c("p = 5", "p = 10")) ~ factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")),  
                      scales = "free_y", switch="y") +
  # manually specify the x-axis break
  scale_x_continuous(breaks=c(seq(50, 10000, by = 1000),10000)) +
  ## other x-axis labeling options that could be considered
  # 1) log10 transformation with unequally spaced labels
  # scale_x_continuous(trans= "log10", breaks=scales::breaks_log(n=15)) +
  # 2) log 10 transformation, trying to get the equal spaced labels with log transformation
  # scale_x_continuous(trans= "log10", breaks=c(50, 100, 250, 500, 1000, 2500, 5000, 10000)) +
  # 3) real N values with equally spaced labels
  # scale_x_continuous(breaks=c(seq(1000, 10000, by = 1000),10000)) +
  # 4) N as a factor and equally spaced labels: 
  # then need to change N as factor: ggplot(aes(x= factor(N, levels= c(50, 150, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000)) ...
  ggtitle("SHD") 

## save the plot
# ggsave(filename = "figures/SHD_alpha0.01_label10.pdf", width = 25, height = 13, dpi = 300, units = "cm")



## Precision figure
precision_plot <- pre_rec %>% 
  # extract the average precision
  filter(grepl("average_precision", metric)) %>% 
  # convert it to a wide format
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  # create a ggplot object
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
  # create facets
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")),  switch="y") +
  # manually specify the x-axis break
  scale_x_continuous(breaks=c(seq(50, 10000, by = 1000),10000)) +
  labs(title = "(a) Precision", x = "", y = "")



## Recall figure
recall_plot <- pre_rec %>% 
  filter(grepl("average_recall", metric)) %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
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
  # create facets
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")),  switch="y") +
  # manually specify the x-axis break
  scale_x_continuous(breaks=c(seq(50, 10000, by = 1000),10000)) +
  labs(title = "(b) Recall", x = "", y = "")



## Uncertainty figure
uncertainty_plot <- uncertainties %>%
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= as.numeric(N), y=means, group = algorithm, colour = algorithm, fill = algorithm)) +
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
  labs(x="N", y="", title = "") +
  # apply the theme
  theme_minimal() +
  MyTheme + 
  # create facets
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")),  
                      scales = "free_y", switch="y") +
  # manually specify the x-axis break
  scale_x_continuous(breaks=c(seq(50, 10000, by = 1000),10000)) +
  ggtitle("(c) Uncertainty")

# combine the plots
ggpubr::ggarrange(precision_plot, recall_plot, uncertainty_plot, nrow=3, 
                  common.legend = TRUE, legend = "bottom")

## save the plot
# ggsave(filename = "figures/prec-recall-uncer_alpha0.01_10label.pdf", width = 25, height = 30, dpi = 300, units = "cm")


