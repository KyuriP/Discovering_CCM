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
depression <- mcnally[,1:16] %>% apply(., 2, as.numeric) # convert it to numeric (for CCD)
ocd <- mcnally[,17:26] %>% apply(., 2, as.numeric)

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
ccd_mcnally_dep <- ccdKP(df=depression, dataType = "continuous", alpha=0.01)
# create an adjacency matrix for PAG
mat_mcnally_dep <- CreateAdjMat(ccd_mcnally_dep, p = ncol(depression))
# plot the PAG
pag_mcnally_dep <- plotPAG(ccd_mcnally_dep, mat_mcnally_dep)

## estimate the PAG on depression symptoms by running FCI
fci(list(C = cor(depression), n = nrow(depression)), gaussCItest, alpha=0.01, 
    labels = colnames(depression), verbose=TRUE) %>% .@amat %>% plotAG 

## estimate the PAG on depression symptoms by running CCI
cci(list(C = cor(depression), n = nrow(depression)), gaussCItest, alpha=0.01, 
    labels = colnames(depression), p = ncol(depression), verbose=TRUE) %>% .$maag %>% plotAG



## =======================================
## 4. Extra: Transformed data distribution
## =======================================
# specify my custom plotting theme
MyTheme2 <-  theme(plot.title = element_text(family = "Palatino", size = 14, hjust=0.5),
                   plot.subtitle = element_text(face = "italic", family = "Palatino", size = 15, hjust=0.5),
                   axis.text=element_text(face = "bold",family = "Palatino", size = 11),
                   #axis.text.x = element_text(angle = 45, hjust = 1.2, vjust =1.2),
                   axis.title = element_text(face = "bold",family = "Palatino", size = 12),
                   legend.text = element_text(face = "bold", family = "Palatino", size = 12),
                   legend.position="bottom",
                   strip.text = element_text(size=12, family = "Palatino"),
                   strip.background = element_rect(fill="#f0f0f0", linetype = "solid", color="gray"),
                   strip.placement = "outside",
                   panel.border = element_rect(color = "#DCDCDC", fill = NA)
)

# original data distribution
p1 <- depression %>%
  tidyr::pivot_longer(where(is.numeric)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 10) +
  facet_wrap(~name) +
  theme_minimal() +
  ggtitle("(a) Original data") +
  MyTheme2

## transform data
# some are not normal --> makes the data semi-parametric Gaussian using huge package
transformed_dat <- huge::huge.npn(depression) %>% as.data.frame()

## check the distributions of transformed data
p2 <- transformed_dat %>%
  tidyr::pivot_longer(where(is.numeric)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 10) +
  facet_wrap(~name) +
  theme_minimal() +
  ggtitle("(b) Transformed data") +
  MyTheme2

ggpubr::ggarrange(p1,p2)
#ggsave(filename = "results/transformdat.pdf", width = 25, height = 13, dpi = 300, units = "cm")