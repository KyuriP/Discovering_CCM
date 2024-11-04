## =============================================================================
## Description
#
# This script concerns applying the CCD, FCI and CCI algorithms to an 
# empirical data from McNally et al. (2017) in order to test the practical 
# applicability of the algorithms.
#
# It contains the code to create "Figure 18. Estimated statistical 
# network model and PAGs"
#
# The data can be obtained from : https://www.cambridge.org/core/journals/psychological-medicine/article/abs/comorbid-obsessivecompulsive-disorder-and-depression-a-bayesian-network-approach/DAA4E2352A9E26809A4EAE35C366E900#supplementary-materials
## =============================================================================
# The content is as follows.
# 1. Preparation: we source the necessary functions and packages. Also, we load
#    the data from the `data` folder and store each of the depression and OCD 
#    symptoms separately.
# 
# 2. Estimate GGM with GLASSO: we estimate the statistical network model (GGM) 
#    using glasso regularization with the `qgraph` package.
#
# 3. Estimate PAGs using CCD, FCI and CCI: we apply the causal discovery algorithms
#    and correspondingly estimate the PAGs. Here, we use alpha level of 0.01.
# 
# 4. Extra: we experiment with paranormal transformation to examine if it helps 
#    in approximating the distribution to Gaussian. Please note that this analysis
#    is an additional exploration and not part of the paper.
## =============================================================================


## =============================================================================
## 1. Preparation
## =============================================================================
## load necessary packages
library(qgraph)
library(pcalg)
library(ggplot2)
library(dplyr)
# remotes::install_github("KyuriP/CCI_KP")
library(CCI.KP)
library(furrr)

## source all the necessary functions
source("utils/CCD_fnc.R")
source("utils/plot_fnc.R")

## import the empirical data
mcnally <- read.csv("empirical_example/data/McNally.csv")

# separate depression / OCD symptoms
# (original data contains both depression and OCD symptoms)
# (in the paper, we only use depression symptoms)
depression <- mcnally[,1:16] %>% apply(., 2, as.numeric) # convert it to numeric (for CCD)
ocd <- mcnally[,17:26] %>% apply(., 2, as.numeric)


## =============================================================================
## 2. Estimate GGM with GLASSO 
## =============================================================================
## estimate GGM via graphical LASSO on depression symptoms
cordep <- cor(depression)
# found the optimal sparsity with gamma = 1
glassoFitdep <- EBICglasso(cordep, n = nrow(depression), gamma = 1)
qgraph(glassoFitdep, layout = "spring", theme="colorblind",
       nodeNames = colnames(depression), legend.cex = 0.4)


## =============================================================================
## 3. Estimate PAGs using CCD, FCI and CCI
## =============================================================================
# set the seed
set.seed(123)
# set alpha to 0.01
alpha <- 0.01

## estimate the PAG on depression symptoms by running CCD
# run CCD
ccd_mcnally_dep <- ccdKP(df=depression, dataType = "continuous", alpha = alpha)
# create an adjacency matrix for PAG
mat_mcnally_dep <- CreateAdjMat(ccd_mcnally_dep, p = ncol(depression))
# plot the PAG
pag_mcnally_dep <- plotPAG(ccd_mcnally_dep, mat_mcnally_dep)

## estimate the PAG on depression symptoms by running FCI
fci(list(C = cor(depression), n = nrow(depression)), gaussCItest, alpha = alpha, 
    labels = colnames(depression), verbose=TRUE) %>% .@amat %>% plotAG 

## estimate the PAG on depression symptoms by running CCI
cci(list(C = cor(depression), n = nrow(depression)), gaussCItest, alpha = alpha, 
    labels = colnames(depression), p = ncol(depression), verbose=TRUE) %>% .$maag %>% plotAG


## =============================================================================
## 4. Extra: Transformed data distribution
## =============================================================================
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

# distributions of depression symptoms (figureJ1 in the paper)
dist <- as.data.frame(depression) %>%
  # convert it to a long format
  tidyr::pivot_longer(where(is.numeric)) %>%
  # create ggplot object
  ggplot(aes(x = value)) +
  # add histograms
  geom_histogram(aes(y = after_stat(density)),bins = 10) +
  # add density
  geom_density(adjust=1.5)+
  # create a facet
  facet_wrap(~name) +
  # apply themes
  theme_minimal() +
  MyTheme2
# ggsave(dist, filename = "results/dep_dist.pdf", width = 20, height = 13, dpi = 300, units = "cm")

# original data distribution
p1 <- as.data.frame(depression) %>%
  # convert it to a long format
  tidyr::pivot_longer(where(is.numeric)) %>%
  # create ggplot object
  ggplot(aes(x = value)) +
  # add histograms
  geom_histogram(bins = 10) +
  # create a facet
  facet_wrap(~name) +
  # apply themes
  theme_minimal() +
  ggtitle("(a) Original data") +
  MyTheme2

## paranormal transformation using huge package
# some are not normal --> makes the data semi-parametric Gaussian
transformed_dat <- huge::huge.npn(depression) %>% as.data.frame()

## check the distributions of transformed data
p2 <- transformed_dat %>%
  # convert it to a long format
  tidyr::pivot_longer(where(is.numeric)) %>%
  # create ggplot object
  ggplot(aes(x = value)) +
  # add histograms
  geom_histogram(bins = 10) +
  # create a facet
  facet_wrap(~name) +
  # apply themes
  theme_minimal() +
  ggtitle("(b) Transformed data") +
  MyTheme2

# combine figures
ggpubr::ggarrange(p1, p2)
#ggsave(filename = "results/transformdat.pdf", width = 25, height = 13, dpi = 300, units = "cm")