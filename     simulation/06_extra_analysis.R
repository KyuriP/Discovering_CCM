## =============================================================================
## Purpose
#
# To investigate the unexpected patterns (i.e., brief dips and spikes in the 
# performance graphs) in the "B5 dense" conditions from the main simulation study. 
#
#
## Description
#
# This script contains code for investigating the observed phenomenon where the 
# performance in the "B5 dense" cases deteriorates as the sample size N increases.
#
# We presumed that the inducing path between X2 and X5 plays a role in
# achieving the correct orientation when N is relatively small.
#
# To test our hypothesis, we examined the partial correlations between X2 and X5 
# for different sample sizes and looked into the results of conditional
# independence tests.
## =============================================================================
# The content is as follows.
# 0. Preparation: Source and load necessary functions & packages.
# 
# 1. Examining partial correlation: Examine all partial correlations 
#    between X2 and X5. And create histograms.
#
# 2. Testing conditional independencies: Check the results of conditional 
#    independence tests; whether it is significant or not significant as 
#    sample size (N) becomes larger.
## =============================================================================


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
library(ppcor)

# set the seed
set.seed(123)


## =============================================================================
## 1. Examining partial correlations
## =============================================================================
## Check every partial correlations between X2 and X5
# storage
marginal <- c()
partialX1 <- c()
partialX3 <- c()
partialX4 <- c()
partialX1X4 <- c()
partialX1X3 <- c()
partialX3X4 <- c()
partialX1X3X4 <- c()
# extract marginal/partial correlations
for(j in 1:500){
  data <- simdata_5pdense[[10]][[j]] # using simulated data with N = 10000
  marginal[j] <- cor.test(data[,"X2"], data[,"X5"])$estimate
  partialX1[j] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X1")])$estimate
  partialX3[j] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X3")])$estimate
  partialX4[j] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X4")])$estimate
  partialX1X4[j] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X4")])$estimate
  partialX1X3[j] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X3")])$estimate
  partialX3X4[j] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X3", "X4")])$estimate
  partialX1X3X4[j] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X3", "X4")])$estimate
}


## Create histograms of partial correlations
# specify the common figure theme
MyTheme3 <-  theme(
                  axis.text=element_text(face = "bold",family = "Palatino", size = 15),
                  axis.text.x = element_text(hjust = 1.2, vjust =1.2),
                  axis.title = element_text(face = "bold",family = "Palatino", size = 16),
                  strip.text = element_text(face="bold", size=17, family = "Palatino"),
                  strip.background = element_rect(fill="#f0f0f0", linetype = "solid", color="gray"),
                  strip.placement = "outside",
                  panel.border = element_rect(color = "#DCDCDC", fill = NA),
                  panel.spacing.y = unit(4, "mm"),panel.spacing.x = unit(4, "mm")
)

## ggplot histograms
bind_cols("marginal" = marginal, "partial_X1" = partialX1, "partial_X3" = partialX3, 
          "partial_X4" = partialX4, "partial_X1 & X3" = partialX1X3, "partial_X1 & X4" = partialX1X4, 
          "partial_X3 & X4" = partialX3X4, "partial_X1 & X3 & X4" = partialX1X3X4) %>% 
  # convert it to a long format
  tidyr::pivot_longer(cols=everything(), names_to = "id", values_to = "cors") %>% 
  group_by(id) %>% 
  # store the average and minimum values
  mutate(means= round(mean(cors),3), min = min(cors)) %>% 
  # create ggplot object
  ggplot(aes(x=cors)) + 
  # add histogram
  geom_histogram(fill="gray", col="white", bins=25) + 
  # add mean values in text
  geom_text(aes(x = min, y = 60, hjust= -0.01, 
                # label formatting
                label = paste0("M = ", formatC(means, 3, format="f")), 
                family = "Palatino", fontface = "italic"), size=5.5, col='#8b0000',
            # prevent overplotting text
            check_overlap = T) +
  # create faceted plot
  facet_wrap(~ factor(id, levels=c("marginal", "partial_X1", "partial_X3", 
                                   "partial_X4", "partial_X1 & X3", "partial_X1 & X4", 
                                   "partial_X3 & X4", "partial_X1 & X3 & X4")), 
             nrow=2, ncol=4, scales = "free_x") +
  # adjust the x-axis ticks
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  # add lines for the mean values
  geom_vline(aes(xintercept = means, group = id), colour = '#8b0000', linetype=2, linewidth=0.4 ) + 
  # apply the theme
  theme_minimal() + MyTheme3 +
  labs(x ="", y="count")
# ggsave(filename = "figures/partialcorr.pdf", width = 35, height = 18, dpi = 300, units = "cm")


## =============================================================================
## 2. Testing conditional independencies
## =============================================================================
## Create a table counting how many times it was significant given alpha = 0.05
# storage
count <- matrix(0, 8, 10)
indep <- matrix(NA, 8, 10)
# extract p-values
for(j in 1:500){
  for(i in 1:10){
    data <- simdata_5pdense[[i]][[j]]
    indep[1,i] <- cor.test(data[,"X2"], data[,"X5"])$p.value
    indep[2,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X1")])$p.value
    indep[3,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X3")])$p.value
    indep[4,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X4")])$p.value
    indep[5,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X4")])$p.value
    indep[6,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X3")])$p.value
    indep[7,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X3", "X4")])$p.value
    indep[8,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X3", "X4")])$p.value
  }
  # count how many p-value > 0.05
  count <- ifelse(indep > 0.05, count+1, count)
}
# specify dim names
rownames(count) <- c("marginal", "conditional on X1", "conditional on X3", 
                     "conditional on X4", "conditional on X1 & X4", "conditional on X1 & X3", 
                     "conditional on X3 & X4", "conditional on X1 & X3 & X4")
colnames(count) <- c("N=50" , "N=150", "N=500", "N=1000", "N=2000", "N=3000", 
                     "N=4000", "N=5000", "N=7500", "N=10000")

# show the table of proportion
count/500 
