## =============================================================================
## Purpose: to investigate the unexpected patterns in the B5 dense conditions 
# from the main simulation study. 
#
#
## Description
#
# This script contains code for investigating the observed phenomenon 
# where the performance in the "B5 dense" case worsens as the sample size N increases.
# We presumed that the inducing path between X2 and X5 plays a role in
# achieving the correct orientation when N is relatively small.
#
# To test our hypothesis, we examined the partial correlations between
# X2 and X5 per different sample sizes and looked into what was happening in
# each step of the algorithm.
#
# The first part concerns examining the partial correlations between X2 and x5
# and the second part concerns conditional independence test; 
# whether it is significant or not significant as sample size (N) becomes larger.
## =============================================================================


# source the simulation study results
source("simulation/01_main_simulation.R")

# load packages
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ppcor)


## ======================================
## PART I: examining partial correlations
## ======================================
## Check every partial correlations between X2 and X5
marginal <- c()
partialX1 <- c()
partialX3 <- c()
partialX4 <- c()
partialX1X4 <- c()
partialX1X3 <- c()
partialX3X4 <- c()
partialX1X3X4 <- c()
for(j in 1:500){
  data <- simdata_5pdense[[10]][[j]]
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
                  axis.text=element_text(face = "bold",family = "Palatino", size = 11),
                  axis.text.x = element_text(hjust = 1.2, vjust =1.2),
                  axis.title = element_text(face = "bold",family = "Palatino", size = 12),
                  strip.text = element_text(face="bold", size=13, family = "Palatino"),
                  strip.background = element_rect(fill="#f0f0f0", linetype = "solid", color="gray"),
                  strip.placement = "outside",
                  panel.border = element_rect(color = "#DCDCDC", fill = NA),
                  panel.spacing.y = unit(4, "mm"),panel.spacing.x = unit(4, "mm")
)

# ggplot histograms
bind_cols("marginal" = marginal, "partial_X1" = partialX1, "partial_X3" = partialX3, 
          "partial_X4" = partialX4, "partial_X1&X3" = partialX1X3, "partial_X1&X4" = partialX1X4, 
          "partial_X3&X4" = partialX3X4, "partial_X1&X3&X4" = partialX1X3X4) %>% 
  tidyr::pivot_longer(cols=everything(), names_to = "id", values_to = "cors") %>% 
  group_by(id) %>% 
  mutate(means=mean(cors), min = min(cors)) %>% 
  ggplot(aes(x=cors)) + 
  geom_histogram(fill="gray", col="white") + 
  # adjust the x-axis ticks
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  # add mean values in text
  geom_text(aes(x = min, y = 50, hjust= -0.1, 
                label = paste0("mean = ", round(means,2)), family = "Palatino"), 
            size=5, col='#8b0000') +
  facet_wrap(~ factor(id, levels=c("marginal", "partial_X1", "partial_X3", 
                                   "partial_X4", "partial_X1&X3", "partial_X1&X4", 
                                   "partial_X3&X4", "partial_X1&X3&X4")), 
             nrow=2, ncol=4, scales = "free_x") +
  # add lines for the mean values
  geom_vline(aes(xintercept = means, group = id), colour = '#8b0000', linetype=2, linewidth=0.4 ) + 
  # apply the theme
  theme_minimal() + MyTheme3 +
  labs(x ="", y="count")
ggsave(filename = "results/partialcorr.pdf", width = 28, height = 15, dpi = 300, units = "cm")


## ============================================
## PART II: testing conditional independencies
## ============================================
## Create a table counting how many times it was significant given alpha = 0.05
count <- matrix(0, 8, 10)
indep <- matrix(NA, 8, 10)
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
  count <- ifelse(indep > 0.05, count+1, count)
}
# specify dim names
rownames(count) <- c("marginal", "conditional on X1", "conditional on X3", 
                     "conditional on X4", "conditional on X1 & X4", "conditional on X1 & X3", 
                     "conditional on X3 & X4", "conditional on X1 & X3 & X4")
colnames(count) <- c("N=50" , "N=150", "N=500", "N=1000", "N=2000", "N=3000", 
                     "N=4000", "N=5000", "N=7500", "N=10000")

# show table
count
