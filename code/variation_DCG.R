
## =============================================================================
## Description
#
# This script contains the code to create the Figure 2 and Figure 3 in the manuscript "Follow-up Study".
# Figure 2: Distribution of densities per each set of DCGs.
# Figure 3: Average degree for every node per each set of DCGs with 95% confidence interval.
## =============================================================================


## =============================================================================
# Preparation
## =============================================================================
## load necessary packages
library(qgraph)
library(pcalg)
library(dplyr)
library(ggplot2)

## source all the necessary functions
source("code/R/variation_fnc.R")
source("code/simulation.R")

## load all the equivalence set of DCGs from the original simulation study
load("data/equiv4p.RData")
load("data/equiv4p_high.RData")
load("data/equiv5p.RData")
load("data/equiv5p_high.RData")
load("data/equiv6p.RData")
load("data/equiv6p_high.RData")


## =============================================================================
## Examine overall density variation (Figure 2)
#
# 1. Compute the density for DCGs.
# 2. Compute the density for true models.
# 3. Plot the distribution of densities of DCGs
#    and add a line for true density as well as for average DCG density.
## =============================================================================

## 1) compute the density for DCGs per model
# 4 node-sparse
denvar4p <-DCGdensities(equiv4p)
# 4 node-dense
denvar4phigh <-DCGdensities(equiv4p_high)
# 5 node-sparse
denvar5p <-DCGdensities(equiv5p)
# 5 node-dense
denvar5phigh <-DCGdensities(equiv5p_high)
# 6 node-sparse
denvar6p <-DCGdensities(equiv6p)
# 6 node-dense
denvar6phigh <-DCGdensities(equiv6p_high)
# put them all together in a list
modeldensities <- list(denvar4p, denvar4phigh, denvar5p, denvar5phigh, denvar6p, denvar6phigh)
names(modeldensities) <- c("(a) 4 nodes - sparse", "(b) 4 nodes - dense", "(c) 5 nodes - sparse", "(d) 5 nodes - dense", "(e) 6 nodes - sparse", "(f) 6 nodes - dense")

## 2) compute the true model densities
truemodels <- list(B4, B4_high, B5, B5_high, B6, B6_high) %>%
  purrr::map(~truemoddensity(.)) %>% unlist()

## 3) plot the distributions of density per model
# specify bin numbers
bins <- c(15, 15, 15, 15, 20, 20)
# storage for plots
plotlist1 <- list()
for(i in 1:length(modeldensities)){
  # get the size of equivalence class
  N <- modeldensities[[i]]$class_size[1]
  plotlist1[[i]] <- as.data.frame(modeldensities[[i]]) %>%
    mutate(trueden = truemodels[i]) %>%
    # plot density distribution using histogram
    ggplot(aes(x=densities)) +
    geom_histogram(color="gray", fill="lightblue", bins = bins[i]) +
    xlim(0.25, 1.01) +
    # add a line indicating average density of DCGs
    geom_vline(aes(xintercept = avg_density, color = "Average Density of DCGs"),
               lwd = 0.3, linetype=2) +
    # add a line indicating true density
    geom_vline(aes(xintercept = trueden, color = "Density of True Model"),
               lwd = 0.3, linetype=2) +
    # add extra aesthetics for the plot
    theme_classic() + labs(title = names(modeldensities)[i],
                           subtitle = paste("Size of equivalence class = ", N),
                           x = "density") +
    scale_color_manual(name = "", values = c("Density of True Model" = "red",
                                             "Average Density of DCGs" = "darkblue"))+
    theme(plot.title = element_text(size = 12, face="bold"),
          plot.subtitle=element_text(size=12, face="italic"),
          legend.text = element_text(size = 11))

}
# put all plots together (results in Figure2)
ggpubr::ggarrange(plotlist=plotlist1, nrow=3, ncol=2, common.legend = T, legend = "bottom")
# save the plot
#ggsave("manuscript/img/densityvariation.png", dpi=600)


## =============================================================================
## Examine degree centrality variation (Figure 3)
#
# 1. Compute the degree for all nodes in DCGs.
# 2. Compute the standard errors (SE) per model.
# 3. Plot the average degree centrality with the 95% confidence interval.
## =============================================================================

## 1) obtain degrees per node
## 2) and compute the SE  for all DCGs
# 4 node-sparse
degvar4p <- DCGdegrees(equiv4p) %>%
  # compute SE
  mutate(se = deg_sd * qnorm(0.975) / length(equiv4p), name = "(a) 4 nodes - sparse")
# 4 node-dense
degvar4phigh <- DCGdegrees(equiv4p_high) %>%
  mutate(se = deg_sd * qnorm(0.975) / length(equiv4p_high), name = "(b) 4 nodes - dense")
# 5 node-sparse
degvar5p <- DCGdegrees(equiv5p) %>%
  mutate(se = deg_sd * qnorm(0.975) / length(equiv5p), name = "(c) 5 nodes - sparse")
# 5 node-dense
degvar5phigh <- DCGdegrees(equiv5p_high) %>%
  mutate(se = deg_sd * qnorm(0.975) / length(equiv5p_high), name = "(d) 5 nodes - dense")
# 6 node-sparse
degvar6p <- DCGdegrees(equiv6p) %>%
  mutate(se = deg_sd * qnorm(0.975) / length(equiv6p), name = "(e) 6 nodes - sparse")
# 6 node-dense
degvar6phigh <- DCGdegrees(equiv6p_high) %>%
  mutate(se = deg_sd * qnorm(0.975) / length(equiv6p_high), name = "(f) 6 nodes - dense")
# put all degrees per model together in a list
modeldegrees <- list(degvar4p, degvar4phigh, degvar5p, degvar5phigh, degvar6p, degvar6phigh)

## 3) create the degree centrality plots for each model with 95% confidence interval
# storage for plots
plotlist2 <- modeldegrees %>%
  purrr::map(~
               ggplot(data =., aes(x=node, y=avg_degree)) +
               # we multiply SEs by 100 to scale up (hard to visualize them otherwise)
               # add the errorbars (95% CIs)
               geom_errorbar(aes(ymin=avg_degree-se*100, ymax=avg_degree+se*100), width=0.1,
                             color = "darkblue") +
               geom_line(group=1, color = "darkblue") + geom_point(size=0.7) +
               # add extra aesthetics for the plot
               theme_minimal() + ylim(0, 5) +
               labs(title = .$name,  y = "average degree") +
               theme(plot.title = element_text(size = 10, face="bold"))
  )
# put all plots together (results in Figure3)
ggpubr::ggarrange(plotlist=plotlist2, nrow=3, ncol=2, common.legend = T, legend = "bottom")
# save the plot
#ggsave("manuscript/img/degreevariation.png", dpi=600)
