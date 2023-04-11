## =============================================================================
## Description
# 
# This script contains all the code for the secondary analysis with 
# randomly sampled coefficients for the regression matrix B "only positive values".
#
# Purpose: to assess the dependence of the results from the main analysis on 
# the specified weights by randomly sampling coefficients for the regression matrix 
# B with "only positive values" 
# (we simulate a positive manifold to examine if there is any impact on the results).
#
# The rest set up is exactly the same as `sampling_Bs.R`.
## =============================================================================

## ======================
## 0. Preparation
## ======================
# source the simulation study results
source("code/simulation_code.R")
# load packages
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggh4x)

## set the seed
set.seed(123)

## ======================
## 1. Generate data 
## ======================

# specify the sample sizes
N <- c(50, 150, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000)
# specify replication number
n <- 500
# specify alpha level
alpha <- 0.01



## function used to generate data from randomly sampled B weights (only positive)
sampleranB2 <- function(B, LV=NULL, seed=123){
  set.seed(seed)
  # variable number
  p <- ncol(B)
  # data storage
  simdat <- list()
  # number of diff sample size
  for(i in 1:length(N)){
    simdat[[i]] <- list()
    # number of repetition
    for(j in 1:n){
      # sample random weights for B from uniform distribution: unif[0.1, 0.8]
      ranB <- matrix((0.7*runif(p^2)+0.1), p, p)
      # assign 0 where B = 0
      ind <- which(B == 0, arr.ind = TRUE)
      ranB[ind] <- 0
      if(is.null(LV)){
        simdat[[i]][[j]] <- gen_dat(ranB, N = N[i])
        colnames(simdat[[i]][[j]]) <- colnames(B)
      } else {
        simdat[[i]][[j]] <- gen_dat(ranB, N = N[i])[,-LV]
        colnames(simdat[[i]][[j]]) <- colnames(B)[-LV]
      }
    }
  }
  return(simdat)
}

## create all simulated data using random B 
simdata_woLV <- list(B5sparse = B5sparse, B5dense = B5dense, B10sparse =  B10sparse, B10dense = B10dense) %>% 
  map(~
        sampleranB2(.x)
  )

simdata_5pwLV <- list(B5_lvsparse = B5_lvsparse, B5_lvdense = B5_lvdense) %>% 
  map(~
        sampleranB2(.x, LV = 6)
  )

simdata_10pwLV <- list(B10_lvsparse = B10_lvsparse, B10_lvdense = B10_lvdense) %>% 
  map(~
        sampleranB2(.x, LV = c(11, 12))
  )

simdatalist <- append(simdata_woLV, append(simdata_5pwLV, simdata_10pwLV))


## ======================
## 2. Run algorithms
## ======================
CCDres_pos <- simdatalist %>% 
  map_depth(3, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)) %>% 
  map_depth(3, ~CreateAdjMat(.x, length(.x$nodes)))


FCIres_pos <- simdatalist %>% 
  map_depth(3, ~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
  )

CCIres_pos <- simdatalist %>% 
  map_depth(3, ~ cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag 
  )

# true adj.matrices
truemods <- list(trueag_5psparse, trueag_5pdense, trueag_10psparse, trueag_10pdense, trueag_5psparseLV, trueag_5pdenseLV, trueag_10psparseLV, trueag_10pdenseLV)

# load results
load("data/largedata_posB/CCDres3_pos.Rdata") # CCDres
load("data/largedata_posB/FCIres3_pos.Rdata") # FCIres
load("data/largedata_posB/CCIres3_pos.Rdata") # CCIres




## =========================
## 3. Evaluate performance
## =========================

## SHD
CCDshd <- list()
for(i in 1:length(CCDres)){
  CCDshd[[i]] <- CCDres[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N))  %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCDshd) <- names(CCDres)


FCIshd <- list()
for(i in 1:length(FCIres)){
  FCIshd[[i]] <- FCIres[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(FCIshd) <- names(FCIres)


CCIshd <- list()
for(i in 1:length(CCIres)){
  CCIshd[[i]] <- CCIres[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCIshd) <- names(CCIres)

## combine the SHDs
SHD_ranB <- bind_rows(CCD = CCDshd, FCI = FCIshd, CCI = CCIshd, .id="id") %>% 
  tidyr::pivot_longer(cols = -c(id), names_to = "condition", values_to = "value") %>% 
  mutate(
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LC", "without LC"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>%  
  tidyr::unnest(value) %>% 
  tidyr::pivot_longer(cols = starts_with("N ="), names_to = "n", values_to = "value") %>% 
  mutate(statistics = stringr::str_split(n, "_", simplify=T)[,2],
         n = as.numeric(stringr::str_extract_all(n, "[0-9]+")))  %>% 
  relocate(where(is.character), .before = where(is.numeric))



## Precision
CCDprec <- list()
for(i in 1:length(CCDres)){
  CCDprec[[i]] <- CCDres[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCDprec) <- names(CCDres)

FCIprec <- list()
for(i in 1:length(FCIres)){
  FCIprec[[i]] <- FCIres[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(FCIprec) <- names(FCIres)

CCIprec <- list()
for(i in 1:length(CCIres)){
  CCIprec[[i]] <- CCIres[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCIprec) <- names(CCIres)

## combine the precisions
prec_ranB <- bind_rows(CCD = CCDprec, FCI = FCIprec, CCI = CCIprec, .id="id") %>% 
  tidyr::pivot_longer(cols = -c(id), names_to = "condition", values_to = "value") %>% 
  mutate(
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LC", "without LC"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>%  
  tidyr::unnest(value) %>% 
  tidyr::pivot_longer(cols = starts_with("N ="), names_to = "n", values_to = "value") %>% 
  mutate(statistics = stringr::str_split(n, "_", simplify=T)[,2],
         n = as.numeric(stringr::str_extract_all(n, "[0-9]+")))  %>% 
  relocate(where(is.character), .before = where(is.numeric))


## Recall
CCDrec <- list()
for(i in 1:length(CCDres)){
  CCDrec[[i]] <- CCDres[[i]] %>% 
    map_depth(2, ~recall2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCDrec) <- names(CCDres)

FCIrec <- list()
for(i in 1:length(FCIres)){
  FCIrec[[i]] <- FCIres[[i]] %>% 
    map_depth(2, ~ recall2(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(FCIrec) <- names(FCIres)

CCIrec <- list()
for(i in 1:length(CCIres)){
  CCIrec[[i]] <- CCIres[[i]] %>% 
    map_depth(2, ~ recall2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCIrec) <- names(CCIres)

## combine the recalls
rec_ranB <- bind_rows(CCD = CCDrec, FCI = FCIrec, CCI = CCIrec, .id="id") %>% 
  tidyr::pivot_longer(cols = -c(id), names_to = "condition", values_to = "value") %>% 
  mutate(
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LC", "without LC"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>%  
  tidyr::unnest(value) %>% 
  tidyr::pivot_longer(cols = starts_with("N ="), names_to = "n", values_to = "value") %>% 
  mutate(statistics = stringr::str_split(n, "_", simplify=T)[,2],
         n = as.numeric(stringr::str_extract_all(n, "[0-9]+")))  %>% 
  relocate(where(is.character), .before = where(is.numeric))


## Uncertainty
CCDunc <- list()
for(i in 1:length(CCDres)){
  CCDunc[[i]] <- CCDres[[i]] %>% 
    map_depth(2, ~uncertainty(.x)) %>% 
    do.call("cbind", .) %>% 
    apply(., 2, unlist) %>% 
    as_tibble %>% 
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means = mean, sds = sd))  %>% t()
}
names(CCDunc) <- names(CCDres)


FCIunc <- list()
for(i in 1:length(FCIres)){
  FCIunc[[i]] <- FCIres[[i]] %>% 
    map_depth(2, ~uncertainty(.x))  %>% 
    do.call("cbind", .) %>% 
    apply(., 2, unlist) %>% 
    as_tibble %>% 
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means = mean, sds = sd)) %>% t()
}
names(FCIunc) <- names(FCIres)


CCIunc <- list()
for(i in 1:length(CCIres)){
  CCIunc[[i]] <- CCIres[[i]] %>% 
    map_depth(2, ~uncertainty(.x)) %>% 
    do.call("cbind", .) %>% 
    apply(., 2, unlist) %>% 
    as_tibble %>% 
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means = mean, sds = sd)) %>% t() 
}
names(CCIunc) <- names(CCIres)


## combine the uncertainties
unc_ranB <- bind_rows(CCD = CCDunc, FCI = FCIunc, CCI = CCIunc, .id="id") %>% 
  # n = repeat N: 3 algo * 2 (mean & sd), repeat c(means,sd) by length of N * 3 algo
  mutate(n = rep(N,6), statistics = rep(c("means", "sds"), each = length(N), times = 3)) %>% 
  tidyr::pivot_longer(cols = -c(id, n, statistics), names_to = "condition", values_to = "value") %>% 
  mutate(
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LC", "without LC"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>% 
  relocate(where(is.character), .before = where(is.numeric))



## ===================================
## 4. Plot the figures for comparison
## ===================================
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


## plot SHD
shdplot_ranB_pos <- SHD_ranB %>%
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= as.numeric(n), y=means, group = id, colour = id, fill = id)) +
  geom_line(aes(group = id)) +
  geom_point(size=1) + 
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="", y="", title = "") +
  theme_minimal() +
  MyTheme + 
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), labels=c("p = 5", "p = 10")) ~ factor(latentvar, levels = c("without LC", "with LC")) + factor(densities, levels=c("sparse", "dense")), scales = "free_y", switch="y") +
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  ggtitle("(a) SHD")  +
  guides(color = "none", fill = "none")
# save the plot
# ggsave(shdplot_ranB_pos, filename = "results/samplingbeta_pos_shd.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")



## plot precision
precisionplot_ranB_pos <- prec_ranB %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= as.numeric(n), y=means, group = id, colour = id, fill=id)) +
  geom_line(aes(group = id)) +
  geom_point(size=1) +
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  theme_minimal() +
  MyTheme + 
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), labels=c("p = 5", "p = 10")) ~ factor(latentvar, levels = c("without LC", "with LC")) + factor(densities, levels=c("sparse", "dense")), switch="y") +
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  labs(title = "(b) Precision", x = "", y = "") +
  guides(color = "none", fill = "none")
# save the plot
# ggsave(precisionplot_ranB_pos, filename = "results/samplingbeta_pos_prec.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")


## plot recall
recallplot_ranB_pos <- rec_ranB %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= as.numeric(n), y=means, group = id, colour = id, fill=id)) +
  geom_line(aes(group = id)) +
  geom_point(size=1) +
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  theme_minimal() +
  MyTheme + 
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), labels=c("p = 5", "p = 10")) ~ factor(latentvar, levels = c("without LC", "with LC")) + factor(densities, levels=c("sparse", "dense")), switch="y") +
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  labs(title = "(c) Recall", x = "", y = "")+
  guides(color = "none", fill = "none")
# save the plot
# ggsave(recallplot_ranB_pos, filename = "results/samplingbeta_pos_recall.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")

## plot uncertainty 
uncertaintyplot_ranB_pos <- unc_ranB %>%  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= as.numeric(n), y=means, group = id, colour = id, fill = id)) +
  geom_line(aes(group = id)) +
  geom_point(size=1) + 
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x="N", y="", title = "") +
  theme_minimal() +
  MyTheme + 
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), labels=c("p = 5", "p = 10")) ~ factor(latentvar, levels = c("without LC", "with LC")) + factor(densities, levels=c("sparse", "dense")),  scales = "free_y", switch="y") +
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  ggtitle("(d) Uncertainty") 
# save the plot
# ggsave(uncertaintyplot_ranB_pos, filename = "results/samplingbeta_pos_unc.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")

# combine the plots
ggpubr::ggarrange(shdplot_ranB_pos, precisionplot_ranB_pos, 
                  recallplot_ranB_pos, uncertaintyplot_ranB_pos, 
                  nrow=4, common.legend = TRUE, legend = "bottom")
#ggsave(filename = "results/samplingbeta_pos_result.pdf", width = 25, height = 35, dpi = 300, units = "cm")

ggpubr::ggarrange(shdplot_ranB_pos, precisionplot_ranB_pos, 
                  nrow=2, common.legend = TRUE, legend = "bottom")
#ggsave(filename = "results/samplingbeta_pos_result1.pdf", width = 25, height = 35, dpi = 300, units = "cm")
ggpubr::ggarrange(recallplot_ranB_pos, uncertaintyplot_ranB_pos, 
                  nrow=2, common.legend = TRUE, legend = "bottom")
#ggsave(filename = "results/samplingbeta_pos_result2.pdf", width = 25, height = 35, dpi = 300, units = "cm")
