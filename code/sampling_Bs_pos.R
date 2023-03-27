## =============================================================================
## Description
# 
# This script contains all the code for the secondary analysis with 
# randomly sampled coefficients for the regression matrix B.
#
# Purpose: to gauge how much the results of the main analysis are dependent on 
# the specific weights we specified.
#
# As is the case with the main analysis, there are in total 8 models and
# we generate 500 datasets from each model.
#
# The content is as follows:
# 0. Preparation: we source and load necessary functions & packages and generate data.
# 1. Generate data: we generate data based on randomly sampled B matrix at every iteration. 
# 2. Run algorithms: we again run three algorithms CCD, FCI, and CCI then estimate PAGs.
# 3. Evaluate performance: we compute structural Hamming distance, precision, recall,
# and uncertainty rate for each condition.
# 4. Create figures: we create figures for each evaluation metric comparing the 
# performance of each algorithm per condition.
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
N <- c(50, 150, 500, 1000, 1500, 2000, 3000, 4000, 5000, 10000)
# specify replication number
n <- 1e3
# specify alpha level
alpha <- 0.05


sampleranB <- function(B, seed=123){
  # variable number
  p <- ncol(B)
  # sample random weight for B from uniform distribution (unif[-0.8, -0.1] U unif[0.1, 0.8])
  ranB <- matrix((0.8*runif(p^2)+0.1)*sample(c(-1,1), p^2,replace=TRUE), p, p)
  # assign 0 where B = 0
  ind <- which(B == 0, arr.ind = TRUE)
  ranB[ind] <- 0
  # generate data n times for each N
  simdat <- N %>% future_map(function(z) {
    replicate(n = n,
              expr = gen_dat(ranB, N = z),  
              simplify = FALSE)
  }, .options = furrr_options(seed=seed))
  return(simdata = simdat)
}



### if randomB is used for each iteration..
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
CCDres <- simdatalist %>% 
  map_depth(3, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)) %>% 
  map_depth(3, ~CreateAdjMat(.x, length(.x$nodes)))


FCIres <- simdatalist %>% 
  map_depth(3, ~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat 
  )

CCIres <- simdatalist %>% 
  map_depth(3, ~ cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag 
  )


# true adj.matrices
truemods <- list(trueag_5psparse, trueag_5pdense, trueag_10psparse, trueag_10pdense, trueag_5psparseLV, trueag_5pdenseLV, trueag_10psparseLV, trueag_10pdenseLV)


# load results
load("data/largedata_n500/CCDres2_randomB.Rdata")
load("data/largedata_n500/FCIres2_randomB.Rdata")
load("data/largedata_n500/CCIres2_randomB.Rdata")

load("data/largedata_randomB/CCDres_randomB.Rdata")
load("data/largedata_randomB/FCIres_randomB.Rdata")
load("data/largedata_randomB/CCIres_randomB.Rdata")



CCDshd1000 <- list()
for(i in 1:length(CCDres)){
  
}

## =========================
## 3. Evaluate performance
## =========================

## SHD
CCDshd <- list()
for(i in 1:length(CCDres2)){
  CCDshd[[i]] <- CCDres2[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N))  %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCDshd) <- names(CCDres2)


FCIshd <- list()
for(i in 1:length(FCIres2)){
  FCIshd[[i]] <- FCIres2[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(FCIshd) <- names(FCIres2)


CCIshd <- list()
for(i in 1:length(CCIres2)){
  CCIshd[[i]] <- CCIres2[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCIshd) <- names(CCIres2)

## combine the SHDs
SHD_ranB <- bind_rows(CCD = CCDshd, FCI = FCIshd, CCI = CCIshd, .id="id") %>% 
  tidyr::pivot_longer(cols = -c(id), names_to = "condition", values_to = "value") %>% 
  mutate(
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LV", "without LV"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>%  
  tidyr::unnest(value) %>% 
  tidyr::pivot_longer(cols = starts_with("N ="), names_to = "n", values_to = "value") %>% 
  mutate(statistics = stringr::str_split(n, "_", simplify=T)[,2],
         n = as.numeric(stringr::str_extract_all(n, "[0-9]+")))  %>% 
  relocate(where(is.character), .before = where(is.numeric))



## Precision
CCDprec <- list()
for(i in 1:length(CCDres2)){
  CCDprec[[i]] <- CCDres2[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCDprec) <- names(CCDres2)

FCIprec <- list()
for(i in 1:length(FCIres2)){
  FCIprec[[i]] <- FCIres2[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(FCIprec) <- names(FCIres2)

CCIprec <- list()
for(i in 1:length(CCIres2)){
  CCIprec[[i]] <- CCIres2[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCIprec) <- names(CCIres2)

## combine the precisions
prec_ranB <- bind_rows(CCD = CCDprec, FCI = FCIprec, CCI = CCIprec, .id="id") %>% 
  tidyr::pivot_longer(cols = -c(id), names_to = "condition", values_to = "value") %>% 
  mutate(
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LV", "without LV"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>%  
  tidyr::unnest(value) %>% 
  tidyr::pivot_longer(cols = starts_with("N ="), names_to = "n", values_to = "value") %>% 
  mutate(statistics = stringr::str_split(n, "_", simplify=T)[,2],
         n = as.numeric(stringr::str_extract_all(n, "[0-9]+")))  %>% 
  relocate(where(is.character), .before = where(is.numeric))


## Recall
CCDrec <- list()
for(i in 1:length(CCDres2)){
  CCDrec[[i]] <- CCDres2[[i]] %>% 
    map_depth(2, ~recall2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCDrec) <- names(CCDres2)

FCIrec <- list()
for(i in 1:length(FCIres2)){
  FCIrec[[i]] <- FCIres2[[i]] %>% 
    map_depth(2, ~ recall2(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(FCIrec) <- names(FCIres2)

CCIrec <- list()
for(i in 1:length(CCIres2)){
  CCIrec[[i]] <- CCIres2[[i]] %>% 
    map_depth(2, ~ recall2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCIrec) <- names(CCIres2)

## combine the recalls
rec_ranB <- bind_rows(CCD = CCDrec, FCI = FCIrec, CCI = CCIrec, .id="id") %>% 
  tidyr::pivot_longer(cols = -c(id), names_to = "condition", values_to = "value") %>% 
  mutate(
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LV", "without LV"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>%  
  tidyr::unnest(value) %>% 
  tidyr::pivot_longer(cols = starts_with("N ="), names_to = "n", values_to = "value") %>% 
  mutate(statistics = stringr::str_split(n, "_", simplify=T)[,2],
         n = as.numeric(stringr::str_extract_all(n, "[0-9]+")))  %>% 
  relocate(where(is.character), .before = where(is.numeric))


## Uncertainty
CCDunc <- list()
for(i in 1:length(CCDres2)){
  CCDunc[[i]] <- CCDres2[[i]] %>% 
    map_depth(2, ~uncertainty(.x)) %>% 
    do.call("cbind", .) %>% 
    apply(., 2, unlist) %>% 
    as_tibble %>% 
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means = mean, sds = sd))  %>% t()
}
names(CCDunc) <- names(CCDres2)


FCIunc <- list()
for(i in 1:length(FCIres2)){
  FCIunc[[i]] <- FCIres2[[i]] %>% 
    map_depth(2, ~uncertainty(.x))  %>% 
    do.call("cbind", .) %>% 
    apply(., 2, unlist) %>% 
    as_tibble %>% 
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means = mean, sds = sd)) %>% t()
}
names(FCIunc) <- names(FCIres2)


CCIunc <- list()
for(i in 1:length(CCIres2)){
  CCIunc[[i]] <- CCIres2[[i]] %>% 
    map_depth(2, ~uncertainty(.x)) %>% 
    do.call("cbind", .) %>% 
    apply(., 2, unlist) %>% 
    as_tibble %>% 
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means = mean, sds = sd)) %>% t() 
}
names(CCIunc) <- names(CCIres2)


## combine the uncertainties
unc_ranB <- bind_rows(CCD = CCDunc, FCI = FCIunc, CCI = CCIunc, .id="id") %>% 
  # n = repeat N: 3 algo * 2 (mean & sd), repeat c(means,sd) by length of N * 3 algo
  mutate(n = rep(N,6), statistics = rep(c("means", "sds"), each = length(N), times = 3)) %>% 
  tidyr::pivot_longer(cols = -c(id, n, statistics), names_to = "condition", values_to = "value") %>% 
  mutate(
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LV", "without LV"),
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
shdplot_ranB <- SHD_ranB %>%
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(n, levels = c("50", "150", "500", "1000", "1500", "2000","3000", "4000", "5000", "10000")), y=means, group = id, colour = id, fill = id)) +
  geom_line(aes(group = id)) +
  geom_point(size=1) + 
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by 3)
  geom_ribbon(aes(ymin=means-qnorm(0.975)*sds/sqrt(n)*3, ymax=means+qnorm(0.975)*sds/sqrt(n)*3), alpha=0.2, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x="", y="", title = "") +
  theme_minimal() +
  MyTheme + 
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  scales = "free_y", switch="y") +
  ggtitle("(a) SHD") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


## plot precision
precisionplot_ranB <- prec_ranB %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(n, levels = c("50", "150", "500", "1000", "1500", "2000", "3000", "4000", "5000", "10000")), y=means, group = id, colour = id, fill=id)) +
  geom_line(aes(group = id)) +
  geom_point(size=1) +
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by 2)
  geom_ribbon(aes(ymin=means - qnorm(0.975)*sds/sqrt(n)*2, ymax=means + qnorm(0.975)*sds/sqrt(n)*2), alpha=0.2, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  theme_minimal() +
  MyTheme + 
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  switch="y") +
  labs(title = "(b) Precision", x = "", y = "") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

## plot recall
recallplot_ranB <- rec_ranB %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(n, levels = c("50", "150", "500", "1000", "1500", "2000", "3000", "4000", "5000", "10000")), y=means, group = id, colour = id, fill=id)) +
  geom_line(aes(group = id)) +
  geom_point(size=1) +
  geom_ribbon(aes(ymin=means - qnorm(0.975)*sds/sqrt(n)*2, ymax=means + qnorm(0.975)*sds/sqrt(n)*2), alpha=0.2, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  theme_minimal() +
  MyTheme + 
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  switch="y") +
  labs(title = "(c) Recall", x = "", y = "") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


## plot uncertainty 
uncertaintyplot_ranB <- unc_ranB %>%  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(n, levels = c("50", "150", "500", "1000", "1500", "2000",  "3000", "4000", "5000", "10000")), y=means, group = id, colour = id, fill = id)) +
  geom_line(aes(group = id)) +
  geom_point(size=1) + 
  #geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), width=0.1) +
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by )
  geom_ribbon(aes(ymin=means-qnorm(0.975)*sds/sqrt(as.numeric(N)), ymax=means+qnorm(0.975)*sds/sqrt(as.numeric(N))), alpha=0.2, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x="N", y="", title = "") +
  theme_minimal() +
  MyTheme + 
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  scales = "free_y", switch="y") +
  ggtitle("(d) Uncertainty") 
# theme(axis.title.x=element_blank(),
#       axis.text.x=element_blank(),
#       axis.ticks.x=element_blank())


# combine the plots
ggpubr::ggarrange(shdplot_ranB, precisionplot_ranB, recallplot_ranB, uncertaintyplot_ranB, nrow=4, common.legend = TRUE, legend = "bottom")
#ggsave(filename = "results/samplingbeta_result.pdf", width = 25, height = 35, dpi = 300, units = "cm")

ggpubr::ggarrange(shdplot_ranB, precisionplot_ranB, nrow=2, common.legend = TRUE, legend = "bottom")
#ggsave(filename = "results/samplingbeta_result1.pdf", width = 25, height = 35, dpi = 300, units = "cm")
ggpubr::ggarrange(recallplot_ranB, uncertaintyplot_ranB, nrow=2, common.legend = TRUE, legend = "bottom")
#ggsave(filename = "results/samplingbeta_result2.pdf", width = 25, height = 35, dpi = 300, units = "cm")
