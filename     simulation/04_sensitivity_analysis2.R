## =============================================================================
## Description
# 
# This script contains the code to perform the secondary analysis with randomly 
# sampled coefficients for the regression matrix B with "only positive values".
#
# Purpose: to examine the impact of the specified weights on the results of the 
# main analysis, and to determine if there are any differences when randomly 
# sampling B weights from both negative and positive values.
#
# The remaining setup is identical to that of `03_sensitivity_analysis1.R`.
## =============================================================================
# The content is as follows.
# 0. Preparation: Source and load necessary functions and packages.
#
# 1. Simulate data: Generate data based on randomly sampled B matrix (with only
#    positive values) at every iteration. 
#
# 2. Run algorithms: Apply three algorithms -- CCD, FCI, and CCI -- to each 
#    of the simulated datasets then estimate PAGs.
#
# 3. Evaluate performance: Compute structural Hamming distance, precision, 
#    recall, and uncertainty rate for each of the estimated PAGs.
#
# 4. Create figures: Plot figures for each evaluation metric comparing the 
#    performance of each algorithm per condition.
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
library(ggh4x)

## set the seed
set.seed(123)


## =============================================================================
## 1. Generate data 
## =============================================================================

# specify the sample sizes
N <- c(50, 150, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000)
# specify replication number
n <- 500
# specify alpha level
alpha <- 0.01


## create all simulated data using random B 
simdata_woLV <- list(B5sparse = B5sparse, B5dense = B5dense, 
                     B10sparse =  B10sparse, B10dense = B10dense) %>% 
  map(~
        sampleranB2_pos(.x)
  )

simdata_5pwLV <- list(B5_lvsparse = B5_lvsparse, B5_lvdense = B5_lvdense) %>% 
  map(~
        sampleranB2_pos(.x, LV = 6)
  )

simdata_10pwLV <- list(B10_lvsparse = B10_lvsparse, B10_lvdense = B10_lvdense) %>% 
  map(~
        sampleranB2_pos(.x, LV = c(11, 12))
  )

# append all datasets in a single list
simdatalist <- append(simdata_woLV, append(simdata_5pwLV, simdata_10pwLV))


## =============================================================================
## 2. Run algorithms
## =============================================================================
# run CCD
CCDres_pos <- simdatalist %>% 
  map_depth(3, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)) %>% 
  map_depth(3, ~CreateAdjMat(.x, length(.x$nodes)))

# run FCI
FCIres_pos <- simdatalist %>% 
  map_depth(3, ~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
                     alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
                     labels = colnames(.x)) %>% .@amat 
  )

# run CCI
CCIres_pos <- simdatalist %>% 
  map_depth(3, ~ cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=alpha, 
                     labels = colnames(.x), p = ncol(.x)) %>% .$maag 
  )

# append all true adj.matrices into a single list
truemods <- list(trueag_5psparse, trueag_5pdense, trueag_10psparse, 
                 trueag_10pdense, trueag_5psparseLV, trueag_5pdenseLV, 
                 trueag_10psparseLV, trueag_10pdenseLV)

# load results
load("    simulation/output/randomB_pos/CCDres3_pos.Rdata") # CCDres
load("    simulation/output/randomB_pos/FCIres3_pos.Rdata") # FCIres
load("    simulation/output/randomB_pos/CCIres3_pos.Rdata") # CCIres


## =============================================================================
## 3. Evaluate performance
## =============================================================================

## compute SHD
# SHD values for CCD
CCDshd <- list()
for(i in 1:length(CCDres)){
  CCDshd[[i]] <- CCDres[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N))  %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCDshd) <- names(CCDres)

# SHD values for FCI
FCIshd <- list()
for(i in 1:length(FCIres)){
  FCIshd[[i]] <- FCIres[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(FCIshd) <- names(FCIres)

# SHD values for CCI
CCIshd <- list()
for(i in 1:length(CCIres)){
  CCIshd[[i]] <- CCIres[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCIshd) <- names(CCIres)

## combine the SHD values
SHD_ranB <- bind_rows(CCD = CCDshd, FCI = FCIshd, CCI = CCIshd, .id="id") %>% 
  # convert it to a long format
  tidyr::pivot_longer(cols = -c(id), names_to = "condition", values_to = "value") %>% 
  mutate(
    # create variables (conditions)
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LC", "without LC"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>%  
  # unnest a list-column
  tidyr::unnest(value) %>% 
  # convert it to a long format
  tidyr::pivot_longer(cols = starts_with("N ="), names_to = "n", values_to = "value") %>% 
  # create variables (statistics & sample size)
  mutate(statistics = stringr::str_split(n, "_", simplify=T)[,2],
         n = as.numeric(stringr::str_extract_all(n, "[0-9]+")))  %>% 
  # bring the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric))



## compute precision
# precision for CCD
CCDprec <- list()
for(i in 1:length(CCDres)){
  CCDprec[[i]] <- CCDres[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCDprec) <- names(CCDres)

# precision for FCI
FCIprec <- list()
for(i in 1:length(FCIres)){
  FCIprec[[i]] <- FCIres[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(FCIprec) <- names(FCIres)

# precision for CCI
CCIprec <- list()
for(i in 1:length(CCIres)){
  CCIprec[[i]] <- CCIres[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCIprec) <- names(CCIres)

## combine the precision values
prec_ranB <- bind_rows(CCD = CCDprec, FCI = FCIprec, CCI = CCIprec, .id="id") %>% 
  # convert it to a long format
  tidyr::pivot_longer(cols = -c(id), names_to = "condition", values_to = "value") %>% 
  mutate(
    # create variables (conditions)
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LC", "without LC"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>%  
  # unnest a list-column
  tidyr::unnest(value) %>% 
  # convert it to a long format
  tidyr::pivot_longer(cols = starts_with("N ="), names_to = "n", values_to = "value") %>% 
  # create variables (statistics & sample size)
  mutate(statistics = stringr::str_split(n, "_", simplify=T)[,2],
         n = as.numeric(stringr::str_extract_all(n, "[0-9]+")))  %>% 
  # bring the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric))


## compute recall
# recall for CCD
CCDrec <- list()
for(i in 1:length(CCDres)){
  CCDrec[[i]] <- CCDres[[i]] %>% 
    map_depth(2, ~recall2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCDrec) <- names(CCDres)

# recall for FCI
FCIrec <- list()
for(i in 1:length(FCIres)){
  FCIrec[[i]] <- FCIres[[i]] %>% 
    map_depth(2, ~ recall2(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(FCIrec) <- names(FCIres)

# recall for CCI
CCIrec <- list()
for(i in 1:length(CCIres)){
  CCIrec[[i]] <- CCIres[[i]] %>% 
    map_depth(2, ~ recall2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%
    rename_with(~ paste0("N = ", N)) %>% summarize_all(list(means=mean, sds=sd))
}
names(CCIrec) <- names(CCIres)

## combine the recall values
rec_ranB <- bind_rows(CCD = CCDrec, FCI = FCIrec, CCI = CCIrec, .id="id") %>% 
  # convert it to a long format
  tidyr::pivot_longer(cols = -c(id), names_to = "condition", values_to = "value") %>% 
  mutate(
    # create variables (conditions)
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LC", "without LC"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>%  
  # unnest a list-column
  tidyr::unnest(value) %>% 
  # convert it to a long format
  tidyr::pivot_longer(cols = starts_with("N ="), names_to = "n", values_to = "value") %>% 
  # create variables (statistics & sample size)
  mutate(statistics = stringr::str_split(n, "_", simplify=T)[,2],
         n = as.numeric(stringr::str_extract_all(n, "[0-9]+")))  %>% 
  # bring the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric))


## compute uncertainty
# uncertainty for CCD
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

# uncertainty for FCI
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

# uncertainty for CCI
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


## combine the uncertainty values
unc_ranB <- bind_rows(CCD = CCDunc, FCI = FCIunc, CCI = CCIunc, .id="id") %>% 
  # n = repeat N: 3 algo * 2 (mean & sd), repeat c(means,sd) by length of N * 3 algo
  mutate(n = rep(N,6), statistics = rep(c("means", "sds"), each = length(N), times = 3)) %>% 
  # convert it to a long format
  tidyr::pivot_longer(cols = -c(id, n, statistics), names_to = "condition", values_to = "value") %>% 
  mutate(
    # create variables (conditions)
    netsize = paste0(stringr::str_match_all(condition, "[0-9]+"), "p"),
    latentvar = ifelse(stringr::str_detect(condition, "lv")==TRUE, "with LC", "without LC"),
    densities = ifelse(stringr::str_detect(condition, "dense")==TRUE, "dense", "sparse") 
  ) %>% 
  # bring the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric))



## =============================================================================
## 4. Plot the figures for comparison
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


## plot SHD figure
shdplot_ranB_pos <- SHD_ranB %>%
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  # create a ggplot object
  ggplot(aes(x= as.numeric(n), y=means, group = id, colour = id, fill = id)) +
  # add line graphs
  geom_line(aes(group = id)) +
  # add scattered points
  geom_point(size=1) + 
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="", y="", title = "") +
  # apply themes
  theme_minimal() +
  MyTheme + 
  # create facets
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ 
                        factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")), scales = "free_y", switch="y") +
  # specify custom breaks
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  ggtitle("(a) SHD")  +
  guides(color = "none", fill = "none")
# save the plot
# ggsave(shdplot_ranB_pos, filename = "results/samplingbeta_pos_shd.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")



## plot precision figure
precisionplot_ranB_pos <- prec_ranB %>% 
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  # create a ggplot object
  ggplot(aes(x= as.numeric(n), y=means, group = id, colour = id, fill=id)) +
  # add line graphs
  geom_line(aes(group = id)) +
  # add scattered points
  geom_point(size=1) +
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  # apply themes
  theme_minimal() +
  MyTheme + 
  # create facets
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ 
                        factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")), switch="y") +
  # specify custom breaks
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  labs(title = "(b) Precision", x = "", y = "") +
  guides(color = "none", fill = "none")
# save the plot
# ggsave(precisionplot_ranB_pos, filename = "results/samplingbeta_pos_prec.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")


## plot recall figure
recallplot_ranB_pos <- rec_ranB %>% 
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  # create a ggplot object
  ggplot(aes(x= as.numeric(n), y=means, group = id, colour = id, fill=id)) +
  # add line graphs
  geom_line(aes(group = id)) +
  # add scattered points
  geom_point(size=1) +
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  # apply themes
  theme_minimal() +
  MyTheme + 
  # create facets
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ 
                        factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")), switch="y") +
  # specify custom breaks
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  labs(title = "(c) Recall", x = "", y = "")+
  guides(color = "none", fill = "none")
# save the plot
# ggsave(recallplot_ranB_pos, filename = "results/samplingbeta_pos_recall.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")

## plot uncertainty figure
uncertaintyplot_ranB_pos <- unc_ranB %>%  
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  # create a ggplot object
  ggplot(aes(x= as.numeric(n), y=means, group = id, colour = id, fill = id)) +
  # add line graphs
  geom_line(aes(group = id)) +
  # add scattered points
  geom_point(size=1) + 
  # add interquartile range (IQR)
  geom_ribbon(aes(ymin=means+qnorm(0.25)*sds, ymax=means+qnorm(0.75)*sds), alpha=0.15, color=NA) +
  # specify custom colors
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="N", y="", title = "") +
  # apply themes
  theme_minimal() +
  MyTheme + 
  # create facets
  ggh4x::facet_nested(factor(netsize, levels = c("5p", "10p"), 
                             labels=c("p = 5", "p = 10")) ~ 
                        factor(latentvar, levels = c("without LC", "with LC")) + 
                        factor(densities, levels=c("sparse", "dense")),  scales = "free_y", switch="y") +
  # specify custom breaks
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  ggtitle("(d) Uncertainty") 
# save the plot
# ggsave(uncertaintyplot_ranB_pos, filename = "results/samplingbeta_pos_unc.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")


## combine the plots
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
