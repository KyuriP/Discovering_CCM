## =============================================================================
## Description
# 
# This script contains code for the secondary analysis with 
# randomly sampled coefficients for the regression matrix B.
#
# Purpose: to determine the extent to which the results of the main analysis 
# depend on the specific weights that were specified. 
#
# As is the case with the main analysis, there are in total 8 models considered
# and 500 datasets are generated from each model.
## =============================================================================
# The content is as follows.
# 0. Preparation: Source and load necessary functions and packages. 
#                  
# 1. Simulate data: Generate data based on randomly sampled B matrix 
#    at every iteration. 
#
# 2. Run algorithms: Apply three algorithms -- CCD, FCI, and CCI -- to each 
#    of the simulated datasets then estimate PAGs.
#    *NOTE*: the code for running the algorithms is currently commented out due 
#    to the long running time. If you are interested in running the algorithms, 
#    simply uncomment the relevant code lines and execute the code.
#
# 3. Evaluate performance: Compute structural Hamming distance, precision, 
#    recall, and uncertainty rate for each of the estimated PAGs.
#
# 4. Create figures: Plot figures for each evaluation metric comparing the 
#    performance of each algorithm per condition.
#
# 5. Create an extra sub-figure: Show the performance under 5p dense 
#    without a latent confounder condition (Figure 18 in the paper).
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

# set the seed
set.seed(123)



## =============================================================================
## 1. Simulate data 
## =============================================================================

# specify the sample sizes
N <- c(50, 150, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000)
# specify replication number
n <- 500
# specify alpha level
alpha <- 0.01

## simulate datasets using randomly sampled B matrices
simdata_woLV <- list(B5sparse = B5sparse, B5dense = B5dense, 
                     B10sparse = B10sparse, B10dense = B10dense) %>% 
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

# append all datasets in a single list
simdatalist <- append(simdata_woLV, append(simdata_5pwLV, simdata_10pwLV))



## =============================================================================
## 2. Run algorithms
## =============================================================================
# run CCD 
# CCDres <- simdatalist %>% 
#   map_depth(3, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)) %>% 
#   map_depth(3, ~CreateAdjMat(.x, length(.x$nodes)))
            
# run FCI
# FCIres <- simdatalist %>% 
#   map_depth(3, ~ fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
#                                    alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, 
#                      labels = colnames(.x)) %>% .@amat 
#   )

# run CCI
# CCIres <- simdatalist %>% 
#   map_depth(3, ~ cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=alpha, 
#                      labels = colnames(.x), p = ncol(.x)) %>% .$maag 
#   )

# append all true adj.matrices into a single list
truemods <- list(trueag_5psparse, trueag_5pdense, trueag_10psparse, 
                 trueag_10pdense, trueag_5psparseLV, trueag_5pdenseLV, 
                 trueag_10psparseLV, trueag_10pdenseLV)

# save results
# save(CCDres, file = "simulation/output/randomB_n500/CCDres2_randomB.Rdata")
# save(FCIres, file = "simulation/output/randomB_n500/FCIres2_randomB.Rdata")
# save(CCIres, file = "simulation/output/randomB_n500/CCIres2_randomB.Rdata")

# load results
load("simulation/output/randomB_n500/CCDres2_randomB.Rdata")
load("simulation/output/randomB_n500/FCIres2_randomB.Rdata")
load("simulation/output/randomB_n500/CCIres2_randomB.Rdata")



## =============================================================================
## 3. Evaluate performance
## =============================================================================

## compute SHD
# SHD values for CCD
CCDshd <- list()
for(i in 1:length(CCDres2)){
  CCDshd[[i]] <- CCDres2[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N))  %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCDshd) <- names(CCDres2)

# SHD values for FCI
FCIshd <- list()
for(i in 1:length(FCIres2)){
  FCIshd[[i]] <- FCIres2[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(FCIshd) <- names(FCIres2)

# SHD values for CCI
CCIshd <- list()
for(i in 1:length(CCIres2)){
  CCIshd[[i]] <- CCIres2[[i]] %>% 
    map_depth(2, ~SHD(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>% apply(., 2, unlist) %>%  
    as.data.frame %>% rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCIshd) <- names(CCIres2)

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
for(i in 1:length(CCDres2)){
  CCDprec[[i]] <- CCDres2[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCDprec) <- names(CCDres2)

# precision for FCI
FCIprec <- list()
for(i in 1:length(FCIres2)){
  FCIprec[[i]] <- FCIres2[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(FCIprec) <- names(FCIres2)

# precision for CCI
CCIprec <- list()
for(i in 1:length(CCIres2)){
  CCIprec[[i]] <- CCIres2[[i]] %>% 
    map_depth(2, ~precision2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCIprec) <- names(CCIres2)

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
for(i in 1:length(CCDres2)){
  CCDrec[[i]] <- CCDres2[[i]] %>% 
    map_depth(2, ~recall2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCDrec) <- names(CCDres2)

# recall for FCI
FCIrec <- list()
for(i in 1:length(FCIres2)){
  FCIrec[[i]] <- FCIres2[[i]] %>% 
    map_depth(2, ~ recall2(truemods[[i]], .x))  %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%     
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(FCIrec) <- names(FCIres2)

# recall for CCI
CCIrec <- list()
for(i in 1:length(CCIres2)){
  CCIrec[[i]] <- CCIres2[[i]] %>% 
    map_depth(2, ~ recall2(truemods[[i]], .x)) %>% 
    do.call("cbind", .) %>%  apply(., 2, unlist) %>% as.data.frame()  %>%
    rename_with(~ paste0("N = ", N)) %>% 
    summarize_all(list(means=mean, sds=sd))
}
names(CCIrec) <- names(CCIres2)

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

# uncertainty for FCI
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

# uncertainty for CCI
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
shdplot_ranB <- SHD_ranB %>%
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
## save the plot
# ggsave(filename = "figures/samplingbeta_shd.pdf", width = 25, height = 10, dpi = 300, units = "cm")


## plot precision figure
precisionplot_ranB <- prec_ranB %>% 
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
# ggsave(filename = "figures/samplingbeta_prec.pdf", width = 25, height = 10, dpi = 300, units = "cm")


## plot recall figure
recallplot_ranB <- rec_ranB %>% 
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
# ggsave(filename = "figures/samplingbeta_rec.pdf", width = 25, height = 10, dpi = 300, units = "cm")


## plot uncertainty figure
uncertaintyplot_ranB <- unc_ranB %>%  
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
# ggsave(filename = "figures/samplingbeta_unc.pdf", width = 25, height = 10.5, dpi = 300, units = "cm")


## combine the figures
ggpubr::ggarrange(shdplot_ranB, precisionplot_ranB, recallplot_ranB, uncertaintyplot_ranB, nrow=4, common.legend = TRUE, legend = "bottom")
#ggsave(filename = "figures/samplingbeta_result.pdf", width = 25, height = 35, dpi = 300, units = "cm")
ggpubr::ggarrange(shdplot_ranB, precisionplot_ranB, nrow=2, common.legend = TRUE, legend = "bottom")
# ggsave(filename = "figures/samplingbeta_result1.pdf", width = 25, height = 22, dpi = 300, units = "cm")
ggpubr::ggarrange(recallplot_ranB, uncertaintyplot_ranB, nrow=2, common.legend = TRUE, legend = "bottom")
# ggsave(filename = "figures/samplingbeta_result2.pdf", width = 25, height = 22, dpi = 300, units = "cm")



## =============================================================================
## 5. Extract only 5p dense cases (Figure 18)
## =============================================================================
# shd plot
dense5pshd <- SHD_ranB |> filter(condition == "B5dense") |>
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) |>
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
  # specify custom breaks
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  # scale_x_continuous(breaks=c(seq(50, 10000, by = 1000),10000)) +
  labs(x="N", y="", title = "") +
  # apply themes
  theme_minimal() +
  MyTheme + 
  ggtitle("(a) SHD") 

# precision plot
dense5pprec <-prec_ranB |> filter(condition == "B5dense") |>
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) |>
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
  # specify custom breaks
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  labs(x="N", y="", title = "") +
  # apply themes
  theme_minimal() +
  MyTheme + 
  ggtitle("(b) Precision") 

# recall plot
dense5prec <- rec_ranB |> filter(condition == "B5dense") |>
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) |>
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
  # specify custom breaks
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  labs(x="N", y="", title = "") +
  # apply themes
  theme_minimal() +
  MyTheme + 
  ggtitle("(c) Recall") 

# uncertainty plot
dense5punc <- unc_ranB |> filter(condition == "B5dense") |>
  # convert it to a wide format
  tidyr::pivot_wider(names_from = statistics, values_from=value) |>
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
  # specify custom breaks
  scale_x_continuous(breaks=c(50, 2500, 5000, 7500, 10000)) +
  labs(x="N", y="", title = "") +
  # apply themes
  theme_minimal() +
  MyTheme + 
  ggtitle("(d) Uncertainty") 


## combine the plots
ggpubr::ggarrange(dense5pshd, dense5pprec, dense5prec, dense5punc, ncol=2, nrow=2, common.legend = TRUE, legend = "bottom")
# save the plot
# ggsave(filename = "figures/samplingbeta_dense5p.pdf", width = 17, height = 13, dpi = 300, units = "cm")
# ggsave(filename = "figures/samplingbeta_dense5p_label10.pdf", width = 17, height = 13, dpi = 300, units = "cm")
