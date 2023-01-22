source("../code/simulation_code.R")
library(dplyr)
library(purrr)

## Precision & Recall
# specify rownames
# rows <- c(apply(expand.grid(c("ccd", "fci", "cci"), c("5p", "10p"), c("sparse", "dense")),1,paste, collapse="."),  apply(expand.grid(c("ccd", "fci", "cci"), c("5p_LV", "10p_LV")),1,paste, collapse="."))

# compute the average precision and sd
results <- list(res_ccd5psparse, res_fci5psparse, res_cci5psparse, res_ccd10psparse, res_fci10psparse, res_cci10psparse, res_ccd5pdense, res_fci5pdense, res_cci5pdense, res_ccd10pdense, res_fci10pdense, res_cci10pdense, res_ccd5pLV2, res_fci5pLV2, res_cci5pLV2, res_ccd10pLV,  res_fci10pLV, res_cci10pLV) %>% 
  map(~
  # think about how to deal with NAs or do I want to define sth. else instead of NAs.
  na.omit(.x) %>% 
  summarise(across(everything(.x), list(mean = mean, sd = sd)))
) %>% bind_rows() %>% 
  mutate(algorithm = rep(c("ccd", "fci", "cci"), 6),
         condition = rep(c("5p_sparse", "10p_sparse", "5p_dense", "10p_dense", "5p_LV", "10p_LV"), each=3)) %>% 
  # brings the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric))

## Uncertainty
uncertainties <- data.frame(uncer_ccd5psparse, uncer_fci5psparse, uncer_cci5psparse, uncer_ccd10psparse, uncer_fci10psparse, uncer_cci10psparse, uncer_ccd5pdense, uncer_fci5pdense, uncer_cci5pdense, uncer_ccd10pdense, uncer_fci10pdense, uncer_cci10pdense, uncer_ccd5pLV2, uncer_fci5pLV2, uncer_cci5pLV2, uncer_ccd10pLV, uncer_fci10pLV, uncer_cci10pLV) %>% 
  summarise(across(everything(), list(mean = mean, sd = sd))) %>% tidyr::pivot_longer(everything()) %>%
  mutate(algorithm = substr(stringr::str_split(name, "_", simplify = T)[,2], 1, 3),
         condition = substring(stringr::str_split(name, "_", simplify = T)[,2], 4, last = 1000000L), 
         statistics = stringr::str_split(name, "_", simplify = T)[,3]
  )%>% dplyr::select(-name) %>% tidyr::pivot_wider(names_from= statistics, values_from=value)


## SHD
SHDs <- data.frame(SHD_ccd5psparse, SHD_fci5psparse, SHD_cci5psparse, SHD_ccd10psparse, SHD_fci10psparse, SHD_cci10psparse, SHD_ccd5pdense, SHD_fci5pdense, SHD_cci5pdense, SHD_ccd10pdense, SHD_fci10pdense, SHD_cci10pdense, SHD_ccd5pLV2, SHD_fci5pLV2, SHD_cci5pLV2, SHD_ccd10pLV, SHD_fci10pLV, SHD_cci10pLV) %>% 
  summarise(across(everything(), list(mean = mean, sd = sd))) %>%  tidyr::pivot_longer(cols = everything()) %>%
  mutate(algorithm = substr(stringr::str_split(name, "_", simplify = T)[,2], 1, 3),
         condition = substring(stringr::str_split(name, "_", simplify = T)[,2], 4, last = 1000000L), 
         statistics = stringr::str_split(name, "_", simplify = T)[,3]
  ) %>% dplyr::select(-name) %>% tidyr::pivot_wider(names_from= statistics, values_from=value)



## Plot the results
library(ggplot2)
## Specify my custom theme
MyTheme <-  theme(plot.title = element_blank(),
                  plot.subtitle = element_text( face = "italic"),
                  axis.text=element_text(face = "bold"),
                  legend.text = element_text(face = "bold"))

## ========================
## WITHOUT LV CONDITION
## ========================

## WITHOUT LV CONDITION: PRECISION
results %>% 
  # exclude LV conditions
  filter(!grepl("LV", condition)) %>% 
ggplot(aes(x= factor(condition, levels = c("5p_sparse", "5p_dense", "10p_sparse", "10p_dense")), y=average_precision_mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(1000), ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(1000)), width=0.1) +
  labs(x="", y="", title = "", subtitle = "Without Latent Variable_PRECISION") +
  theme_classic() + MyTheme

## WITHOUT LV CONDITION: RECALL
results %>% 
  # exclude LV conditions
  filter(!grepl("LV", condition)) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p_sparse", "5p_dense", "10p_sparse", "10p_dense")), y=average_recall_mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(1000), ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(1000)), width=0.1) +
  labs(x="", y="", title = "", subtitle = "Without Latent Variable_RECALL") +
  theme_classic() + MyTheme

## WITHOUT LV CONDITION: UNCERTAINTY
uncertainties %>% 
  # exclude LV conditions
  filter(!grepl("LV", condition)) %>% 
  ggplot(aes(x= factor(condition, levels = c("5psparse", "5pdense", "10psparse", "10pdense")), y=mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(1000), ymax=mean+qnorm(0.975)*sd/sqrt(1000)), width=0.1) +
  labs(x="", y="", title = "", subtitle = "Without Latent Variable_UNCERTAINTY") +
  theme_classic() + MyTheme

## WITHOUT LV CONDITION: SHD
SHDs %>% 
  # exclude LV conditions
  filter(!grepl("LV", condition)) %>% 
  ggplot(aes(x= factor(condition, levels = c("5psparse", "5pdense", "10psparse", "10pdense")), y=mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000)) 
  geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(1000), ymax=mean+qnorm(0.975)*sd/sqrt(1000)), width=0.1) +
  labs(x="", y="", title = "", subtitle = "Without Latent Variable_SHD") +
  theme_classic() + MyTheme


## ========================
## WITH LV CONDITION
## ========================

## WITH LV CONDITION: PRECISION
results %>% 
  # exclude LV conditions
  filter(grepl("LV", condition)) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p_LV", "10p_LV")), y=average_precision_mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(1000), ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(1000)), width=0.1) +
  labs(x="", y="", title = "", subtitle = "With a Latent Variable_PRECISION") +
  theme_classic() + MyTheme

## WITH LV CONDITION: RECALL
results %>% 
  # exclude LV conditions
  filter(grepl("LV", condition)) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p_LV", "10p_LV")), y=average_recall_mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(1000), ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(1000)), width=0.1) +
  labs(x="", y="", title = "", subtitle = "With a Latent Variable_RECALL") +
  theme_classic() + MyTheme

## WITH LV CONDITION: UNCERTAINTY
uncertainties %>% 
  # exclude LV conditions
  filter(grepl("LV", condition)) %>% 
  ggplot(aes(x= factor(condition, levels = c("5pLV2", "10pLV")), y=mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(1000), ymax=mean+qnorm(0.975)*sd/sqrt(1000)), width=0.1) +
  labs(x="", y="", title = "", subtitle = "With a Latent Variable_UNCERTAINTY") +
  theme_classic() + MyTheme

## WITH LV CONDITION: SHD
SHDs %>% 
  # exclude LV conditions
  filter(grepl("LV", condition)) %>% 
  ggplot(aes(x= factor(condition, levels = c("5pLV2", "10pLV")), y=mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(1000), ymax=mean+qnorm(0.975)*sd/sqrt(1000)), width=0.1) +
  labs(x="", y="", title = "", subtitle = "With a Latent Variable_SHD") +
  theme_classic() + MyTheme




