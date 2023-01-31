source("code/simulation_code.R")
library(dplyr)
library(purrr)
library(ggplot2)

## Precision & Recall

# compute the average precision and sd
results <- list(res_ccd5psparse, res_fci5psparse, res_cci5psparse, res_ccd10psparse, res_fci10psparse, res_cci10psparse, res_ccd5pdense, res_fci5pdense, res_cci5pdense, res_ccd10pdense, res_fci10pdense, res_cci10pdense, res_ccd5pLV2, res_fci5pLV2, res_cci5pLV2, res_ccd10pLV,  res_fci10pLV, res_cci10pLV) %>% 
  # transpose df
  map(~ sjmisc::rotate_df(.x) %>%
  # add sample size (N) info
  rename_with(~paste0(.x, "N = ", rep(N, each=8)))  %>%
  # think about how to deal with NAs or do I want to define sth. else instead of NAs.
  #na.omit(.x) %>% 
  summarise(across(everything(.), list(mean = ~mean(., na.rm=T), sd = ~sd(., na.rm=T))))) %>% 
  bind_rows() %>% 
  mutate(algorithm = rep(c("ccd", "fci", "cci"), 6),
         condition = rep(c("5p_sparse", "10p_sparse", "5p_dense", "10p_dense", "5p_LV", "10p_LV"), each=3)) %>%
  # brings the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric)) %>% 
  # convert itto a long format
  tidyr::pivot_longer(!c(algorithm, condition), names_to = "metric", values_to = "value") %>% 
  # Add sample size column (N) & clean up the column name
  mutate(N = stringr::str_extract(metric, "(?<=[N =])\\d+"),
         metric = stringr::str_replace_all(metric, "[0-9.]+|[N =]", "")) 

## Uncertainty
uncertainties <- bind_rows("ccd_5p-sparse" = uncer_ccd5psparse, "fci_5p-sparse" = uncer_fci5psparse, "cci_5p-sparse"=uncer_cci5psparse, "ccd_10p-sparse"=uncer_ccd10psparse, "fci_10p-sparse" = uncer_fci10psparse, "cci_10p-sparse" = uncer_cci10psparse, "ccd_5p-dense"=uncer_ccd5pdense, "fci_5p-dense"=uncer_fci5pdense, "cci_5p-dense"=uncer_cci5pdense, "ccd_10p-dense"=uncer_ccd10pdense, "fci_10p-dense"=uncer_fci10pdense, "cci_10p-dense"=uncer_cci10pdense, "ccd_5p-LV"=uncer_ccd5pLV2, "fci_5p-LV"=uncer_fci5pLV2, "cci_5p-LV"=uncer_cci5pLV2, "ccd_10p-LV"=uncer_ccd10pLV, "fci_10p-LV"=uncer_fci10pLV, "cci_10p-LV"=uncer_cci10pLV, .id="id") %>% 
  group_by(id) %>% 
  summarise_all(list(mean = mean, sd = sd)) %>%  
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2]) %>% 
  tidyr::pivot_longer(!c(algorithm, condition, id), names_to = "name", values_to = "value") %>% 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  relocate(where(is.character), .before = where(is.numeric))



## SHD
SHDs <- bind_rows("ccd_5p-sparse" = SHD_ccd5psparse, "fci_5p-sparse" = SHD_fci5psparse, "cci_5p-sparse"=SHD_cci5psparse, "ccd_10p-sparse"= SHD_ccd10psparse, "fci_10p-sparse" = SHD_fci10psparse, "cci_10p-sparse" = SHD_cci10psparse, "ccd_5p-dense"= SHD_ccd5pdense, "fci_5p-dense"=SHD_fci5pdense, "cci_5p-dense"=SHD_cci5pdense, "ccd_10p-dense"= SHD_ccd10pdense, "fci_10p-dense"=SHD_fci10pdense, "cci_10p-dense"=SHD_cci10pdense, "ccd_5p-LV"=SHD_ccd5pLV2, "fci_5p-LV"=SHD_fci5pLV2, "cci_5p-LV"=SHD_cci5pLV2, "ccd_10p-LV"=SHD_ccd10pLV, "fci_10p-LV"=SHD_fci10pLV, "cci_10p-LV"=SHD_cci10pLV, .id="id") %>% 
  group_by(id) %>% 
  summarise_all(list(mean = mean, sd = sd)) %>%  
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2]) %>% 
  tidyr::pivot_longer(!c(algorithm, condition, id), names_to = "name", values_to = "value") %>% 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  relocate(where(is.character), .before = where(is.numeric)) 
  
  


## Plot the results
library(ggplot2)
## Specify my custom theme
MyTheme <-  theme(plot.title = element_blank(),
                  plot.subtitle = element_text( face = "italic"),
                  axis.text=element_text(face = "bold"),
                  legend.text = element_text(face = "bold"))

## ========================
## precision plots
## ========================
precision_plots <- c(unique(results$condition)) %>% 
  map(~
results %>% 
  filter(condition == .x & grepl("average_precision", metric)) %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "5000")), y=average_precision_mean, group = algorithm, colour = algorithm, fill=algorithm)) +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  #geom_errorbar(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N)), ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))), width=0.1) +
  geom_ribbon(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N)), ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))), alpha=0.2,lwd=0) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="N", y="", title = "", subtitle = .x) +
  theme_classic() + MyTheme
)

## ========================
## recall plots
## ========================
recall_plots <- c(unique(results$condition)) %>% 
  map(~
results %>% 
  filter(condition == .x & grepl("average_recall", metric)) %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "5000")), y=average_recall_mean, group = algorithm, colour = algorithm, fill= algorithm)) +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  #geom_errorbar(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N)), ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))), width=0.1) +
  geom_ribbon(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N)), ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))), alpha=0.2,lwd=0) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="N", y="", title = "", subtitle = .x) +
  theme_classic() + MyTheme
)

## ========================
## uncertainty plots
## ========================
## uncertainty (uncertainty rate of fci and cci are exactly the same!)
uncertainty_plots <- c(unique(uncertainties$condition)) %>% 
  map(~
uncertainties %>%
  filter(condition == .x) %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "5000")), y=mean, group = algorithm, colour = algorithm, fill=algorithm)) +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  #geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), width=0.1) +
  geom_ribbon(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), alpha=0.2,lwd=0) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="N", y="", title = "", subtitle = paste0("UNCERTAINTY -", .x)) +
  theme_classic() + MyTheme
)

## ========================
## SHD plots
## ========================
SHD_plots <- c(unique(SHDs$condition)) %>% 
  map(~
SHDs %>%
  filter(condition == .x) %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "5000")), y=mean, group = algorithm, colour = algorithm, fill = algorithm)) +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  #geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), width=0.1) +
  geom_ribbon(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), alpha=0.2,lwd=0) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(x="N", y="", title = "", paste0("SHD -", .x)) +
  theme_classic() + MyTheme
)
################################################################################
## WITHOUT LV CONDITION: PRECISION
results %>% 
  # exclude LV conditions
  filter(!grepl("LV", condition), N == "1000") %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
ggplot(aes(x= factor(condition, levels = c("5p_sparse", "5p_dense", "10p_sparse", "10p_dense")), y=average_precision_mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N)), ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))), width=0.1) +
  labs(x="", y="", title = "", subtitle = "Without Latent Variable_PRECISION") +
  theme_classic() + MyTheme

## WITHOUT LV CONDITION: RECALL
results %>% 
  # exclude LV conditions
  filter(!grepl("LV", condition), N == "1000") %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p_sparse", "5p_dense", "10p_sparse", "10p_dense")), y=average_recall_mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N)), ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))), width=0.1) +
  labs(x="", y="", title = "", subtitle = "Without Latent Variable_RECALL") +
  theme_classic() + MyTheme

## WITHOUT LV CONDITION: UNCERTAINTY
uncertainties %>% 
  # exclude LV conditions
  filter(!grepl("LV", condition), N == "1000") %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p-sparse", "5p-dense", "10p-sparse", "10p-dense")), y=mean, group = algorithm, colour = algorithm)) +
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
  filter(!grepl("LV", condition), N == "1000") %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p-sparse", "5p-dense", "10p-sparse", "10p-dense")), y=mean, group = algorithm, colour = algorithm)) +
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
  filter(grepl("LV", condition), N == "1000") %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p_LV", "10p_LV")), y=average_precision_mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N)), ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))), width=0.1) +
  labs(x="", y="", title = "", subtitle = "With a Latent Variable_PRECISION") +
  theme_classic() + MyTheme

## WITH LV CONDITION: RECALL
results %>% 
  # exclude LV conditions
  filter(grepl("LV", condition), N == "1000") %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p_LV", "10p_LV")), y=average_recall_mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N)), ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))), width=0.1) +
  labs(x="", y="", title = "", subtitle = "With a Latent Variable_RECALL") +
  theme_classic() + MyTheme

## WITH LV CONDITION: UNCERTAINTY
uncertainties %>% 
  # exclude LV conditions
  filter(grepl("LV", condition), N == "1000") %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p-LV", "10p-LV")), y=mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), width=0.1) +
  labs(x="", y="", title = "", subtitle = "With a Latent Variable_UNCERTAINTY") +
  theme_classic() + MyTheme

## WITH LV CONDITION: SHD
SHDs %>% 
  # exclude LV conditions
  filter(grepl("LV", condition), N == "1000") %>% 
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(condition, levels = c("5p-LV", "10p-LV")), y=mean, group = algorithm, colour = algorithm)) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  geom_line(aes(group = algorithm)) +
  geom_point() +
  # in this case, N = 1000 (hence, sqrt(1000))
  geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), width=0.1) +
  labs(x="", y="", title = "", subtitle = "With a Latent Variable_SHD") +
  theme_classic() + MyTheme





