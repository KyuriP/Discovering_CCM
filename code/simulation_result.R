source("code/simulation_code.R")
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)

## Precision & Recall
# compute the average precision and sd
results <- list(res_ccd5psparse, res_fci5psparse, res_cci5psparse, res_ccd10psparse, res_fci10psparse, res_cci10psparse, res_ccd5pdense, res_fci5pdense, res_cci5pdense, res_ccd10pdense, res_fci10pdense, res_cci10pdense, res_ccd5pLVsparse, res_fci5pLVsparse, res_cci5pLVsparse, res_ccd5pLVdense, res_fci5pLVdense, res_cci5pLVdense, res_ccd10pLVsparse, res_fci10pLVsparse, res_cci10pLVsparse, res_ccd10pdense,  res_fci10pLVdense, res_cci10pLVdense) %>% 
  # transpose df
  map(~ sjmisc::rotate_df(.x) %>%
        # add sample size (N) info
        rename_with(~paste0(.x, "N = ", rep(N, each=8)))  %>%
        # think about how to deal with NAs or do I want to define sth. else instead of NAs.
        #na.omit(.x) %>% 
        summarise(across(everything(.), list(mean = ~mean(., na.rm=T), sd = ~sd(., na.rm=T))))) %>% 
  bind_rows() %>% 
  mutate(algorithm = rep(c("ccd", "fci", "cci"), 8),
         condition = rep(c("5p_sparse", "10p_sparse", "5p_dense", "10p_dense", "5p_LVsparse", "5p_LVdense", "10p_LVsparse", "10p_LVdense"), each=3)) %>%
  # brings the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric)) %>% 
  # convert it to a long format
  tidyr::pivot_longer(!c(algorithm, condition), names_to = "metric", values_to = "value") %>% 
  # Add sample size column (N) & clean up the column name 
  mutate(N = stringr::str_extract(metric, "(?<=[N =])\\d+"),
         metric = stringr::str_replace_all(metric, "[0-9.]+|[N =]", "")) 

## Uncertainty
uncertainties <- bind_rows("ccd_5p-sparse" = uncer_ccd5psparse, "fci_5p-sparse" = uncer_fci5psparse, "cci_5p-sparse"=uncer_cci5psparse, "ccd_10p-sparse"=uncer_ccd10psparse, "fci_10p-sparse" = uncer_fci10psparse, "cci_10p-sparse" = uncer_cci10psparse, "ccd_5p-dense"=uncer_ccd5pdense, "fci_5p-dense"=uncer_fci5pdense, "cci_5p-dense"=uncer_cci5pdense, "ccd_10p-dense"=uncer_ccd10pdense, "fci_10p-dense"=uncer_fci10pdense, "cci_10p-dense"=uncer_cci10pdense, "ccd_5p-LVsparse"=uncer_ccd5pLVsparse, "fci_5p-LVsparse"=uncer_fci5pLVsparse, "cci_5p-LVsparse"=uncer_cci5pLVsparse, "ccd_10p-LVsparse"=uncer_ccd10pLVsparse, "fci_10p-LVsparse"=uncer_fci10pLVsparse, "cci_10p-LVsparse"=uncer_cci10pLVsparse,
                           "ccd_5p-LVdense"=uncer_ccd5pLVdense, "fci_5p-LVdense"=uncer_fci5pLVdense, "cci_5p-LVdense"=uncer_cci5pLVdense, "ccd_10p-LVdense"=uncer_ccd10pLVdense, "fci_10p-LVdense"=uncer_fci10pLVdense, "cci_10p-LVdense"=uncer_cci10pLVdense,.id="id") %>% 
  group_by(id) %>% 
  summarise_all(list(mean = mean, sd = sd)) %>%  
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2]) %>% 
  tidyr::pivot_longer(!c(algorithm, condition, id), names_to = "name", values_to = "value") %>% 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  relocate(where(is.character), .before = where(is.numeric))



## SHD
SHDs <- bind_rows("ccd_5p-sparse" = SHD_ccd5psparse, "fci_5p-sparse" = SHD_fci5psparse, "cci_5p-sparse"=SHD_cci5psparse, "ccd_10p-sparse"= SHD_ccd10psparse, "fci_10p-sparse" = SHD_fci10psparse, "cci_10p-sparse" = SHD_cci10psparse, "ccd_5p-dense"= SHD_ccd5pdense, "fci_5p-dense"=SHD_fci5pdense, "cci_5p-dense"=SHD_cci5pdense, "ccd_10p-dense"= SHD_ccd10pdense, "fci_10p-dense"=SHD_fci10pdense, "cci_10p-dense"=SHD_cci10pdense, "ccd_5p-LVsparse"=SHD_ccd5pLVsparse, "fci_5p-LVsparse"=SHD_fci5pLVsparse, "cci_5p-LVsparse"=SHD_cci5pLVsparse, "ccd_10p-LVsparse"=SHD_ccd10pLVsparse, "fci_10p-LVsparse"=SHD_fci10pLVsparse, "cci_10p-LVsparse"=SHD_cci10pLVsparse, 
                  "ccd_5p-LVdense"=SHD_ccd5pLVdense, "fci_5p-LVdense"=SHD_fci5pLVdense, "cci_5p-LVdense"=SHD_cci5pLVdense, "ccd_10p-LVdense"=SHD_ccd10pLVdense, "fci_10p-LVdense"=SHD_fci10pLVdense, "cci_10p-LVdense"=SHD_cci10pLVdense, .id="id") %>% 
  group_by(id) %>% 
  summarise_all(list(means = mean, sds = sd)) %>%  
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2]) %>% 
  tidyr::pivot_longer(!c(algorithm, condition, id), names_to = "name", values_to = "value") %>% 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  relocate(where(is.character), .before = where(is.numeric)) 



## Specify my custom theme
MyTheme <-  theme(plot.title = element_blank(),
                  plot.subtitle = element_text(face = "italic", family = "Palatino", size = 15, hjust=0.5),
                  axis.text=element_text(face = "bold",family = "Palatino", size = 13),
                  axis.text.x = element_text(angle = 45, hjust = 1.2, vjust =1.2),
                  legend.text = element_text(face = "bold", family = "Palatino", size = 13))

## ========================
## precision plots
## ========================
precision_plots <- c("5p_sparse", "5p_dense", "10p_sparse", "10p_dense", "5p_LVsparse", "5p_LVdense", "10p_LVsparse", "10p_LVdense" ) %>% 
  map(~
        results %>% 
        filter(condition == .x & grepl("average_precision", metric)) %>% 
        tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
        ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "5000")), y=average_precision_mean, group = algorithm, colour = algorithm, fill=algorithm)) +
        geom_line(aes(group = algorithm)) +
        geom_point() +
        #geom_errorbar(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N)), ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))), width=0.1) +
        geom_ribbon(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N)), ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))), alpha=0.2, color=NA) +
        scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
        scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
        labs(x="N", y="", title = "", subtitle = .x) +
        scale_y_continuous(limits = c(0.4, 1)) +
        theme_classic() + MyTheme
  )

## ========================
## recall plots
## ========================
recall_plots <-  c("5p_sparse", "5p_dense", "10p_sparse", "10p_dense", "5p_LVsparse", "5p_LVdense", "10p_LVsparse", "10p_LVdense" ) %>% 
  map(~
        results %>% 
        filter(condition == .x & grepl("average_recall", metric)) %>% 
        tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
        ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "5000")), y=average_recall_mean, group = algorithm, colour = algorithm, fill= algorithm)) +
        geom_line(aes(group = algorithm)) +
        geom_point() +
        #geom_errorbar(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N)), ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))), width=0.1) +
        geom_ribbon(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N)), ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))), alpha=0.2, color=NA) +
        scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
        scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
        labs(x="N", y="", title = "", subtitle = .x) +
        scale_y_continuous(limits = c(0.4, 1)) +
        theme_classic() + MyTheme
  )

## ========================
## uncertainty plots
## ========================
## uncertainty (uncertainty rate of fci and cci are exactly the same!)
uncertainty_plots <-  c("5p-sparse", "5p-dense", "10p-sparse", "10p-dense", "5p-LVsparse", "5p-LVdense",  "10p-LVsparse", "10p-LVdense") %>% 
  map(~
        uncertainties %>%
        filter(condition == .x) %>% 
        tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
        ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "5000")), y=mean, group = algorithm, colour = algorithm, fill=algorithm)) +
        geom_line(aes(group = algorithm)) +
        geom_point() +
        #geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), width=0.1) +
        geom_ribbon(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), alpha=0.2, color=NA) +
        scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
        scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
        labs(x="N", y="", title = "", subtitle = .x) +
        scale_y_continuous(limits = c(0.4, 1)) +
        theme_classic() + MyTheme
  )


## ========================
## SHD plots 
## ========================
SHD_plots <- c("5p-sparse", "5p-dense", "10p-sparse", "10p-dense", "5p-LVsparse", "5p-LVdense",  "10p-LVsparse", "10p-LVdense") %>% 
  map(~ 
        SHDs %>%
        filter(condition == .x) %>% 
        tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
        ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "5000")), y=means, group = algorithm, colour = algorithm, fill = algorithm)) +
        geom_line(aes(group = algorithm)) +
        geom_point() + 
        #geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), width=0.1) +
        geom_ribbon(aes(ymin=means-qnorm(0.975)*sds/sqrt(as.numeric(N)), ymax=means+qnorm(0.975)*sds/sqrt(as.numeric(N))), alpha=0.2, color=NA) +
        scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
        scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
        labs(x="N", y="", title = "", subtitle = .x) +
        #scale_y_continuous(limits = ~c((min(means)-5), (max(means)+5))) +
        theme_classic() + MyTheme
  )




# combine plots
# precision plot
ggarrange(plotlist = precision_plots,
                            ncol = 4, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(top = text_grob("Precision", face = "bold", size = 18, family = "Palatino"))
# save the figure
ggsave("figure/precisionplots.pdf", width = 10, height=6, units = "in")

# recall plot
ggarrange(plotlist = recall_plots,
                         ncol = 4, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(top = text_grob("Recall", face = "bold", size = 18, family = "Palatino"))
# save the figure
ggsave("figure/recallplots.pdf", width = 10, height=6, units = "in")

# uncertainty plot
ggarrange(plotlist = uncertainty_plots,
                        ncol = 4, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(top = text_grob("Uncertainty", face = "bold", size = 18, family = "Palatino"))
# save the figure
ggsave("figure/uncertaintyplots.pdf", width = 10, height=6, units = "in")

# shd plot
ggarrange(plotlist = SHD_plots,
                      ncol = 4, nrow = 2, common.legend = TRUE, legend = "bottom") %>% 
  annotate_figure(top = text_grob("SHD", face = "bold", size = 15, family = "Palatino"))
# save the figure
ggsave("figure/shdplots.pdf", width = 10, height=6, units = "in")


## ========================
## TIME plots 
## ========================
times %>%
  ggplot(aes(x=factor(condition, levels= c("5psparse", "5pdense", "10psparse", "10pdense", "5pLVsparse","5pLVdense", "10pLVsparse","10pLVdense")), y = log(time), col= factor(algorithm))) +
  geom_boxplot(position = "dodge",   outlier.size = 0.8, outlier.alpha = 0.2) + theme_classic() +
  # scale_x_discrete(name ="Condition",
  #                  labels=c("", "5p-sparse", "", "","5p-dense","","", "10p-sparse","","","10p-dense","","","5p-LV","","","10p-LV","")) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  labs(y = " log(ms)", x = "", title = "Algorithm Running Time", subtitle = "Time in milliseconds (ms)") +
  theme(axis.text = element_text(face = "bold", family =  "Palatino", margin = margin(t = 13), size=15),
        legend.position="bottom",
        plot.subtitle = element_text(face = "italic", family = "Palatino", size=16),
        plot.title = element_text(family = "Palatino", size=18, face="bold"),
        axis.title = element_text(family = "Palatino", size=15),
        legend.text = element_text(family =  "Palatino", size=13, face="bold"))

# save the figure
ggsave("figure/timeplots.pdf", width = 11, height=4, units = "in")


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





