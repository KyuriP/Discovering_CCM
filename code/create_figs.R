
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
  mutate(algorithm = rep(c("CCD", "FCI", "CCI"), 8),
         condition = rep(c("5p_sparse", "10p_sparse", "5p_dense", "10p_dense", "5p_LVsparse", "5p_LVdense", "10p_LVsparse", "10p_LVdense"), each=3),
         latentvar = ifelse(stringr::str_detect(condition, "LV") == TRUE, "with LV", "without LV"),
         netsize = stringr::str_split(condition, "_", simplify = T)[,1],
         densities = stringr::str_remove(stringr::str_split(condition, "_", simplify = T)[,2], "LV")
           ) %>%
  # brings the algorithm and condition names first
  relocate(where(is.character), .before = where(is.numeric)) %>% 
  # convert it to a long format
  tidyr::pivot_longer(!c(algorithm, condition, latentvar, netsize, densities), names_to = "metric", values_to = "value") %>% 
  # Add sample size column (N) & clean up the column name 
  mutate(N = stringr::str_extract(metric, "(?<=[N =])\\d+"),
         metric = stringr::str_replace_all(metric, "[0-9.]+|[N =]", "")) 

## Uncertainty
uncertainties <- bind_rows("ccd_5p-sparse" = uncer_ccd5psparse, "fci_5p-sparse" = uncer_fci5psparse, "cci_5p-sparse"=uncer_cci5psparse, "ccd_10p-sparse"=uncer_ccd10psparse, "fci_10p-sparse" = uncer_fci10psparse, "cci_10p-sparse" = uncer_cci10psparse, "ccd_5p-dense"=uncer_ccd5pdense, "fci_5p-dense"=uncer_fci5pdense, "cci_5p-dense"=uncer_cci5pdense, "ccd_10p-dense"=uncer_ccd10pdense, "fci_10p-dense"=uncer_fci10pdense, "cci_10p-dense"=uncer_cci10pdense, "ccd_5p-LVsparse"=uncer_ccd5pLVsparse, "fci_5p-LVsparse"=uncer_fci5pLVsparse, "cci_5p-LVsparse"=uncer_cci5pLVsparse, "ccd_10p-LVsparse"=uncer_ccd10pLVsparse, "fci_10p-LVsparse"=uncer_fci10pLVsparse, "cci_10p-LVsparse"=uncer_cci10pLVsparse,
                           "ccd_5p-LVdense"=uncer_ccd5pLVdense, "fci_5p-LVdense"=uncer_fci5pLVdense, "cci_5p-LVdense"=uncer_cci5pLVdense, "ccd_10p-LVdense"=uncer_ccd10pLVdense, "fci_10p-LVdense"=uncer_fci10pLVdense, "cci_10p-LVdense"=uncer_cci10pLVdense,.id="id") %>% 
  group_by(id) %>% 
  summarise_all(list(means = mean, sds = sd)) %>%  
  mutate(algorithm = toupper(stringr::str_split(id, "_", simplify = T)[,1]),
         condition = stringr::str_split(id, "_", simplify = T)[,2],
         latentvar = ifelse(stringr::str_detect(condition, "LV") == TRUE, "with LV", "without LV"),
         netsize = stringr::str_split(condition, "-", simplify = T)[,1],
         densities = stringr::str_remove(stringr::str_split(condition, "-", simplify = T)[,2], "LV")
         ) %>% 
  tidyr::pivot_longer(!c(algorithm, condition, id, netsize, latentvar, densities), names_to = "name", values_to = "value") %>% 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  relocate(where(is.character), .before = where(is.numeric))



## SHD
SHDs <- bind_rows("ccd_5p-sparse" = SHD_ccd5psparse, "fci_5p-sparse" = SHD_fci5psparse, "cci_5p-sparse"=SHD_cci5psparse, "ccd_10p-sparse"= SHD_ccd10psparse, "fci_10p-sparse" = SHD_fci10psparse, "cci_10p-sparse" = SHD_cci10psparse, "ccd_5p-dense"= SHD_ccd5pdense, "fci_5p-dense"=SHD_fci5pdense, "cci_5p-dense"=SHD_cci5pdense, "ccd_10p-dense"= SHD_ccd10pdense, "fci_10p-dense"=SHD_fci10pdense, "cci_10p-dense"=SHD_cci10pdense, "ccd_5p-LVsparse"=SHD_ccd5pLVsparse, "fci_5p-LVsparse"=SHD_fci5pLVsparse, "cci_5p-LVsparse"=SHD_cci5pLVsparse, "ccd_10p-LVsparse"=SHD_ccd10pLVsparse, "fci_10p-LVsparse"=SHD_fci10pLVsparse, "cci_10p-LVsparse"=SHD_cci10pLVsparse, 
                  "ccd_5p-LVdense"=SHD_ccd5pLVdense, "fci_5p-LVdense"=SHD_fci5pLVdense, "cci_5p-LVdense"=SHD_cci5pLVdense, "ccd_10p-LVdense"=SHD_ccd10pLVdense, "fci_10p-LVdense"=SHD_fci10pLVdense, "cci_10p-LVdense"=SHD_cci10pLVdense, .id="id") %>% 
  group_by(id) %>% 
  summarise_all(list(means = mean, sds = sd)) %>%  
  mutate(algorithm = stringr::str_split(id, "_", simplify = T)[,1],
         condition = stringr::str_split(id, "_", simplify = T)[,2],
         latentvar = ifelse(stringr::str_detect(condition, "LV") == TRUE, "with LV", "without LV"),
         netsize = stringr::str_split(condition, "-", simplify = T)[,1],
         densities = stringr::str_remove(stringr::str_split(condition, "-", simplify = T)[,2], "LV")
         ) %>% 
  tidyr::pivot_longer(!c(algorithm, condition, id, latentvar, netsize, densities), names_to = "name", values_to = "value") %>% 
  mutate(N = stringr::str_extract(stringr::str_split(name, "_", simplify = T)[,1], "(\\d)+"),
         statistics = stringr::str_split(name, "_", simplify = T)[,2]) %>% 
  dplyr::select(-id, -name) %>%  relocate(where(is.character), .before = where(is.numeric)) 





library(ggh4x)


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
                  panel.border = element_rect(color = "#DCDCDC", fill = NA)
                  )


## ========================
## SHD plots
## ========================
SHDs %>%
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("100", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000", "20000")), y=means, group = algorithm, colour = algorithm, fill = algorithm)) +
  geom_line(aes(group = algorithm)) +
  geom_point(size=1) + 
  #geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), width=0.1) +
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by 3)
  geom_ribbon(aes(ymin=means-qnorm(0.975)*sds/sqrt(as.numeric(N))*3, ymax=means+qnorm(0.975)*sds/sqrt(as.numeric(N))*3), alpha=0.2, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x="N", y="", title = "") +
  theme_minimal() +
  MyTheme + 
  facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  scales = "free_y", switch="y") +
  ggtitle("SHD")

## save shd plot
ggsave(filename = "results/SHD.pdf", width = 25, height = 13, dpi = 300, units = "cm")


## ========================
## precision & recall plots
## ========================
precision_plot <- results %>% 
  filter(grepl("average_precision", metric)) %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("100", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000", "20000")), y=average_precision_mean, group = algorithm, colour = algorithm, fill=algorithm)) +
  geom_line(aes(group = algorithm)) +
  geom_point(size=1) +
  #geom_errorbar(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N)), ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))), width=0.1) +
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by 2)
  geom_ribbon(aes(ymin=average_precision_mean-qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))*2, ymax=average_precision_mean+qnorm(0.975)*average_precision_sd/sqrt(as.numeric(N))*2), alpha=0.2, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  theme_minimal() +
  MyTheme + 
  facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  switch="y") +
  labs(title = "Precision", x = "N", y = "")

# ## save precision plot
# ggsave(filename = "results/precision.pdf", width = 25, height = 13, dpi = 300, units = "cm")


recall_plot <- results %>% 
  filter(grepl("average_recall", metric)) %>% 
  tidyr::pivot_wider(names_from = metric, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("100", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000", "20000")), y=average_recall_mean, group = algorithm, colour = algorithm, fill= algorithm)) +
  geom_line(aes(group = algorithm)) +
  geom_point(size=1) +
  #geom_errorbar(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N)), ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))), width=0.1) +
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by 2)
  geom_ribbon(aes(ymin=average_recall_mean-qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))*2, ymax=average_recall_mean+qnorm(0.975)*average_recall_sd/sqrt(as.numeric(N))*2), alpha=0.2, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  theme_minimal() +
  MyTheme + 
  facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  switch="y") +
  labs(title = "Recall", x = "N", y = "")
  
ggpubr::ggarrange(precision_plot, recall_plot, nrow=2, common.legend = TRUE, legend = "bottom")

## save precision&recall plot
ggsave(filename = "results/prec-recall.pdf", width = 25, height = 26, dpi = 300, units = "cm")


## ========================
## Uncertainty plots
## ========================
uncertainties %>%
  tidyr::pivot_wider(names_from = statistics, values_from=value) %>% 
  ggplot(aes(x= factor(N, levels = c("100", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000", "20000")), y=means, group = algorithm, colour = algorithm, fill = algorithm)) +
  geom_line(aes(group = algorithm)) +
  geom_point(size=1) + 
  #geom_errorbar(aes(ymin=mean-qnorm(0.975)*sd/sqrt(as.numeric(N)), ymax=mean+qnorm(0.975)*sd/sqrt(as.numeric(N))), width=0.1) +
  # exaggerate the intervals a bit to ensure they are visible in the plot (times by )
  geom_ribbon(aes(ymin=means-qnorm(0.975)*sds/sqrt(as.numeric(N))*2, ymax=means+qnorm(0.975)*sds/sqrt(as.numeric(N))*2), alpha=0.2, color=NA) +
  scale_colour_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00"), name= "") +
  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x="N", y="", title = "") +
  theme_minimal() +
  MyTheme + 
  facet_nested(factor(netsize, levels = c("5p", "10p")) ~ factor(latentvar, levels = c("without LV", "with LV")) + factor(densities, levels=c("sparse", "dense")),  scales = "free_y", switch="y") +
  ggtitle("Uncertainty")

## save shd plot
ggsave(filename = "results/uncertainty.pdf", width = 25, height = 13, dpi = 300, units = "cm")

