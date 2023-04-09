library(ggh4x)


# specify the figure theme
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
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000")), y=means, group = algorithm, colour = algorithm, fill = algorithm)) +
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
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000")), y=average_precision_mean, group = algorithm, colour = algorithm, fill=algorithm)) +
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
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000")), y=average_recall_mean, group = algorithm, colour = algorithm, fill= algorithm)) +
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
  ggplot(aes(x= factor(N, levels = c("50", "150", "500", "1000", "1500", "2000", "2500", "3000", "4000", "5000", "10000")), y=means, group = algorithm, colour = algorithm, fill = algorithm)) +
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
