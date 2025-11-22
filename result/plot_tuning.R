###### plot tuning setting
library(dplyr)
library(ggplot2)
library(ggpubr)

plot_tuning = function(trait){
  result = read.table(paste0("~/Documents/PRS_Bridge/result_sd/", trait, "/Bridge_small_ukbb_all_result.txt"), header = TRUE)
  result = result %>% filter(percent != 0 )
  result_sd = result %>% group_by(alpha, percent) %>% summarise( mean_R2 = mean(tuning_R2), se_R2   = sd(tuning_R2) / sqrt(n()), .groups = "drop") %>% as.data.frame()
  
  if (trait == "resting_heart_rate") {
    trait = "RHR"
  }
  plot1 = ggplot(result_sd, aes(x = factor(percent), y = mean_R2, group = factor(alpha), color = factor(alpha))) +
    geom_point(position = position_dodge(width = 0.3), size = 2) +
    geom_line(position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = mean_R2 - 1.96 * se_R2, ymax = mean_R2 + 1.96 * se_R2),
                  width = 0.2, position = position_dodge(width = 0.3)) + theme_bw() +
    labs(x = "Percent", y = expression(bold(R^2)), color = expression(bold(alpha)), title = trait) +
    theme(axis.text.x=element_text(face="bold", size = 14), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 20),
          legend.text = element_text(face = "bold", size = 14), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = c("0.125" = "#1b9e77", "0.25"  = "#d95f02", "0.5" = "#7570b3"))
  return(plot1)
}
BMI = plot_tuning("BMI") + theme(axis.text.x=element_blank())
resting_heart_rate = plot_tuning("resting_heart_rate") + theme(axis.text.x=element_blank())
HDL = plot_tuning("HDL") + theme(axis.text.x=element_blank())
LDL = plot_tuning("LDL")
APOEA = plot_tuning("APOEA")
APOEB = plot_tuning("APOEB")

plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL, APOEA + rremove("ylab"), APOEB+ rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="right")
ggsave("~/Documents/PRS_Bridge/prs-bridge-manuscript/JASA_application/tuning_figure.pdf", plot, width = 8, height = 6)

BC = plot_tuning("BC") + theme(axis.text.x=element_blank())
CAD = plot_tuning("CAD") + theme(axis.text.x=element_blank())
Depression = plot_tuning("Depression") + theme(axis.text.x=element_blank())
RA = plot_tuning("RA")
IBD = plot_tuning("IBD")

plot = ggarrange(BC+ rremove("xlab"), CAD + rremove("ylab") + rremove("xlab"), Depression + rremove("ylab") + rremove("xlab"), RA+ rremove("xlab"), IBD + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="right")
ggsave("~/Documents/PRS_Bridge/prs-bridge-manuscript/JASA_application/tuning_disease_figure.pdf", plot, width = 8, height = 6)





###### plot tuning setting
library(dplyr)
library(ggplot2)
library(ggpubr)

plot_tuning = function(trait){
  result = read.table(paste0("~/Documents/PRS_Bridge/result_sd/", trait, "/Bridge_large_ukbb_all_result.txt"), header = TRUE)
  result = result %>% filter(percent != 0)
  result_sd = result %>% group_by(alpha, percent) %>% summarise( mean_R2 = mean(tuning_R2), se_R2 = sd(tuning_R2), .groups = "drop") %>% as.data.frame()
  best_setting = result_sd %>% filter(mean_R2 == max(mean_R2))
  best_alpha = best_setting$alpha
  result_sd = result_sd %>% filter(alpha == best_alpha)
  result_sd$ref = "ukbb"
  
  result = read.table(paste0("~/Documents/PRS_Bridge/result_sd/", trait, "/Bridge_small_1kg_all_result.txt"), header = TRUE)
  result = result %>% filter(percent != 0)
  result_sd1 = result %>% group_by(alpha, percent) %>% summarise( mean_R2 = mean(tuning_R2), se_R2   = sd(tuning_R2), .groups = "drop") %>% as.data.frame()
  best_setting = result_sd1 %>% filter(mean_R2 == max(mean_R2))
  best_alpha = best_setting$alpha
  result_sd1 = result_sd1 %>% filter(alpha == best_alpha)
  result_sd1$ref = "1kg"
  result_sd = rbind(result_sd, result_sd1)
  
  if (trait == "resting_heart_rate") {
    trait = "RHR"
  }
  plot1 = ggplot(result_sd, aes(x = factor(percent), y = mean_R2, group = factor(ref), color = factor(ref))) +
    geom_point(position = position_dodge(width = 0.3), size = 1) +
    geom_line(position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = mean_R2 - 1.96 * se_R2, ymax = mean_R2 + 1.96 * se_R2),
                  width = 0.3, position = position_dodge(width = 0.3)) + theme_bw() +
    labs(x = "Percent", y = "Transformed AUC", color = "Reference", title = trait) +
    theme(axis.text.x=element_text(face="bold", size = 14), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_blank(),
          legend.text = element_text(face = "bold", size = 14), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = c("ukbb" = "#1b9e77", "1kg"  = "#d95f02"))
  return(plot1)
}
BMI = plot_tuning("BMI") + theme(axis.text.x=element_blank())
resting_heart_rate = plot_tuning("resting_heart_rate") + theme(axis.text.x=element_blank())
HDL = plot_tuning("HDL") + theme(axis.text.x=element_blank())
LDL = plot_tuning("LDL")
APOEA = plot_tuning("APOEA")
APOEB = plot_tuning("APOEB")

plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL, APOEA + rremove("ylab"), APOEB+ rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="right")
ggsave("~/Documents/PRS_Bridge/prs-bridge-manuscript/JASA_application/tuning_figure.pdf", plot, width = 8, height = 6)

BC = plot_tuning("BC") + theme(axis.text.x=element_blank())
CAD = plot_tuning("CAD") + theme(axis.text.x=element_blank())
Depression = plot_tuning("Depression") + theme(axis.text.x=element_blank())
RA = plot_tuning("RA")
IBD = plot_tuning("IBD")

plot = ggarrange(BC+ rremove("xlab"), CAD + rremove("ylab") + rremove("xlab"), Depression + rremove("ylab") + rremove("xlab"), RA+ rremove("xlab"), IBD + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="right")
ggsave("~/Documents/PRS_Bridge/prs-bridge-manuscript/JASA_application/tuning_disease_figure.pdf", plot, width = 8, height = 6)

