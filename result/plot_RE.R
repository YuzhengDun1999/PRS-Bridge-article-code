library(dplyr)
library(ggplot2)
library(ggpubr)

##################### exclude lassosum ##################
method = c('LDpred2 (Banded)', 'LDpred2 (Large-block)', 'PRS-CS (Small-block)', 'PRS-Bridge (Small-block)', 'PRS-Bridge (Large-block)', 'ldpred2_banded_1kg', 'PRScs_1kg', 'PRSBridge_small_block_1kg', 'PRSBridge_large_block_1kg')  

my_plot_new = function(pheno_name){
  method = c('1', '2', '3', '4', '5', '6', '7', '8', '9')
  result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/RR.txt'), header = TRUE)
  #result$RR = result$RR - 1
  result = result %>% filter(trait == pheno_name)
  result$method = method
  result$method = factor(result$method, levels = method)
  
  result$ref = c(rep('ukbb', 5), rep('1kg', 4))
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (pheno_name == "resting_heart_rate") {
    pheno_name = "RHR"
  }
  plot1 = ggplot(data = result, mapping = aes(x = method, y = RR, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black')+
    geom_errorbar(aes(ymin = RR - 1.96 * std, ymax = RR + 1.96 * std), width = 0.15) + theme_bw() +  
    geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1) +
    labs(x = "Method", title = paste0(pheno_name)) +
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 1),
          legend.text = element_text(face = "bold", size = 1), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c("#99aecf", "#0250c9", "#0ecc41", "#c2a482", '#c26802', '#99aecf','#0ecc41',"#c2a482", '#c26802')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

BMI = my_plot_new("BMI") + labs(y = expression(bold("Ratio of "~R^2)))
resting_heart_rate = my_plot_new("resting_heart_rate") + labs(y = expression(bold("Ratio of "~R^2)))
HDL = my_plot_new("HDL") + labs(y = expression(bold("Ratio of "~R^2)))
LDL = my_plot_new("LDL") + labs(y = expression(bold("Ratio of "~R^2)))
APOEA = my_plot_new("APOEA") + labs(y = expression(bold("Ratio of "~R^2)))
APOEB = my_plot_new("APOEB") + labs(y = expression(bold("Ratio of "~R^2)))
plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("~/Documents/PRS_Bridge/figure/continuous_figure_RE.pdf", plot, width = 8, height = 6)

BC = my_plot_new("BC") + labs(y = "Ratio of Transformed AUC")
CAD = my_plot_new("CAD") + labs(y = "Ratio of Transformed AUC")
Depression = my_plot_new("Depression") + labs(y = "Ratio of Transformed AUC")
RA = my_plot_new("RA") + labs(y = "Ratio of Transformed AUC")
IBD = my_plot_new("IBD") + labs(y = "Ratio of Transformed AUC")
plot = ggarrange(BC+ rremove("xlab"), CAD + rremove("ylab") + rremove("xlab"), Depression + rremove("ylab") + rremove("xlab"), RA+ rremove("xlab"), IBD + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")

ggsave("~/Documents/PRS_Bridge/figure/disease_figure_RE.pdf", plot, width = 8, height = 6)
