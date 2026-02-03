library(dplyr)
library(ggplot2)
library(ggpubr)

### To generate this plot, set alpha_list = c("auto") in run-methods/PRSBridge.R and use default percent_list
### For PRS-CS, use phi_list = c(0), name_list = c('_auto'), name_list_prs = c('auto')


### input parameter trait is the same as trait used in run-methods and evaluation

my_plot_auto = function(trait){
  method = method = c('PRS-CS-auto (Small-block)', 'PRS-Bridge (Large-block)', 'PRS-Bridge-auto (Large-block)', 'PRS-CS-auto-1kg (Small-block)', 'PRS-Bridge (Small-block)', 'PRS-Bridge-auto (Small-block)')  
  result = read.table(paste0(trait, '/PRScs_auto_ukbb_validation_result_std.txt'), header = TRUE)
  result = data.frame(method = method[1], R2 = result$R2, std = result$std)
  result_i = 2
  validation_result = read.table(paste0(trait, '/Bridge_large_ukbb_validation_result_std.txt'), header = TRUE)
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/Bridge_large_auto_ukbb_validation_result_std.txt'), header = TRUE) %>% filter(percent == 0.8); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/PRScs_auto_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/Bridge_small_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/Bridge_small_auto_1kg_validation_result_std.txt'), header = TRUE) %>% filter(percent == 0.8); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  result$method = factor(result$method, levels = result$method)
  
  result$ref = c(rep('ukbb', 3), rep('1kg', 3))
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (trait == "resting_heart_rate") {
    trait = "RHR"
  }
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black') + labs(title = paste0(trait))+
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw()+
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 1),
          legend.text = element_text(face = "bold", size = 1), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c('#0ecc41', '#c26802','#c4046f', '#0ecc41', "#c2a482", "#c2709e")) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  #return(plot1)
  return((result[3, 2] - result[1, 2]) / result[3, 2])
}
BMI = my_plot_auto("BMI") + labs(x = "Method", y = expression(bold(R^2))) 
resting_heart_rate = my_plot_auto("resting_heart_rate") + labs(x = "Method", y = expression(bold(R^2))) 
HDL = my_plot_auto("HDL") + labs(x = "Method", y = expression(bold(R^2))) 
LDL = my_plot_auto("LDL") + labs(x = "Method", y = expression(bold(R^2))) 
APOEA = my_plot_auto("APOEA") + labs(x = "Method", y = expression(bold(R^2))) 
APOEB = my_plot_auto("APOEB") + labs(x = "Method", y = expression(bold(R^2))) 
plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("continuous_figure_auto.pdf", plot, width = 8, height = 6)

(my_plot_auto("BMI") + my_plot_auto("HDL") + my_plot_auto("LDL") + my_plot_auto("APOEA") + my_plot_auto("APOEB")) / 6

BC = my_plot_auto("BC") + labs(x = "Method", y = "Transformed AUC") 
CAD = my_plot_auto("CAD") + labs(x = "Method", y = "Transformed AUC") 
Depression = my_plot_auto("Depression") + labs(x = "Method", y = "Transformed AUC") 
RA = my_plot_auto("RA") + labs(x = "Method", y = "Transformed AUC") 
IBD = my_plot_auto("IBD") + labs(x = "Method", y = "Transformed AUC") 
plot = ggarrange(BC+ rremove("xlab"), CAD + rremove("ylab") + rremove("xlab"), Depression + rremove("ylab") + rremove("xlab"), RA+ rremove("xlab"), IBD + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("disease_figure_auto.pdf", plot, width = 8, height = 6)
