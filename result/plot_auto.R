library(dplyr)
library(ggplot2)
library(ggpubr)

###############
my_plot_explore = function(pheno_name){
  result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_auto_ukbb_validation_result_std.txt'), header = TRUE)
  result = rbind(result, read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_auto_1kg_validation_result_std.txt'), header = TRUE))
  result$method = c('1', '2', '3', '4', '5', '6', '7', '8')
  result$method = factor(result$method, levels = result$method)
  
  result$ref = c(rep('ukbb', 4), rep('1kg', 4))
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (pheno_name == "resting_heart_rate") {
    pheno_name = "RHR"
  }
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black')+
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw() +
    labs(x = "Method", y = expression(bold(R^2)), title = paste0(pheno_name)) +
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 1),
          legend.text = element_text(face = "bold", size = 1), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c('#99aecf','#0ecc41',"#c2a482", '#c26802', '#99aecf','#0ecc41',"#c2a482", '#c26802')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
  #return()
}

BMI = my_plot_explore("BMI")
resting_heart_rate = my_plot_explore("resting_heart_rate")
HDL = my_plot_explore("HDL")
LDL = my_plot_explore("LDL")
APOEA = my_plot_explore("APOEA")
APOEB = my_plot_explore("APOEB")
ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")


BC = my_plot_explore("BC")
CAD = my_plot_explore("CAD")
Depression = my_plot_explore("Depression")
RA = my_plot_explore("RA")
IBD = my_plot_explore("IBD")
ggarrange(BC+ rremove("xlab"), CAD + rremove("ylab") + rremove("xlab"), Depression + rremove("ylab") + rremove("xlab"), RA+ rremove("xlab"), IBD + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")


my_plot_auto = function(pheno_name){
  method = c('1', '2', '3', '4', '5', '6')
  result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/PRScs_auto_ukbb_validation_result_std.txt'), header = TRUE)
  result = data.frame(method = method[1], R2 = result$R2, std = result$std)
  result_i = 2
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_ukbb_validation_result_std.txt'), header = TRUE)
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_auto_ukbb_validation_result_std.txt'), header = TRUE) %>% filter(percent == 0.8); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/PRScs_auto_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_auto_1kg_validation_result_std.txt'), header = TRUE) %>% filter(percent == 0.8); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  result$method = factor(result$method, levels = result$method)
  
  result$ref = c(rep('ukbb', 3), rep('1kg', 3))
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (pheno_name == "resting_heart_rate") {
    pheno_name = "RHR"
  }
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black') + labs(title = paste0(pheno_name))+
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
ggsave("~/Documents/PRS_Bridge/figure/continuous_figure_auto.pdf", plot, width = 8, height = 6)

(my_plot_auto("BMI") + my_plot_auto("HDL") + my_plot_auto("LDL") + my_plot_auto("APOEA") + my_plot_auto("APOEB")) / 6

BC = my_plot_auto("BC") + labs(x = "Method", y = "Transformed AUC") 
CAD = my_plot_auto("CAD") + labs(x = "Method", y = "Transformed AUC") 
Depression = my_plot_auto("Depression") + labs(x = "Method", y = "Transformed AUC") 
RA = my_plot_auto("RA") + labs(x = "Method", y = "Transformed AUC") 
IBD = my_plot_auto("IBD") + labs(x = "Method", y = "Transformed AUC") 
plot = ggarrange(BC+ rremove("xlab"), CAD + rremove("ylab") + rremove("xlab"), Depression + rremove("ylab") + rremove("xlab"), RA+ rremove("xlab"), IBD + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("~/Documents/PRS_Bridge/figure/disease_figure_auto.pdf", plot, width = 8, height = 6)



plot_legend = function(pheno_name){
  method = c('PRS-CS-auto (Small-block)', 'PRS-Bridge (Large-block)', 'PRS-Bridge-auto (Large-block)', 'PRS-CS-auto-1kg (Small-block)', 'PRS-Bridge (Small-block)', 'PRS-Bridge-auto (Small-block)')  
  result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_ukbb_validation_result_std.txt'), header = TRUE)
  result = data.frame(method = method[1], R2 = result$R2, std = result$std)
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_auto_ukbb_validation_result_std.txt'), header = TRUE) %>% filter(percent == 0.8)
  result_i = 2
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_auto_1kg_validation_result_std.txt'), header = TRUE) %>% filter(percent == 0.8); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  result$method = factor(result$method, levels = result$method)
  
  result$ref = c(rep('ukbb', 2), rep('1kg', 2))
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (pheno_name == "resting_heart_rate") {
    pheno_name = "RHR"
  }
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black') + labs(title = paste0(pheno_name))+
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw()+
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(face = "bold", size = 12), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c('#c26802','#c4046f',"#c2a482", "#c2709e")) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}
plot_legend("BMI")
