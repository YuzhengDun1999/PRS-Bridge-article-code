library(dplyr)
library(ggplot2)
library(ggpubr)

##################### plot lassosum ##################
method = c('Lassosum (Small-block)', 'LDpred2 (Banded)', 'LDpred2 (Large-block)', 'PRS-CS (Small-block)', 'PRS-Bridge (Small-block)', 'PRS-Bridge (Large-block)', 'Lassosum2_1kg', 'ldpred2_banded_1kg', 'PRScs_1kg', 'PRSBridge_small_block_1kg', 'PRSBridge_large_block_1kg')  

my_plot_new = function(pheno_name){
  method = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11')
  result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Lassosum_N5000_validation_result_std.txt'), header = TRUE)
  result = data.frame(method = method[1], R2 = result$R2, std = result$std)
  result_i = 1
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/ldpred2_N5000_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/ldpred2_Nblk_provided_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/PRScs_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Lassosum_N1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/ldpred2_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/PRScs_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  
  result$method = factor(result$method, levels = method)
  
  result$ref = c(rep('ukbb', 6), rep('1kg', 5))
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (pheno_name == "resting_heart_rate") {
    pheno_name = "RHR"
  }
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black')+
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw() +
    labs(x = "Method", title = paste0(pheno_name)) +
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 1),
          legend.text = element_text(face = "bold", size = 1), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c("#e3c278", "#99aecf", "#0250c9", "#0ecc41", "#c2a482", '#c26802', '#e3c278', '#99aecf','#0ecc41',"#c2a482", '#c26802')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return((result[6, 2] - result[4, 2]) / result[6, 2])
  #return(plot1)
}

BMI = my_plot_new("BMI") + labs(y = expression(bold(R^2)))
resting_heart_rate = my_plot_new("resting_heart_rate") + labs(y = expression(bold(R^2)))
HDL = my_plot_new("HDL") + labs(y = expression(bold(R^2)))
LDL = my_plot_new("LDL") + labs(y = expression(bold(R^2)))
APOEA = my_plot_new("APOEA") + labs(y = expression(bold(R^2)))
APOEB = my_plot_new("APOEB") + labs(y = expression(bold(R^2)))
(my_plot_new("BMI") + my_plot_new("resting_heart_rate") + my_plot_new("HDL") + 
    my_plot_new("LDL") + my_plot_new("APOEA") + my_plot_new("APOEB")) / 6

plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("~/Documents/PRS_Bridge/figure/continuous_figure.pdf", plot, width = 8, height = 6)

BC = my_plot_new("BC") + labs(y = "Transformed AUC")
CAD = my_plot_new("CAD") + labs(y = "Transformed AUC")
Depression = my_plot_new("Depression") + labs(y = "Transformed AUC")
RA = my_plot_new("RA") + labs(y = "Transformed AUC")
IBD = my_plot_new("IBD") + labs(y = "Transformed AUC")
plot = ggarrange(BC+ rremove("xlab"), CAD + rremove("ylab") + rremove("xlab"), Depression + rremove("ylab") + rremove("xlab"), RA+ rremove("xlab"), IBD + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
(my_plot_new("BC") + my_plot_new("CAD") + my_plot_new("Depression") + 
    my_plot_new("RA") + my_plot_new("IBD") + my_plot_new("APOEB")) / 5

ggsave("~/Documents/PRS_Bridge/figure/disease_figure.pdf", plot, width = 8, height = 6)



##################### exclude lassosum ##################
method = c('LDpred2 (Banded)', 'LDpred2 (Large-block)', 'PRS-CS (Small-block)', 'PRS-Bridge (Small-block)', 'PRS-Bridge (Large-block)', 'ldpred2_banded_1kg', 'PRScs_1kg', 'PRSBridge_small_block_1kg', 'PRSBridge_large_block_1kg')  

my_plot_new = function(pheno_name){
  method = c('1', '2', '3', '4', '5', '6', '7', '8', '9')
  result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/ldpred2_N5000_validation_result_std.txt'), header = TRUE)
  result = data.frame(method = method[1], R2 = result$R2, std = result$std)
  result_i = 1
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/ldpred2_Nblk_provided_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/PRScs_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/ldpred2_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/PRScs_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  
  result$method = factor(result$method, levels = method)
  
  result$ref = c(rep('ukbb', 5), rep('1kg', 4))
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (pheno_name == "resting_heart_rate") {
    pheno_name = "RHR"
  }
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black')+
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw() +
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
  #return((result[6, 2] - result[4, 2]) / result[6, 2])
  return(plot1)
}

BMI = my_plot_new("BMI") + labs(y = expression(bold(R^2)))
resting_heart_rate = my_plot_new("resting_heart_rate") + labs(y = expression(bold(R^2)))
HDL = my_plot_new("HDL") + labs(y = expression(bold(R^2)))
LDL = my_plot_new("LDL") + labs(y = expression(bold(R^2)))
APOEA = my_plot_new("APOEA") + labs(y = expression(bold(R^2)))
APOEB = my_plot_new("APOEB") + labs(y = expression(bold(R^2)))
(my_plot_new("BMI") + my_plot_new("resting_heart_rate") + my_plot_new("HDL") + 
    my_plot_new("LDL") + my_plot_new("APOEA") + my_plot_new("APOEB")) / 6

plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("~/Documents/PRS_Bridge/figure/continuous_figure.pdf", plot, width = 8, height = 6)

BC = my_plot_new("BC") + labs(y = "Transformed AUC")
CAD = my_plot_new("CAD") + labs(y = "Transformed AUC")
Depression = my_plot_new("Depression") + labs(y = "Transformed AUC")
RA = my_plot_new("RA") + labs(y = "Transformed AUC")
IBD = my_plot_new("IBD") + labs(y = "Transformed AUC")
plot = ggarrange(BC+ rremove("xlab"), CAD + rremove("ylab") + rremove("xlab"), Depression + rremove("ylab") + rremove("xlab"), RA+ rremove("xlab"), IBD + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
(my_plot_new("BC") + my_plot_new("CAD") + my_plot_new("Depression") + 
    my_plot_new("RA") + my_plot_new("IBD") + my_plot_new("APOEB")) / 5

ggsave("~/Documents/PRS_Bridge/figure/disease_figure.pdf", plot, width = 8, height = 6)


plot_legend = function(pheno_name){
  method = c('Lassosum2 (Banded)', 'Lassosum2 (Large-block)', 'LDpred2 (Banded)', 'LDpred2 (Large-block)', 'PRS-CS (Small-block)', 'PRS-Bridge (Small-block)', 'PRS-Bridge (Large-block)', 'Lassosum2_1kg', 'ldpred2_banded_1kg', 'PRScs_1kg', 'PRSBridge_small_block_1kg', 'PRSBridge_large_block_1kg')  
  result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Lassosum2_N5000_validation_result_std.txt'), header = TRUE)
  result = data.frame(method = method[1], R2 = result$R2, std = result$std)
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Lassosum2_Nblk_provided_validation_result_std.txt'), header = TRUE)
  result_i = 2
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/ldpred2_N5000_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/ldpred2_Nblk_provided_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/PRScs_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Lassosum2_N1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/ldpred2_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/PRScs_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_small_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0('~/Documents/PRS_Bridge/result_sd/', pheno_name, '/Bridge_large_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  
  result$ref = c(rep('UKBB', 7), rep('1000G',5))
  result$ref = factor(result$ref, levels = c('UKBB', '1000G'))
  result$method = factor(result$method, levels = method)
  
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black') + theme_bw() +
    labs(x = "Method", y = "R2", title = paste0(pheno_name)) +
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(face = "bold", size = 12), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5)) + 
    scale_fill_manual(values=c("#e3c278", '#E69F00', "#99aecf", "#0250c9", "#0ecc41", "#c2a482", '#c26802', '#e3c278', '#99aecf','#0ecc41',"#c2a482", '#c26802')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  
  return(plot1)
}
plot_legend("BMI")


#################### plot PRS-CS ###########
###### summarize result
result = data.frame(matrix(ncol = 5))
names(result) = c("trait", "a", "ref", "R2", "std"); result_i = 1
for (trait in c("BMI", "HDL", "LDL", "resting_heart_rate", "APOEA", "APOEB")) {
  for (ref in c("ukbb", "1kg")) {
    PRScs_a0.5 = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_a0.5_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, 0.5, ref, PRScs_a0.5$R2, PRScs_a0.5$std); result_i = result_i + 1
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, 1, ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    PRScs_a1.5 = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_a1.5_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, 1.5, ref, PRScs_a1.5$R2, PRScs_a1.5$std); result_i = result_i + 1
  }
}
readr::write_tsv(result, "/users/ydun/PRScs_a.txt")

my_plot_new = function(pheno_name){
  method = c('a=0.5', 'a=1', 'a=1.5', 'a=0.5 1kg', 'a=1 1kg', 'a=1.5 1kg')  
  result = bigreadr::fread2("~/Documents/PRS_Bridge/result_sd/PRScs_a.txt")
  result = result %>% filter(trait == pheno_name)
  result = result %>% mutate(method = ifelse(ref == "ukbb", paste0("a=", a), paste0("a=", a, " 1kg")))
  result$method = factor(result$method, levels = method)
  if (pheno_name == "resting_heart_rate") {
    pheno_name = "RHR"
  }
  result$ref = c(rep('ukbb', 3), rep('1kg', 3))
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black')+
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw() +
    labs(x = "Method", title = paste0(pheno_name)) +
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(face = "bold", size = 14), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c("#cc9d0e", "#0ecc41", "#0eb3cc", "#cc9d0e", "#0ecc41", "#0eb3cc")) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

BMI = my_plot_new("BMI") + labs(y = expression(bold(R^2)))
resting_heart_rate = my_plot_new("resting_heart_rate") + labs(y = expression(bold(R^2)))
HDL = my_plot_new("HDL") + labs(y = expression(bold(R^2)))
LDL = my_plot_new("LDL") + labs(y = expression(bold(R^2))) 
APOEA = my_plot_new("APOEA") + labs(y = expression(bold(R^2)))
APOEB = my_plot_new("APOEB") + labs(y = expression(bold(R^2)))
plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("~/Documents/PRS_Bridge/figure/PRScs_a.pdf", plot, width = 8, height = 6)
