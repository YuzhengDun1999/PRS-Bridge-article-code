######## summarize result on JHPCE #######
result = data.frame(matrix(ncol = 5))
names(result) = c("trait", "method", "ref", "R2", "std")
result_i = 1
for (trait in c("BMI", "HDL", "LDL", "resting_heart_rate", "APOEA", "APOEB")) {
  for (ref in c("ukbb", "1kg")) {
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs", ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    PRScs_auto = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_auto_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_auto", ref, PRScs_auto$R2, PRScs_auto$std); result_i = result_i + 1
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_diff_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_proj_null", ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    PRScs_auto = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_diff_auto_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_proj_null_auto", ref, PRScs_auto$R2, PRScs_auto$std); result_i = result_i + 1
  }
}
readr::write_tsv(result, "/users/ydun/PRScs_proj_result.txt")


result = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_result.txt", header = TRUE)

library(dplyr)
library(ggplot2)
library(ggpubr)

############### PRScs-proj-null
my_plot_explore = function(pheno_name){
  result = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_result.txt", header = TRUE)
  result_threshold = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_threshold.txt", header = TRUE)
  names(result_threshold)[2] = "method"
  result_a = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_a.txt", header = TRUE)
  result_a = result_a %>% filter(a == 1.0) %>% select(trait, a, ref, R2, std)
  names(result_a)[2] = "method"
  result = rbind(result, result_threshold, result_a)
  result = result %>% filter(trait == pheno_name) %>% filter(method != "PRScs_auto" & method != "PRScs_proj_null_auto")
  result$method = c('PRScs', 'PRScs_auto', 'PRScs_proj_null', 'PRScs_proj_null_auto', '5', '6', '7', '8')
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
          legend.text = element_text(face = "bold", size = 12), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c('#99aecf','#004fcc',"#c2a482", '#c26802', '#99aecf','#004fcc',"#c2a482", '#c26802')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

BMI = my_plot_explore("BMI")
resting_heart_rate = my_plot_explore("resting_heart_rate")
HDL = my_plot_explore("HDL")
LDL = my_plot_explore("LDL")
APOEA = my_plot_explore("APOEA")
APOEB = my_plot_explore("APOEB")
ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")




########################### PRScs-proj-diag #########
######## summarize result on JHPCE #######
result = data.frame(matrix(ncol = 5))
names(result) = c("trait", "method", "ref", "R2", "std")
result_i = 1
for (trait in c("BMI", "HDL", "LDL", "resting_heart_rate", "APOEA", "APOEB")) {
  for (ref in c("ukbb", "1kg")) {
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs", ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    PRScs_auto = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_auto_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_auto", ref, PRScs_auto$R2, PRScs_auto$std); result_i = result_i + 1
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_diag_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_proj_diag", ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    PRScs_auto = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_diag_auto_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_proj_diag_auto", ref, PRScs_auto$R2, PRScs_auto$std); result_i = result_i + 1
  }
}
result[2, 4] = 0.1114434; result[2, 5] = 0.003105511
readr::write_tsv(result, "/users/ydun/PRScs_proj_diag_result.txt")


result = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_diag_result.txt", header = TRUE)

library(dplyr)
library(ggplot2)
library(ggpubr)

############### PRScs-proj-diag
my_plot_explore = function(pheno_name){
  result = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_diag_result.txt", header = TRUE)
  result = result %>% filter(trait == pheno_name) %>% filter(method != 'PRScs_auto') %>% filter(method != 'PRScs_proj_diag_auto')
  result$method = c('PRScs', 'PRScs_proj_diag', '3', '4')
  result$method = factor(result$method, levels = result$method)
  
  result$ref = c(rep('ukbb', 2), rep('1kg', 2))
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
          legend.text = element_text(face = "bold", size = 12), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c('#99aecf', '#004fcc', '#99aecf', '#004fcc')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

BMI = my_plot_explore("BMI")
resting_heart_rate = my_plot_explore("resting_heart_rate")
HDL = my_plot_explore("HDL")
LDL = my_plot_explore("LDL")
APOEA = my_plot_explore("APOEA")
APOEB = my_plot_explore("APOEB")
ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")


############## for disease #############
######## summarize result on JHPCE #######
result = data.frame(matrix(ncol = 5))
names(result) = c("trait", "method", "ref", "R2", "std")
result_i = 1
for (trait in c("BC", "RA", "IBD", "CAD", "Depression")) {
  for (ref in c("ukbb", "1kg")) {
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs", ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    PRScs_auto = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_auto_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_auto", ref, PRScs_auto$R2, PRScs_auto$std); result_i = result_i + 1
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_proj", ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    if (file.exists(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_auto_", ref, "_a1.0_validation_result_std.txt"))){
      PRScs_auto = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_auto_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
      result[result_i,] = c(trait, "PRScs_proj_auto", ref, PRScs_auto$R2, PRScs_auto$std); result_i = result_i + 1
    }
    else {
      result[result_i,] = c(trait, "PRScs_proj_auto", ref, 0, 0); result_i = result_i + 1
    }
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_diff_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_proj_null", ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    PRScs_auto = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_diff_auto_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_proj_null_auto", ref, PRScs_auto$R2, PRScs_auto$std); result_i = result_i + 1
  }
}
readr::write_tsv(result, "/users/ydun/PRScs_proj_disease_result.txt")


result = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_disease_result.txt", header = TRUE)

library(dplyr)
library(ggplot2)
library(ggpubr)

############### PRScs-proj-null
my_plot_explore = function(pheno_name){
  result = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_disease_result.txt", header = TRUE)
  result = result %>% filter(trait == pheno_name)
  result$method = c('PRScs', 'PRScs_auto', 'PRScs_proj', 'PRScs_proj_auto', 'PRScs_proj_null', 'PRScs_proj_null_auto', '7', '8', '9', '10', '11', '12')
  result$method = factor(result$method, levels = result$method)
  
  result$ref = c(rep('ukbb', 6), rep('1kg', 6))
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
          legend.text = element_text(face = "bold", size = 10), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c('#99aecf','#004fcc',"#c2a482", '#c26802', '#b05fb8', '#b202c2','#99aecf','#004fcc',"#c2a482", '#c26802', '#b05fb8', '#b202c2')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

############### PRScs-proj-null-exclude auto
my_plot_explore = function(pheno_name){
  result = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_disease_result.txt", header = TRUE)
  result = result %>% filter(trait == pheno_name) %>% filter(method != 'PRScs_auto') %>% filter(method != 'PRScs_proj_auto') %>% filter(method != 'PRScs_proj_null_auto')
  result$method = c('PRScs', 'PRScs_proj', 'PRScs_proj_null', '4', '5', '6')
  result$method = factor(result$method, levels = result$method)
  
  result$ref = c(rep('ukbb', 3), rep('1kg', 3))
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
          legend.text = element_text(face = "bold", size = 10), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c('#99aecf','#004fcc',"#c2a482", '#99aecf','#004fcc',"#c2a482")) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

BC = my_plot_explore("BC")
IBD = my_plot_explore("IBD")
CAD = my_plot_explore("CAD")
Depression = my_plot_explore("Depression")
RA = my_plot_explore("RA")
ggarrange(BC + rremove("xlab"), IBD+ rremove("ylab") + rremove("xlab"), CAD+ rremove("ylab") + rremove("xlab"), Depression + rremove("xlab"), RA + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")



############ plot PRS-CS-constraints with proj and regularized LD ###########
##### summarize the result of constraints
result = data.frame(matrix(ncol = 5, nrow = 0))
names(result) = c("trait", "threshold", "ref", "R2", "std")
for (trait in c("BMI", "HDL", "LDL", "resting_heart_rate", "APOEA", "APOEB")) {
  for (ref in c("ukbb", "1kg")) {
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_threshold_", ref, "_a1_validation_result_std.txt"), header = TRUE)
    result = rbind(result, data.frame(trait = trait, threshold = PRScs$threshold, ref = ref, R2 = PRScs$R2, std = PRScs$std))
  }
}
readr::write_tsv(result, "/users/ydun/PRScs_threshold.txt")

result = data.frame(matrix(ncol = 4, nrow = 0))
names(result) = c("trait", "ref", "R2", "std")
for (trait in c("BMI", "HDL", "LDL", "resting_heart_rate", "APOEA", "APOEB")) {
  for (ref in c("ukbb", "1kg")) {
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
    result = rbind(result, data.frame(trait = trait, ref = ref, R2 = PRScs$R2, std = PRScs$std))
  }
}
readr::write_tsv(result, "/users/ydun/PRScs_proj_only.txt")

result = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_threshold.txt", header = TRUE)
library(ggpattern)
my_plot_explore = function(pheno_name){
  result_threshold = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_threshold.txt", header = TRUE) %>% filter(trait == pheno_name)
  result_regul = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_diag_result.txt", header = TRUE) %>% filter(trait == pheno_name)
  result_proj = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_only.txt", header = TRUE) %>% filter(trait == pheno_name)
  
  result = data.frame(matrix(nrow = 0, ncol = 4))
  names(result) = c("method", "ref", "R2", "std"); result_i = 1
  for (ref_i in c("ukbb", "1kg")) {
    tmp = result_threshold %>% filter(trait == pheno_name & ref == ref_i)
    result[result_i,] = c("PRS-CS (Bound:0.01)", ref_i, as.numeric(tmp %>% filter(threshold == 0.01) %>% select(R2)), as.numeric(tmp %>% filter(threshold == 0.01) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS (Bound:0.1)", ref_i, as.numeric(tmp %>% filter(threshold == 0.1) %>% select(R2)), as.numeric(tmp %>% filter(threshold == 0.1) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS (Default)", ref_i, as.numeric(result_regul %>% filter(method == "PRScs" & ref == ref_i) %>% select(R2)), as.numeric(result_regul %>% filter(method == "PRScs" & ref == ref_i) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS (Bound:10)", ref_i, as.numeric(tmp %>% filter(threshold == 10) %>% select(R2)), as.numeric(tmp %>% filter(threshold == 10) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS (Bound:100)", ref_i, as.numeric(tmp %>% filter(threshold == 100) %>% select(R2)), as.numeric(tmp %>% filter(threshold == 100) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS-Projection", ref_i, as.numeric(result_proj %>% filter(ref == ref_i) %>% select(R2)), as.numeric(result_proj %>% filter(ref == ref_i) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS-Regularized", ref_i, as.numeric(result_regul %>% filter(method == "PRScs_proj_diag" & ref == ref_i) %>% select(R2)), as.numeric(result_regul %>% filter(method == "PRScs_proj_diag" & ref == ref_i) %>% select(std))); result_i = result_i + 1
  }
  
  result$method = c(result$method[1:7], '8', '9', '10', '11', "12", "13", "14")
  result$method = factor(result$method, levels = result$method)
  
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (pheno_name == "resting_heart_rate") {
    pheno_name = "RHR"
  }
  result$R2 = as.numeric(result$R2); result$std = as.numeric(result$std)
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method,linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.3, colour = 'black') +
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw() +
    labs(x = "Method", y = expression(bold(R^2)), title = paste0(pheno_name)) +
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(face = "bold", size = 12), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c('#a6f3b8','#4ae673',"#0ecc41", '#079933', '#004d1f','#0ea3cc',"#ccae0e", '#a6f3b8','#4ae673',"#0ecc41", '#079933', '#004d1f','#0ea3cc',"#ccae0e"))+
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

BMI = my_plot_explore("BMI")
resting_heart_rate = my_plot_explore("resting_heart_rate")
HDL = my_plot_explore("HDL")
LDL = my_plot_explore("LDL")
APOEA = my_plot_explore("APOEA")
APOEB = my_plot_explore("APOEB")
plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("~/Documents/PRS_Bridge/figure/PRScs_regularized.pdf", plot, width = 8, height = 6)



############ plot PRS-CS proj in disease data ########
######## summarize result on JHPCE #######
result = data.frame(matrix(ncol = 5))
names(result) = c("trait", "method", "ref", "R2", "std")
result_i = 1
for (trait in c("BC", "Depression", "RA", "CAD", "IBD")) {
  for (ref in c("ukbb", "1kg")) {
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs", ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    PRScs = read.table(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/PRScs_proj_", ref, "_a1.0_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, "PRScs_proj", ref, PRScs$R2, PRScs$std); result_i = result_i + 1
  }
}
readr::write_tsv(result, "/users/ydun/PRScs_proj_disease.txt")

############### PRScs-proj-null
PRS_cs_proj = function(pheno_name){
  result = read.table("~/Documents/PRS_Bridge/result_sd/PRScs_proj_disease.txt", header = TRUE)
  result = result %>% filter(trait == pheno_name) %>% select(-trait)
  result_Bridge_large = read.table(paste0("~/Documents/PRS_Bridge/result_sd/", pheno_name, "/Bridge_large_ukbb_validation_result_std.txt"), header = TRUE)
  result_Bridge_small = read.table(paste0("~/Documents/PRS_Bridge/result_sd/", pheno_name, "/Bridge_small_1kg_validation_result_std.txt"), header = TRUE)
  result = rbind( result[1:2, ],
    data.frame(method = "PRS-Bridge", ref = "ukbb", R2 = result_Bridge_large$R2, std = result_Bridge_large$std),
    result[3:4, ],
    data.frame(method = "PRS-Bridge", ref = "1kg", R2 = result_Bridge_small$R2, std = result_Bridge_small$std)
  )
  result$method = c('PRS-CS', 'PRS-CS-Projection', 'PRS-Bridge', '4', '5', '6')
  result$method = factor(result$method, levels = result$method)
  
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black')+
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw() +
    labs(x = "Method", y = "Transformed AUC", title = paste0(pheno_name)) +
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 1),
          legend.text = element_text(face = "bold", size = 1), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c('#0ecc41','#0ea3cc', '#c26802', '#0ecc41','#0ea3cc', '#c26802')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

BC = PRS_cs_proj("BC")
CAD = PRS_cs_proj("CAD")
Depression = PRS_cs_proj("Depression")
RA = PRS_cs_proj("RA")
IBD = PRS_cs_proj("IBD")
plot = ggarrange(BC+ rremove("xlab"), CAD + rremove("ylab") + rremove("xlab"), Depression + rremove("ylab") + rremove("xlab"), RA+ rremove("xlab"), IBD + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("~/Documents/PRS_Bridge/figure/PRScs_proj_disease.pdf", plot, width = 8, height = 6)


