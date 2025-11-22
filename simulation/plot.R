library(dplyr)
library(ggplot2)
library(ggpubr)
result = data.frame()
for (method_i in c("Bridge", "ldpred2", "Lassosum", "PRScs")) {
  for (rho_i in 1:3) {
    for (GA_i in c(1, 4, 5)) {
      result_tmp = read.table(paste0("~/Documents/PRS_Bridge/simulation/result/", method_i, "-rho", rho_i, "-size4-rep1-GA", GA_i, "_validation_result_std.txt"), header = TRUE)
      result_tmp = result_tmp %>% select(R2, std) %>% mutate(method = method_i, rho = rho_i, GA = GA_i)
      result = rbind(result, result_tmp)
    }
  }
}
result[result$GA==1, 'GA'] = 'Strong'
result[result$GA==4, 'GA'] = 'No'
result[result$GA==5, 'GA'] = 'Mild'
result[result$method == "Bridge", 'method'] = "PRS-Bridge"
result[result$method == "ldpred2", 'method'] = "LDpred2"
result[result$method == "PRScs", 'method'] = "PRS-CS"
result[result$method == "Lassosum", 'method'] = "lassosum"
result$method = factor(result$method, levels = c("lassosum", "LDpred2", "PRS-CS", "PRS-Bridge"))
result$GA = factor(result$GA, levels = c('No', 'Mild', 'Strong'))

plot_list = list()
result = result %>% rename(Method = method)
for (rho1 in 1:3) {
  dat = result %>% filter(rho == rho1)
  if(rho1==1) title1 = expression(bold(paste("Causal SNP Proportion: ", 1%*%10^{-2})))
  if(rho1==2) title1 = expression(bold(paste("Causal SNP Proportion: ", 1%*%10^{-3})))
  if(rho1==3) title1 = expression(bold(paste("Causal SNP Proportion: ", 5%*%10^{-4})))
  
  temp = ggplot(data = dat, mapping = aes(x = GA, y = R2, fill = Method)) + geom_bar(stat = 'identity', position = position_dodge(width = 0.7), width = 0.7) + theme_bw() +
    labs(x = "Negative selection", y = expression(bold(R^2)), title = title1)  +
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), position = position_dodge(width = 0.7), width = 0.15) +
    theme(axis.title = element_text(face="bold", size = 12),
          axis.text = element_text(face="bold", size = 10),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(face = "bold", size = 14), 
          plot.title = element_text(face="bold", size = 13, hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    scale_fill_manual(values=c("#e3c278", "#0250c9", "#0ecc41", "#c26802"))
  plot_list[[rho1]] = temp
}

############ no Lassosum ########
result = data.frame()
for (method_i in c("Bridge", "ldpred2", "PRScs")) {
  for (rho_i in 1:3) {
    for (GA_i in c(1, 4, 5)) {
      result_tmp = read.table(paste0("~/Documents/PRS_Bridge/simulation/result/", method_i, "-rho", rho_i, "-size4-rep1-GA", GA_i, "_validation_result_std.txt"), header = TRUE)
      result_tmp = result_tmp %>% select(R2, std) %>% mutate(method = method_i, rho = rho_i, GA = GA_i)
      result = rbind(result, result_tmp)
    }
  }
}
result[result$GA==1, 'GA'] = 'Strong'
result[result$GA==4, 'GA'] = 'No'
result[result$GA==5, 'GA'] = 'Mild'
result[result$method == "Bridge", 'method'] = "PRS-Bridge"
result[result$method == "ldpred2", 'method'] = "LDpred2"
result[result$method == "PRScs", 'method'] = "PRS-CS"
result$method = factor(result$method, levels = c("LDpred2", "PRS-CS", "PRS-Bridge"))
result$GA = factor(result$GA, levels = c('No', 'Mild', 'Strong'))

plot_list = list()
result = result %>% rename(Method = method)


for (rho1 in 1:3) {
  dat = result %>% filter(rho == rho1)
  if(rho1==1) title1 = expression(bold(paste("Causal SNP Proportion: ", 1%*%10^{-2})))
  if(rho1==2) title1 = expression(bold(paste("Causal SNP Proportion: ", 1%*%10^{-3})))
  if(rho1==3) title1 = expression(bold(paste("Causal SNP Proportion: ", 5%*%10^{-4})))
  
  temp = ggplot(data = dat, mapping = aes(x = GA, y = R2, fill = Method)) + geom_bar(stat = 'identity', position = position_dodge(width = 0.7), width = 0.7) + theme_bw() +
    labs(x = "Negative selection", y = expression(bold(R^2)), title = title1)  +
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), position = position_dodge(width = 0.7), width = 0.15) +
    theme(axis.title = element_text(face="bold", size = 12),
          axis.text = element_text(face="bold", size = 10),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(face = "bold", size = 14), 
          plot.title = element_text(face="bold", size = 13, hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    scale_fill_manual(values=c("#0250c9", "#0ecc41", "#c26802"))
  plot_list[[rho1]] = temp
}
plot = ggarrange(plot_list[[1]]+ rremove("xlab"), plot_list[[2]]+ rremove("xlab"), plot_list[[3]], nrow = 3,common.legend = TRUE, legend ="right")
ggsave("~/Documents/PRS_Bridge/prs-bridge-manuscript/JASA_application/simulation_sd.pdf", plot, width = 6, height = 5)


###### plot tuning setting
library(dplyr)
library(ggplot2)
library(ggpubr)

result = data.frame()
for (rho_i in 1:3) {
  for (GA_i in c(1, 4, 5)) {
    result_tmp = read.table(paste0("~/Documents/PRS_Bridge/simulation/result/Bridge-rho", rho_i, "-size4-rep1-GA", GA_i, "_all_result_std.txt"), header = TRUE)
    result_tmp = result_tmp %>% select(R2, std) %>% mutate(method = method_i, rho = rho_i, GA = GA_i)
    result = rbind(result, result_tmp)
  }
}
  
result[result$GA==1, 'GA'] = 'Strong'
result[result$GA==4, 'GA'] = 'No'
result[result$GA==5, 'GA'] = 'Mild'
result[result$method == "Bridge", 'method'] = "PRS-Bridge"
result[result$method == "ldpred2", 'method'] = "LDpred2"
result[result$method == "PRScs", 'method'] = "PRS-CS"
result[result$method == "Lassosum", 'method'] = "lassosum"
result$method = factor(result$method, levels = c("lassosum", "LDpred2", "PRS-CS", "PRS-Bridge"))
result$GA = factor(result$GA, levels = c('No', 'Mild', 'Strong'))

plot_tuning = function(rho_i, GA_i=4){
  result = read.table(paste0("~/Documents/PRS_Bridge/simulation/result/Bridge-rho", rho_i, "-size4-rep1-GA", GA_i, "_all_result_std.txt"), header = TRUE)
  result = result %>% filter(percent != 0 ) %>% filter(alpha != 0.0625)
  result_sd = result %>% group_by(alpha, percent) %>% summarise( mean_R2 = mean(tuning_R2), se_R2   = sd(tuning_R2) / sqrt(n()), .groups = "drop") %>% as.data.frame()
  if(rho_i==1) title1 = expression(bold(paste("Causal SNP Proportion: ", 1%*%10^{-2})))
  if(rho_i==2) title1 = expression(bold(paste("Causal SNP Proportion: ", 1%*%10^{-3})))
  if(rho_i==3) title1 = expression(bold(paste("Causal SNP Proportion: ", 5%*%10^{-4})))
  
  plot1 = ggplot(result_sd, aes(x = factor(percent), y = mean_R2, group = factor(alpha), color = factor(alpha))) +
    geom_point(position = position_dodge(width = 0.3), size = 2) +
    geom_line(position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = mean_R2 - 1.96 * se_R2, ymax = mean_R2 + 1.96 * se_R2),
                  width = 0.2, position = position_dodge(width = 0.3)) + theme_bw() +
    labs(x = "Percent", y = expression(bold(R^2)), color = expression(bold(alpha)), title = title1) +
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
plot1 = plot_tuning(1, 1)
plot2 = plot_tuning(2, 1)
plot3 = plot_tuning(3, 1)
ggarrange(plot1, plot2, plot3, ncol=2, nrow=2, common.legend = TRUE, legend="right")

plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL, APOEA + rremove("ylab"), APOEB+ rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="right")
ggsave("~/Documents/PRS_Bridge/prs-bridge-manuscript/JASA_application/tuning_figure.pdf", plot, width = 8, height = 6)

