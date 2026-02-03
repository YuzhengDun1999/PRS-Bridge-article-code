library(dplyr)
library(ggplot2)
library(ggpubr)
result = data.frame()
for (method_i in c("Bridge_small", "ldpred2", "PRScs")) {
  for (rho_i in 1:3) {
    for (GA_i in c(1, 4, 5)) {
      result_tmp = read.table(paste0(trait, "/1kg/", method, "/result/", "-rho", rho_i, "-size4-rep1-GA", GA_i, "_validation_result_std.txt"), header = TRUE)
      result_tmp = result_tmp %>% select(R2, std) %>% mutate(method = method_i, rho = rho_i, GA = GA_i)
      result = rbind(result, result_tmp)
    }
  }
}
result[result$GA==1, 'GA'] = 'Strong'
result[result$GA==4, 'GA'] = 'No'
result[result$GA==5, 'GA'] = 'Mild'
result[result$method == "Bridge_small", 'method'] = "PRS-Bridge"
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
ggsave("simulation_sd.pdf", plot, width = 6, height = 5)
