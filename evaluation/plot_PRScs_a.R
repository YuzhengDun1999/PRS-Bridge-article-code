#################### plot PRS-CS with different a in supplementary figure S5###########
#################### To generate following plot, please set a_list = c(0.5, 1.0, 1.5), ref = "ukbb" and "1kg" in run-methods and evaluation for PRS-CS methods


result = data.frame(matrix(ncol = 5))
names(result) = c("trait", "a", "ref", "R2", "std"); result_i = 1
for (trait in c("BMI", "HDL", "LDL", "resting_heart_rate", "APOEA", "APOEB")) {
  for (ref in c("ukbb", "1kg")) {
    PRScs_a0.5 = read.table(paste0(trait, "/result/PRScs_a0.5_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, 0.5, ref, PRScs_a0.5$R2, PRScs_a0.5$std); result_i = result_i + 1
    PRScs = read.table(paste0(trait, "/result/PRScs_a1_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, 1, ref, PRScs$R2, PRScs$std); result_i = result_i + 1
    PRScs_a1.5 = read.table(paste0(trait, "/result/PRScs_a1.5_", ref, "_validation_result_std.txt"), header = TRUE)
    result[result_i,] = c(trait, 1.5, ref, PRScs_a1.5$R2, PRScs_a1.5$std); result_i = result_i + 1
  }
}

my_plot_new = function(pheno_name){
  method = c('a=0.5', 'a=1', 'a=1.5', 'a=0.5 1kg', 'a=1 1kg', 'a=1.5 1kg')  
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
ggsave("PRScs_a.pdf", plot, width = 8, height = 6)