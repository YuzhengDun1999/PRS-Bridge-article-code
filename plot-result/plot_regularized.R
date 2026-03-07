###### plot Figure S7 in supplementary material
library(dplyr)
library(ggplot2)
library(ggpubr)

###### To generate this plot, Please run PRScs_threshold.R, PRScs_proj.R, and PRSs_regularized.R in run-methods section using ukbb and 1kg as reference data and do evaluation for each method in evaluation section

####### Input trait is the same as trait used in run-methods and evaluation
my_plot = function(trait){
  result = data.frame(matrix(nrow = 0, ncol = 4))
  names(result) = c("method", "ref", "R2", "std"); result_i = 1
  for (ref_i in c("ukbb", "1kg")) {
    result_threshold = read.table(paste0(trait, "/result/PRScs_threshold_", ref_i, '_a1_validation_result_std.txt'), header = TRUE)
    result_regul = read.table(paste0(trait, "/result/PRScs_regularized_", ref_i, '_a1_validation_result_std.txt'), header = TRUE)
    result_proj = read.table(paste0(trait, "/result/PRScs_proj_", ref_i, '_a1_validation_result_std.txt'), header = TRUE)
    
    result[result_i,] = c("PRS-CS (Bound:0.01)", ref_i, as.numeric(result_threshold %>% filter(threshold == 0.01) %>% select(R2)), as.numeric(tmp %>% filter(threshold == 0.01) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS (Bound:0.1)", ref_i, as.numeric(result_threshold %>% filter(threshold == 0.1) %>% select(R2)), as.numeric(tmp %>% filter(threshold == 0.1) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS (Default)", ref_i, as.numeric(result_threshold %>% filter(threshold == 1) %>% select(R2)), as.numeric(result_regul %>% filter(method == "PRScs" & ref == ref_i) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS (Bound:10)", ref_i, as.numeric(result_threshold %>% filter(threshold == 10) %>% select(R2)), as.numeric(tmp %>% filter(threshold == 10) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS (Bound:100)", ref_i, as.numeric(result_threshold %>% filter(threshold == 100) %>% select(R2)), as.numeric(tmp %>% filter(threshold == 100) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS-Projection", ref_i, as.numeric(result_proj %>% filter(ref == ref_i) %>% select(R2)), as.numeric(result_proj %>% filter(ref == ref_i) %>% select(std))); result_i = result_i + 1
    result[result_i,] = c("PRS-CS-Regularized", ref_i, as.numeric(result_regul %>% filter(method == "PRScs_proj_diag" & ref == ref_i) %>% select(R2)), as.numeric(result_regul %>% filter(method == "PRScs_proj_diag" & ref == ref_i) %>% select(std))); result_i = result_i + 1
  }
  
  result$method = c(result$method[1:7], '8', '9', '10', '11', "12", "13", "14")
  result$method = factor(result$method, levels = result$method)
  
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (trait == "resting_heart_rate") {
    trait = "RHR"
  }
  result$R2 = as.numeric(result$R2); result$std = as.numeric(result$std)
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method,linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.3, colour = 'black') +
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw() +
    labs(x = "Method", y = expression(bold(R^2)), title = paste0(trait)) +
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(face = "bold", size = 12), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c("PRS-CS (Bound:0.01)" = "#a6f3b8", "PRS-CS (Bound:0.1)" = "#4ae673",
                               "PRS-CS (Default)" = "#0ecc41", "PRS-CS (Bound:10)" = "#079933", "PRS-CS (Bound:100)" = "#004d1f",
                               "PRS-CS-Projection" = "#0ea3cc", "PRS-CS-Regularized" = "#ccae0e", "8" = "#a6f3b8", "9" = "#4ae673",
                               "10" = "#0ecc41", "11" = "#079933", "12" = "#004d1f", "13" = "#0ea3cc", "14" = "#ccae0e"
    ), breaks = c('PRS-CS (Bound:0.01)', 'PRS-CS (Bound:0.1)', 'PRS-CS (Default)', 'PRS-CS (Bound:10)', 'PRS-CS (Bound:100)', 'PRS-CS-Projection', 'PRS-CS-Regularized')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

BMI = my_plot("BMI")
resting_heart_rate = my_plot("resting_heart_rate")
HDL = my_plot("HDL")
LDL = my_plot("LDL")
APOEA = my_plot("APOEA")
APOEB = my_plot("APOEB")
plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("PRScs_regularized.pdf", plot, width = 8, height = 6)
