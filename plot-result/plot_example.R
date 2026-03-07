library(dplyr)
library(ggplot2)
library(ggpubr)

######## To generate following plot, please run all files in run-methods and evaluation with suggested default parameters, reference LD needs to be from "ukbb" and "1kg".

##################### plot all methods, generate Figure 4 and 5 in main manuscript ##################
method = c('Lassosum (Small-block)', 'LDpred2 (Banded)', 'LDpred2 (Large-block)', 'PRS-CS (Small-block)', 'PRS-Bridge (Small-block)', 'PRS-Bridge (Large-block)', 'Lassosum2_1kg', 'ldpred2_banded_1kg', 'PRScs_1kg', 'PRSBridge_small_block_1kg', 'PRSBridge_large_block_1kg')  

##### Input trait is the same as trait used in run-methods and evaluation
my_plot_new = function(trait){
  result = read.table(paste0(trait, '/Lassosum_Nukbb_validation_result_std.txt'), header = TRUE)
  result = data.frame(method = method[1], R2 = result$R2, std = result$std)
  result_i = 1
  validation_result = read.table(paste0(trait, '/ldpred2_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/ldpred2_Nblk_provided_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/PRScs_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/Bridge_small_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/Bridge_large_ukbb_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/Lassosum_N1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/ldpred2_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/PRScs_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/Bridge_small_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  validation_result = read.table(paste0(trait, '/Bridge_large_1kg_validation_result_std.txt'), header = TRUE); result_i = result_i + 1
  result[result_i, 1] = method[result_i]; result[result_i, 2] = validation_result$R2; result[result_i, 3] = validation_result$std
  
  result$method = factor(result$method, levels = method)
  
  result$ref = c(rep('ukbb', 6), rep('1kg', 5))
  result$ref = factor(result$ref, levels = c('ukbb', '1kg'))
  if (trait == "resting_heart_rate") {
    trait = "RHR"
  }
  plot1 = ggplot(data = result, mapping = aes(x = method, y = R2, fill = method, linetype = ref)) + 
    geom_bar(stat = 'identity', position = 'dodge', linewidth = 0.7, colour = 'black')+
    geom_errorbar(aes(ymin = R2 - 1.96 * std, ymax = R2 + 1.96 * std), width = 0.15) + theme_bw() +
    labs(x = "Method", title = paste0(trait)) +
    theme(axis.text.x=element_blank(), 
          axis.text.y=element_text(face="bold", size = 14),
          axis.title = element_text(face="bold", size = 14),
          axis.text = element_text(face="bold", size = 14),
          legend.title = element_text(face = "bold", size = 1),
          legend.text = element_text(face = "bold", size = 1), 
          plot.title = element_text(face="bold", size = 18, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c("Lassosum (Small-block)" = "#e3c278", "LDpred2 (Banded)" = "#99aecf", "LDpred2 (Large-block)" = "#0250c9",
                               "PRS-CS (Small-block)" = "#0ecc41", "PRS-Bridge (Small-block)" = "#c2a482", "PRS-Bridge (Large-block)" = "#c26802",
                               "Lassosum2_1kg" = "#e3c278", "ldpred2_banded_1kg" = "#99aecf", "PRScs_1kg" = "#0ecc41",
                               "PRSBridge_small_block_1kg"  = "#c2a482", "PRSBridge_large_block_1kg"  = "#c26802"
    ), breaks = c('Lassosum (Small-block)', 'LDpred2 (Banded)', 'LDpred2 (Large-block)', 'PRS-CS (Small-block)', 'PRS-Bridge (Small-block)', 'PRS-Bridge (Large-block)')) +
    guides(linetype = guide_legend(override.aes = list(fill = NA, col = "black")))
  return(plot1)
}

example1 = my_plot_new("example") + labs(y = expression(bold(R^2)))
example2 = my_plot_new("example") + labs(y = expression(bold(R^2)))
example3 = my_plot_new("example") + labs(y = expression(bold(R^2)))
example4 = my_plot_new("example") + labs(y = expression(bold(R^2)))
example5 = my_plot_new("example") + labs(y = expression(bold(R^2)))
example6 = my_plot_new("example") + labs(y = expression(bold(R^2)))

plot = ggarrange(BMI + rremove("xlab"), resting_heart_rate+ rremove("ylab") + rremove("xlab"), HDL+ rremove("ylab") + rremove("xlab"), LDL + rremove("xlab"), APOEA + rremove("xlab") + rremove("ylab"), APOEB + rremove("xlab") + rremove("ylab"), ncol=3, nrow=2, common.legend = TRUE, legend="left")
ggsave("example_figure.pdf", plot, width = 8, height = 6)
