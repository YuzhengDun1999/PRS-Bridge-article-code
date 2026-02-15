library(dplyr)
library(ggplot2)
library(ggpubr)

#system("python PRS-CS-proj/PRScs_noconverge.py --ref_dir=LD_REFERNCE_DIR --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_noconverge --chrom=22 --write_pst=True --phi=1e-04 --thin=1 --n_iter=800 --n_burnin=1 --seed=1")

############# plot non-convergence-result ###############
PRScs = bigreadr::fread2("data/coef_noconverge_pst_eff_a1_b0.5_phi1e-04_chr22.txt")
coef = as.matrix(PRScs[,6:804]); rownames(coef) = NULL; colnames(coef) = NULL
dat = data.frame(iter = 1:ncol(coef), y = coef[8131,], SNP = PRScs[8131,2])
dat = rbind(dat, data.frame(iter = 1:ncol(coef), y = coef[1622,], SNP = PRScs[1622,2]))
dat = rbind(dat, data.frame(iter = 1:ncol(coef), y = coef[849,], SNP = PRScs[849,2]))

y_range = range(10^{-9},10^13)
dat$SNP = factor(dat$SNP, levels = c(PRScs[8131,2], PRScs[1622,2], PRScs[849,2]))

plot1 = ggplot(data = dat %>% filter(iter <= 700)) + geom_point(mapping = aes(x = iter, y = (abs(y)), color = SNP), size = 0.8) + theme_bw() + 
  theme_bw() + labs(x = "MCMC Iteration", y = "Coefficient Magnitude", title = "Without Projection of Summary Statistics") + 
  theme(plot.title = element_text(face="bold", size = 16, hjust = 0.5), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 18), legend.title = element_text(size = 18),
        axis.title=element_text(size=18), legend.position = "bottom", axis.text = element_text(size = 18)) +
  scale_y_continuous(trans = 'log10', limits = y_range, labels = scales::trans_format('log10', scales::math_format(10^.x))) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  guides(color = guide_legend(override.aes = list(size = 6)))

#system("python PRS-CS-proj/PRScs_proj.py --ref_dir=LD_REFERNCE_DIR--bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_proj --chrom=22 --write_pst=True --phi=1e-04 --thin=1 --n_iter=800 --n_burnin=1 --seed=1 --eigenval_rm=0.2")

############# plot projection result ###############
PRScs_proj = bigreadr::fread2("data/coef_proj_pst_eff_a1_b0.5_phi1e-04_chr22.txt")
coef = as.matrix(PRScs_proj[,6:804]); rownames(coef) = NULL; colnames(coef) = NULL

### Create same plots for different SNPs
dat = data.frame(iter = 1:ncol(coef), y = coef[8131,], SNP = PRScs[8131,2])
dat = rbind(dat, data.frame(iter = 1:ncol(coef), y = coef[1622,], SNP = PRScs[1622,2]))
dat = rbind(dat, data.frame(iter = 1:ncol(coef), y = coef[849,], SNP = PRScs[849,2]))
dat$SNP = factor(dat$SNP, levels = c(PRScs[8131,2], PRScs[1622,2], PRScs[849,2]))

plot2 = ggplot(data = dat %>% filter(iter <= 700)) + geom_point(mapping = aes(x = iter, y = (abs(y)), color = SNP), size = 0.8) + theme_bw() + 
  theme_bw() + labs(x = "MCMC Iteration", y = "Coefficient Magnitude", title = "With Projection of Summary Statistics") + 
  theme(plot.title = element_text(face="bold", size = 16, hjust = 0.5), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 18), legend.title = element_text(size = 18),
        axis.title=element_text(size=18), legend.position = "bottom", axis.text = element_text(size = 18)) +
  scale_y_continuous(trans = 'log10', limits = y_range, labels = scales::trans_format('log10', scales::math_format(10^.x))) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  guides(color = guide_legend(override.aes = list(size = 6)))
plot_all = ggarrange(plot1, plot2 + rremove("ylab"), ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
ggsave("PRScs_traceplot_seed.pdf", plot_all, width = 12, height = 6)


