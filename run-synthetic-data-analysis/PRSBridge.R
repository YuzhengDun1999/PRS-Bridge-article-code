####### This script is used for generating script files for running PRS-Bridge in simulation studies
library(dplyr)
temp <- commandArgs(TRUE)
rho = as.numeric(temp[1]) ###### taking values in 1, 2, 3
GA = as.numeric(temp[2]) ###### taking values in 1, 4, 5

trait = paste0('rho', rho, 'GA', GA)
system(paste0("mkdir -p ", trait, "/1kg/Bridge_small/run_sh/chr{1..22}"))

####### generate run_sh file for PRSBridge using 1kg small ref #######
alpha_list = c(0.5, 0.25, 0.125) # controlling prior sparsity assumption. 'auto' is estimating alpha during MCMC. Recommended for tuning: 0.5, 0.25, 0.125
percent_list = c(0.2, 0.4, 0.6, 0.8) # controlling percent of eigenvalue removed. Recommended: 0.2, 0.4, 0.6, 0.8
for (chr in 1:22) {
  run_sh_list = c()
  for (alpha_i in 1:length(alpha_list)) {
    for (percent_i in 1:length(percent_list)) {
      alpha = alpha_list[alpha_i]
      percent = percent_list[percent_i]
      run_sh = paste0("run_percent", percent, "_alpha", alpha, ".sh")
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sent_dir <- paste0("mkdir -p ", output_dir)
      sumdat = paste0(trait, '/sumdat/chr', chr, '/hm3_sumdat.txt')
      output_dir = paste0(trait, "/1kg/Bridge_small/chr", chr, "/percent", percent, "/alpha", alpha)
      sent_run <- paste0("python PRSBridge/PRSBridge.py --percent ", percent, " --chr ", chr, " --alpha ", alpha,  " --ref ref_1kg/chr",  " --sumdat ", sumdat,  " --h2 ", 0,  " --h2_se ", 0 ##### change your path to LD reference data here
                         ,  " --method ", 'cg',  " --output ", output_dir)
      cat(sent_dir, "\n")
      cat("module load anaconda \n")
      cat(sent_run, "\n")
      sink()
    }
  }
  
  run_all = paste0("run_sh/chr", chr, "/run_all.sh")
  sink(run_all)
  mem = 2
  if(chr<21) mem = 5
  for (run_sh_i in 1:length(run_sh_list)) {
    run_sh = run_sh_list[run_sh_i]
    submit_sh = paste0('sbatch --mem=', mem,'G ', run_sh)
    cat(submit_sh, "\n")
  }
  sink()
  system(paste0("bash ", run_all))
}

