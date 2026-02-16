########### This script is used to generate and submit script files for PRS-CS-threshold
temp <- commandArgs(TRUE)
trait = temp[1] # phenotype name, same for the whole pipeline

a_list = c(1) # controlling prior sparsity assumption. Default: 1, can be 0.5, 1.0, 1.5
phi_list = c(0, 1, 0.01, 0.0001, 0.000001) # controlling prior sparsity assumption. Recommended: 0, 1, 0.01, 0.0001, 0.000001
threshold_list = c(0.01, 0.1, 1, 10, 100) # controlling ad hoc threshold used in PRS-CS. Recommended:0.01, 0.1, 1, 10, 100

name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
for (ref in c('ukbb', '1kg')) {
  system(paste0("mkdir -p ", trait, "/", ref, "/PRScs_threshold/run_sh"))
  system(paste0("mkdir -p ", trait, "/", ref, "/PRScs_threshold/result"))
  for (chr in c(1:22)) {
    for (phi_i in 1:length(phi_list)) {
      for (threshold_i in threshold_list) {
        for (a in a_list) {
          phi = phi_list[phi_i]
          name = name_list[phi_i]
          run_sh = paste0( trait, '/', ref, '/PRScs_threshold/run_sh/run_chr', chr, name, '_threshold', threshold_i, '_a', a, '.sh')
          run_sh_list = c(run_sh_list, run_sh)
          sink(run_sh)
          cat('#!/bin/bash', "\n")
          cat(paste0('mkdir -p ', trait,'/', ref, '/PRScs_threshold/result/threshold', threshold_i), "\n")
          sumdat = paste0(trait, '/data/chr', chr, '/PRScs_sumdat.txt')
          ldsc = paste0(trait, '/data/chr', chr, '/ldsc-hm3.txt')
          sumdat_tmp = bigreadr::fread2(ldsc)
          N = as.integer(floor(mean(sumdat_tmp$N)))
          output_dir = paste0(trait,'/', ref, '/PRScs_threshold/result/threshold', threshold_i,'/chr', chr, '')
          
          sent_run <- paste0('python PRScs/PRScs_threshold.py --ref_dir=data/', ref, '/ldblk_', ref, '_eur --bim_prefix=chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                             output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --threshold=', threshold_i)
          if(phi == 0){
            sent_run <- paste0('python PRScs/PRScs_threshold.py --ref_dir=data/', ref, '/ldblk_', ref, '_eur --bim_prefix=chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                               output_dir, ' --a=', a, ' --chrom=', chr, ' --threshold=', threshold_i)
          }
          cat("module load anaconda \n")
          cat(sent_run, "\n")
          sink()
        }
      }
    }
  }
  
  run_all = paste0(trait, "/", ref, "/PRScs_threshold/run_sh/run_all.sh")
  sink(run_all)
  mem = 3
  for (run_sh_i in 1:length(run_sh_list)) {
    run_sh = run_sh_list[run_sh_i]
    submit_sh = paste0('sbatch  --mem=', mem,'G ', run_sh)
    cat(submit_sh, "\n")
  }
  sink()
  system(paste0("bash ", run_all))
}