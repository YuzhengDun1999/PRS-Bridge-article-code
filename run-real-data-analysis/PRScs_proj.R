########### This script is used to generate and submit script files for PRS-CS-Projection
temp <- commandArgs(TRUE)
trait = temp[1] # phenotype name, same for the whole pipeline

a_list = c(1) # controlling prior sparsity assumption. Default: 1, can be 0.5, 1.0, 1.5
phi_list = c(0, 1, 0.01, 0.0001, 0.000001) # controlling prior sparsity assumption. Recommended: 0, 1, 0.01, 0.0001, 0.000001
percent_list = c(0.2, 0.4, 0.6, 0.8) # controlling percent of eigenvalue removed. Recommended: 0.2, 0.4, 0.6, 0.8

name_list = c('_auto', '_1', '_e2', '_e4', '_e6') # corresponding to the name of phi_list
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06') # corresponding to the name of output file of name_list
run_sh_list = c()
for (ref in c("1kg", "ukbb")) {
  
  system(paste0("mkdir -p ", trait, "/", ref, "/PRScs_proj/run_sh"))
  system(paste0("mkdir -p ", trait, "/", ref, "/PRScs_proj/result"))
  for (chr in c(1:22)) {
    for (phi_i in 1:length(phi_list)) {
      for (percent_i in percent_list) {
        for (a in a_list) {
          phi = phi_list[phi_i]
          name = name_list[phi_i]
          run_sh = paste0( trait, '/', ref, '/PRScs_proj/run_sh/run_chr', chr, name, '_percent', percent_i, '_a', a, '.sh')
          run_sh_list = c(run_sh_list, run_sh)
          sink(run_sh)
          cat('#!/bin/bash', "\n")
          cat(paste0('mkdir -p ', trait,'/', ref, '/PRScs_proj/result/percent', percent_i), "\n")
          sumdat = paste0(trait, '/data/chr', chr, '/PRScs_sumdat.txt')
          ldsc = paste0(trait, '/data/chr', chr, '/ldsc-hm3.txt')
          sumdat_tmp = bigreadr::fread2(ldsc)
          N = as.integer(floor(mean(sumdat_tmp$N)))
          output_dir = paste0(trait,'/', ref, '/PRScs_proj/result/percent', percent_i,'/chr', chr, '')
          
          sent_run <- paste0('python PRScs/PRScs_proj.py --ref_dir=data/', ref, '/ldblk_', ref, '_eur --bim_prefix=chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                             output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --eigenval_rm=', percent_i)
          if(phi == 0){
            sent_run <- paste0('python PRScs/PRScs_proj.py --ref_dir=data/', ref, '/ldblk_', ref, '_eur --bim_prefix=chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                               output_dir, ' --a=', a, ' --chrom=', chr, ' --eigenval_rm=', percent_i)
          }
          cat("module load anaconda \n")
          cat(sent_run, "\n")
          sink()
        }
      }
    }
  }
  
  run_all = paste0(trait, "/", ref, "/PRScs_proj/run_sh/run_all.sh")
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