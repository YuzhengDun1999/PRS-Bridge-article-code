########### This script is used to generate and submit script files for PRS-CS
library(dplyr)
trait = "Example"

system(paste0("mkdir -p ", trait, "/1kg/PRScs/run_sh"))
system(paste0("mkdir -p ", trait, "/ukbb/PRScs/run_sh"))
system(paste0("mkdir -p ", trait, "/1kg/PRScs/result"))
system(paste0("mkdir -p ", trait, "/ukbb/PRScs/result"))

############### generate run_sh file for PRScs using ukbb ref################
a_list = c(1) # controlling prior sparsity assumption. Default: 1
phi_list = c(0, 1, 0.01, 0.0001, 0.000001) # controlling prior sparsity assumption. Recommended: 0.2, 0.4, 0.6, 0.8

name_list = c('_auto', '_1', '_e2', '_e4', '_e6') # corresponding to the name of phi_list
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06') # corresponding to the name of output file of name_list
run_sh_list = c()
for (chr in c(1:22)) {
  for (phi_i in 1:length(phi_list)) {
    for (a in a_list) {
      phi = phi_list[phi_i]
      name = name_list[phi_i]
      run_sh = paste0(trait, '/ukbb/PRScs/run_sh/run_chr', chr, name, '_a', a,'.sh')
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sumdat = paste0(trait, '/data/chr', chr, '/PRScs_sumdat.txt')
      ldsc = paste0(trait, '/data/chr', chr, '/ldsc-hm3.txt')
      sumdat_tmp = bigreadr::fread2(ldsc)
      N = as.integer(floor(mean(sumdat_tmp$N)))
      output_dir = paste0(trait,'/ukbb/PRScs/result/chr', chr)
      
      sent_run <- paste0('python PRScs/PRScs.py --ref_dir=ldblk_ukbb_eur --bim_prefix=chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=', #change the path to the LD reference data and validation bim file here
                         output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi)
      if(phi == 0){
        sent_run <- paste0('python PRScs/PRScs.py --ref_dir=ldblk_ukbb_eur --bim_prefix=chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=', #change the path to the LD reference data and validation bim file here
                           output_dir, ' --a=', a, ' --chrom=', chr)
      }
      cat("module load anaconda \n")
      cat(sent_run, "\n")
      bfile = paste0('chr', chr) #change the path to the validation bim file here
      name1 = name_list_prs[phi_i]
      result_dir = paste0(trait, '/ukbb/PRScs/result/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name1, '_chr', chr)
      prs_output = paste0(result_dir, '_prs.txt')
      prscode = paste(paste0('plink'),
                      paste0('--score ', result_dir, '.txt'),
                      paste0(' 2 5 6 sum --bfile ', bfile),
                      " --threads 1",
                      paste0(' --out ', prs_output))
      cat(prscode, "\n")
      sink()
    }
  }
}

run_all = paste0(trait, "/ukbb/PRScs/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch  --mem=', mem,'G ', run_sh)
  cat(submit_sh, "\n")
}
sink()
system(paste0("bash ", run_all))


############### generate run_sh file for PRScs using 1kg################
a_list = c(0.5, 1.0, 1.5) # controlling prior sparsity assumption. Default: 1
phi_list = c(0, 1, 0.01, 0.0001, 0.000001) # controlling prior sparsity assumption. Recommended: 0.2, 0.4, 0.6, 0.8

name_list = c('_auto', '_1', '_e2', '_e4', '_e6') # corresponding to the name of phi_list
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06') # corresponding to the name of output file of name_list
run_sh_list = c()
for (chr in c(1:22)) {
  for (phi_i in 1:length(phi_list)) {
    for (a in a_list) {
      phi = phi_list[phi_i]
      name = name_list[phi_i]
      run_sh = paste0(trait, '/1kg/PRScs/run_sh/run_chr', chr, name, '_a', a,'.sh')
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sumdat = paste0(trait, '/data/chr', chr, '/PRScs_sumdat.txt')
      ldsc = paste0(trait, '/data/chr', chr, '/ldsc-hm3.txt')
      sumdat_tmp = bigreadr::fread2(ldsc)
      N = as.integer(floor(mean(sumdat_tmp$N)))
      output_dir = paste0(trait,'/1kg/PRScs/result/chr', chr)
      
      sent_run <- paste0('python PRScs/PRScs.py --ref_dir=1kg/ldblk_1kg_eur --bim_prefix=chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=', #change the path to the LD reference data and validation bim file here
                         output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi)
      if(phi == 0){
        sent_run <- paste0('python PRScs/PRScs.py --ref_dir=1kg/ldblk_1kg_eur --bim_prefix=chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=', #change the path to the LD reference data and validation bim file here
                           output_dir, ' --a=', a, ' --chrom=', chr)
      }
      cat("module load anaconda \n")
      cat(sent_run, "\n")
      bfile = paste0('chr', chr)
      name1 = name_list_prs[phi_i]
      result_dir = paste0(trait, '/1kg/PRScs/result/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name1, '_chr', chr)
      prs_output = paste0(result_dir, '_prs.txt')
      prscode = paste(paste0('plink'),
                      paste0('--score ', result_dir, '.txt'),
                      paste0(' 2 5 6 sum --bfile ',bfile),
                      " --threads 1",
                      paste0(' --out ', prs_output))
      cat(prscode, "\n")
      sink()
    }
  }
}

run_all = paste0(trait, "/1kg/PRScs/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch  --mem=', mem,'G ', run_sh)
  cat(submit_sh, "\n")
}
sink()
setwd(paste0(trait, "/1kg/PRScs/run_sh/"))
system(paste0("bash ", run_all))


