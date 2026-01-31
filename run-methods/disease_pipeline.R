library(readr)
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
trait = "IBD"

################### preprocess summary level ################
library(dplyr)
# sumdat = bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/GWAS/sumdat_Rcov.txt'))
sumdat = sumdat %>% mutate(N = round(N))
ldsc_sumdat_all_hm3 = data.frame()
for (chr in 1:22) {
  out_dir = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/')
  sumdat_hm3 = sumdat[sumdat$CHR == chr, ]
  sumdat_hm3 = sumdat_hm3 %>% filter(((BETA/SE)^2/N <= 80)) %>%
    select(CHR, POS, SNP_ID, REF, ALT, REF_FRQ, PVAL, BETA, SE, N) %>%
    filter((REF_FRQ>0.01) & (REF_FRQ<0.99)) %>% filter(!(CHR==6 & POS>26e6 & POS<34e6))#hg37
  sumdat_hm3$PVAL = as.numeric(sumdat_hm3$PVAL)
  sumdat_ldsc = sumdat_hm3 %>% select(SNP_ID, REF, ALT, N, PVAL, BETA) %>% 
    rename(snpid=SNP_ID, a1=REF, a2=ALT, PVAL=PVAL, beta=BETA)
  
  bigreadr::fwrite2(sumdat_hm3 %>% select(SNP_ID, REF, ALT, BETA, SE, N), paste0(out_dir, '/hm3_sumdat.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
  sumdat_ldsc = sumdat_hm3 %>% select(SNP_ID, REF, ALT, N, PVAL, BETA) %>% 
    rename(snpid=SNP_ID, a1=REF, a2=ALT, PVAL=PVAL, beta=BETA)
  bigreadr::fwrite2(sumdat_ldsc, paste0(out_dir, '/ldsc-hm3.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
  ldsc_sumdat_all_hm3 = rbind(ldsc_sumdat_all_hm3, sumdat_ldsc)
  sumdat_PRScs = sumdat_hm3 %>% select(SNP_ID, REF, ALT, BETA, PVAL) %>% 
    rename(SNP=SNP_ID, A1=REF, A2=ALT, BETA=BETA, P=PVAL)
  bigreadr::fwrite2(sumdat_PRScs, paste0(out_dir, '/PRScs_sumdat.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
}
out_dir1 = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/')
bigreadr::fwrite2(ldsc_sumdat_all_hm3, paste0(out_dir1, '/ldsc-hm3.txt'), quote = FALSE, sep = "\t", row.names = FALSE)


############## 1kg small #############
for (chr in c(1:22)) {
  #alpha_list = c(0.5, 0.25, 0.125)
  alpha_list = c("auto")
  percent_list = c(0.2, 0.4, 0.6, 0.8)
  #percent_list = c(0)
  run_sh_list = c()
  for (alpha_i in 1:length(alpha_list)) {
    for (percent_i in 1:length(percent_list)) {
      alpha = alpha_list[alpha_i]
      percent = percent_list[percent_i]
      run_sh = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "/1kg/run_sh_prior/chr", chr, "/run_percent",percent,"_alpha",alpha,".sh")
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      sent_dir <- paste0("mkdir -p /dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "/1kg/result_projected_prior/chr", chr, "/percent",percent,"/alpha",alpha)
      cat('#!/bin/bash', "\n")
      output_dir = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "/1kg/result_projected_prior/chr", chr, "/percent",percent,"/alpha",alpha)
      ref = '/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_1kg/chr'
      sumdat = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/hm3_sumdat.txt')
      hfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/chr', chr, '_1kg_h2.txt')
      h2_txt = readLines(hfile)
      h2_txt = strsplit(h2_txt, ': ')[[1]][2]
      h2 = as.numeric(strsplit(h2_txt, split = " " )[[1]][1])
      h2_se = parse_number(strsplit(h2_txt, split = " " )[[1]][2])
      if (h2 < 0) {
        h2 = 0
      }
      sent_run <- paste0("python /dcs04/nilanjan/data/ydun/tools/PRSBridge/PRSBridge.py --percent ", percent,
                         " --chr ", chr, " --alpha ", alpha,  " --ref ", ref,  " --sumdat ", sumdat,  " --h2 ", h2,  " --h2_se ", h2_se
                         ,  " --method ", 'cg',  " --output ", output_dir)
      
      cat(sent_dir, "\n")
      cat("module load anaconda \n")
      cat(sent_run, "\n")
      sink()
    }
  }
  
  run_all = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "/1kg/run_sh_prior/chr",chr,"/run_all.sh")
  sink(run_all)
  mem = 2
  if(chr<21) mem = 5
  for (run_sh_i in 1:length(run_sh_list)) {
    run_sh = run_sh_list[run_sh_i]
    submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
    cat(submit_sh, "\n")
  }
  sink()
  system(paste0("cd /dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "/1kg/run_sh_prior/chr",chr))
  system(paste0("bash ", run_all))
}

############## UKBB large #############
for (chr in c(1:22)) {
  alpha_list = c(0.5, 0.25, 0.125)
  #alpha_list = c("auto")
  #percent_list = c(0)
  percent_list = c(0.2, 0.4, 0.6, 0.8)
  run_sh_list = c()
  for (alpha_i in 1:length(alpha_list)) {
    for (percent_i in 1:length(percent_list)) {
      alpha = alpha_list[alpha_i]
      percent = percent_list[percent_i]
      run_sh = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "_new/ukbb/run_sh_prior/chr", chr, "/run_percent",percent,"_alpha",alpha,".sh")
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sent_dir <- paste0("mkdir -p /dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "_new/ukbb/result_projected_prior/chr", chr, "/percent",percent,"/alpha",alpha)
      output_dir = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "_new/ukbb/result_projected_prior/chr", chr, "/percent",percent,"/alpha",alpha)
      ref = '/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/chr'
      sumdat = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/hm3_sumdat.txt')
      hfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/chr', chr, '_h2.txt')
      h2_txt = readLines(hfile)
      h2_txt = strsplit(h2_txt, ': ')[[1]][2]
      h2 = as.numeric(strsplit(h2_txt, split = " " )[[1]][1])
      h2_se = parse_number(strsplit(h2_txt, split = " " )[[1]][2])
      if (h2 < 0) {
        h2 = 0
      }
      sent_run <- paste0("python /dcs04/nilanjan/data/ydun/tools/PRSBridge/PRSBridge.py --percent ", percent,
                         " --chr ", chr, " --alpha ", alpha,  " --ref ", ref,  " --sumdat ", sumdat,  " --h2 ", h2,  " --h2_se ", h2_se
                         ,  " --method ", 'cg',  " --output ", output_dir)
      
      cat(sent_dir, "\n")
      cat("module load anaconda \n")
      cat(sent_run, "\n")
      sink()
    }
  }
  
  run_all = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "_new/ukbb/run_sh_prior/chr",chr,"/run_all.sh")
  sink(run_all)
  mem = 5
  if(chr<21) mem = 12
  for (run_sh_i in 1:length(run_sh_list)) {
    run_sh = run_sh_list[run_sh_i]
    submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
    cat(submit_sh, "\n")
  }
  sink()
  system(paste0("cd /dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "/ukbb/run_sh_prior/chr",chr))
  system(paste0("bash ", run_all))
}




############### generate run_sh file for PRScs_proj################
a_list = c("1.0")
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
percent_list = c(0.2, 0.4, 0.6, 0.8)
name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
ref = '1kg'
run_sh_list = c()
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj/result"))
for (chr in c(6:6)) {
  sumdat = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/PRScs_sumdat.txt')
  ldsc = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/ldsc-hm3.txt')
  sumdat_tmp = bigreadr::fread2(ldsc)
  for (phi_i in 1:length(phi_list)) {
    for (percent_i in percent_list) {
      for (a in a_list) {
        phi = phi_list[phi_i]
        name = name_list[phi_i]
        #if (file.exists(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj/result//percent", percent_i, "/chr", chr, "_pst_eff_a", a, "_b0.5_phi", name_list_prs[phi_i], "_chr", chr, ".txt"))) {next}
        run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_proj/run_sh/run_chr', chr, name, '_percent', percent_i, '_a', a, '.sh')
        run_sh_list = c(run_sh_list, run_sh)
        sink(run_sh)
        cat('#!/bin/bash', "\n")
        cat(paste0('mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/', ref, '/PRScs_proj/result/percent', percent_i), "\n")
        N = as.integer(floor(mean(sumdat_tmp$N)))
        output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_proj/result/percent', percent_i,'/chr', chr, '')
        
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --eigenval_rm=', percent_i)
        if(phi == 0){
          sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                             output_dir, ' --a=', a, ' --chrom=', chr, ' --eigenval_rm=', percent_i)
        }
        cat("module load anaconda \n")
        cat(sent_run, "\n")
        bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr)
        name1 = name_list_prs[phi_i]
        result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_proj/result/percent', percent_i, '/chr', chr,
                            '_pst_eff_a', a, '_b0.5_phi', name1, '_chr', chr)
        prs_output = paste0(result_dir, '_prs.txt')
        prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                        paste0('--score ', result_dir, '.txt'),
                        paste0(' 2 5 6 sum --bfile ',bfile),
                        paste0(' --out ', prs_output))
        #cat(prscode, "\n")
        sink()
      }
    }
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
  cat(submit_sh, "\n")
}
sink()
setwd(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj/run_sh/"))
system(paste0("bash ", run_all))


############### generate run_sh file for PRScs_proj_diff using 1kg ref################
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diff/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diff/result"))
a_list = c("1.0")
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
percent_list = c(0, 0.2, 0.4, 0.6, 0.8)
name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
for (chr in c(1:22)) {
  sumdat = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/PRScs_sumdat.txt')
  ldsc = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/ldsc-hm3.txt')
  sumdat_tmp = bigreadr::fread2(ldsc)
  for (phi_i in 1:length(phi_list)) {
    for (percent_i in percent_list) {
      for (a in a_list) {
        phi = phi_list[phi_i]
        name = name_list[phi_i]
        run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_proj_diff/run_sh/run_chr', chr, name, '_percent', percent_i, '_a', a, '.sh')
        run_sh_list = c(run_sh_list, run_sh)
        sink(run_sh)
        cat('#!/bin/bash', "\n")
        cat(paste0('mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/', ref, '/PRScs_proj_diff/result/percent', percent_i), "\n")
        N = as.integer(floor(mean(sumdat_tmp$N)))
        output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/', ref, '/PRScs_proj_diff/result/percent', percent_i,'/chr', chr, '')
        
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj_diff.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --eigenval_rm=', percent_i)
        if(phi == 0){
          sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj_diff.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                             output_dir, ' --a=', a, ' --chrom=', chr, ' --eigenval_rm=', percent_i)
        }
        cat("module load anaconda \n")
        cat(sent_run, "\n")
        bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr)
        name1 = name_list_prs[phi_i]
        result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_proj_diff/result/percent', percent_i, '/chr', chr,
                            '_pst_eff_a', a, '_b0.5_phi', name1, '_chr', chr)
        prs_output = paste0(result_dir, '_prs.txt')
        prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                        paste0('--score ', result_dir, '.txt'),
                        paste0(' 2 5 6 sum --bfile ',bfile),
                        " --threads 1",
                        paste0(' --out ', prs_output))
        #cat(prscode, "\n")
        sink()
      }
    }
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diff/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
  cat(submit_sh, "\n")
}
sink()
setwd(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diff/run_sh/"))
system(paste0("bash ", run_all))


############### generate run_sh file for PRScs_threshold using ukbb ref################
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_threshold/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_threshold/run_sh"))
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
threshold_list = c(0.01, 0.1, 1, 10, 100)
name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
for (chr in c(1:22)) {
  sumdat = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/PRScs_sumdat.txt')
  ldsc = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/', trait, '/data_Rcov/chr', chr, '/ldsc-hm3.txt')
  sumdat_tmp = bigreadr::fread2(ldsc)
  for (phi_i in 1:length(phi_list)) {
    for (threshold_i in threshold_list) {
      for (a in a_list) {
        phi = phi_list[phi_i]
        name = name_list[phi_i]
        run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_threshold/run_sh/run_chr', chr, name, '_threshold', threshold_i, '_a', a, '.sh')
        run_sh_list = c(run_sh_list, run_sh)
        sink(run_sh)
        cat('#!/bin/bash', "\n")
        cat(paste0('mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/', ref, '/PRScs_threshold/result/threshold', threshold_i), "\n")
        N = as.integer(floor(mean(sumdat_tmp$N)))
        output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_threshold/result/threshold', threshold_i,'/chr', chr, '')
        
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/non_convergence/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --threshold=', threshold_i)
        if(phi == 0){
          sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/non_convergence/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                             output_dir, ' --a=', a, ' --chrom=', chr, ' --threshold=', threshold_i)
        }
        cat("module load anaconda \n")
        cat(sent_run, "\n")
        bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr)
        name1 = name_list_prs[phi_i]
        result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_threshold/result/threshold', threshold_i, '/chr', chr,
                            '_pst_eff_a', a, '_b0.5_phi', name1, '_chr', chr)
        prs_output = paste0(result_dir, '_prs.txt')
        prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                        paste0('--score ', result_dir, '.txt'),
                        paste0(' 2 5 6 sum --bfile ',bfile),
                        paste0(' --out ', prs_output))
        #cat(prscode, "\n")
        sink()
      }
    }
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_threshold/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
  cat(submit_sh, "\n")
}
sink()
system(paste0("bash ", run_all))

