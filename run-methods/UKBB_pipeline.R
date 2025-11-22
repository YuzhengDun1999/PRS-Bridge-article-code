library(readr)
trait = "APOEB"
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/data/chr{1..22}"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_small/run_sh/chr{1..22}"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_small/run_sh/chr{1..22}"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_large/run_sh/chr{1..22}"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_large/run_sh/chr{1..22}"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/PRScs/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/PRScs/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/PRScs/result"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/PRScs/result"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ldpred2/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ldpred2/chr{1..22}"))

system(paste0("cp /dcs04/nilanjan/data/ydun/PRS_Bridge/BC_new/ukbb/run_sh_prior/run_all.sh /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_small/run_sh/"))
system(paste0("cp /dcs04/nilanjan/data/ydun/PRS_Bridge/BC_new/ukbb/run_sh_prior/run_all.sh /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_small/run_sh/"))
system(paste0("cp /dcs04/nilanjan/data/ydun/PRS_Bridge/BC_new/ukbb/run_sh_prior/run_all.sh /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_large/run_sh/"))
system(paste0("cp /dcs04/nilanjan/data/ydun/PRS_Bridge/BC_new/ukbb/run_sh_prior/run_all.sh /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_large/run_sh/"))

################### preprocess summary level data of UKBB ################
library(dplyr)
sumdat = bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/GWAS_ukbb/', trait, '/sumdat_Rcov.txt'))
ldsc_sumdat_all_hm3 = data.frame()
ldsc_sumdat_all_mega = data.frame()
for (chr in 1:22) {
  ref_dir = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/chr', chr, '/')
  out_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/')
  sumdat_mega = sumdat[sumdat$CHR == chr, ]
  sumdat_mega = sumdat_mega %>% filter(((BETA/SE)^2/N <= 80)) %>%
    select(CHR, POS, SNP_ID, REF, ALT, REF_FRQ, PVAL, BETA, SE, N) %>%
    filter((REF_FRQ>0.01) & (REF_FRQ<0.99)) %>% filter(!(CHR==6 & POS>26e6 & POS<34e6))#hg37
  sumdat_mega$PVAL = as.numeric(sumdat_mega$PVAL)
  ref = bigreadr::fread2(paste0(ref_dir, 'chr_', chr, '_merged.bim'))
  sumdat_ldsc = sumdat_mega %>% select(SNP_ID, REF, ALT, N, PVAL, BETA) %>% 
    rename(snpid=SNP_ID, a1=REF, a2=ALT, PVAL=PVAL, beta=BETA)
  ldsc_sumdat_all_mega = rbind(ldsc_sumdat_all_mega, sumdat_ldsc)
  
  sumdat_hm3 = sumdat_mega %>% filter(SNP_ID %in% ref[,2])
  
  bigreadr::fwrite2(sumdat_hm3 %>% select(SNP_ID, REF, ALT, BETA, SE, N), paste0(out_dir, '/hm3_sumdat.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
  sumdat_ldsc = sumdat_hm3 %>% select(SNP_ID, REF, ALT, N, PVAL, BETA) %>% 
    rename(snpid=SNP_ID, a1=REF, a2=ALT, PVAL=PVAL, beta=BETA)
  bigreadr::fwrite2(sumdat_ldsc, paste0(out_dir, '/ldsc-hm3.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
  ldsc_sumdat_all_hm3 = rbind(ldsc_sumdat_all_hm3, sumdat_ldsc)
  sumdat_PRScs = sumdat_hm3 %>% select(SNP_ID, REF, ALT, BETA, PVAL) %>% 
    rename(SNP=SNP_ID, A1=REF, A2=ALT, BETA=BETA, P=PVAL)
  bigreadr::fwrite2(sumdat_PRScs, paste0(out_dir, '/PRScs_sumdat.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
}
out_dir1 = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/')
bigreadr::fwrite2(ldsc_sumdat_all_hm3, paste0(out_dir1, '/ldsc-hm3.txt'), quote = FALSE, sep = "\t", row.names = FALSE)


################## ldsc using UKBB ref ################
sink(paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/get_h2.sh'))
cat('module load conda_R',"\n")
cat('source activate ldsc',"\n")
for(chr in 1:22){
  sumstat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc-hm3.txt')
  out = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc_sumdat_1')
  sent_run = paste0('~/ldsc/munge_sumstats.py --sumstats ', sumstat, ' --out ', out)
  cat(sent_run, "\n")
  ref = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/chr', chr, '/chr_', chr, '_merged')
  out1 = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/chr', chr, '_h2')
  sent_run1 = paste0('~/ldsc/ldsc.py --h2 ', out, '.sumstats.gz --ref-ld ', ref, ' --w-ld ', ref, ' --out ', out1)
  sent_run2 = paste0('cat ', out1, '.log',  ' | grep "Observed scale h2:" > ', out1, '.txt')
  cat(sent_run1, "\n")
  cat(sent_run2, "\n")
}
sink()
system(paste0("bash ", '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/get_h2.sh'))


################ ldsc using 1kg ref################
sink(paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/get_h2_1kg.sh'))
cat('module load conda_R',"\n")
cat('source activate ldsc',"\n")
for(chr in 1:22){
  sumstat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc-hm3.txt')
  out = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc_1kg_sumdat_1')
  sent_run = paste0('~/ldsc/munge_sumstats.py --sumstats ', sumstat, ' --out ', out)
  cat(sent_run, "\n")
  ref = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_1kg/chr', chr, '/chr_', chr, '_merged')
  out1 = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/chr', chr, '_1kg_h2')
  sent_run1 = paste0('~/ldsc/ldsc.py --h2 ', out, '.sumstats.gz --ref-ld ', ref, ' --w-ld ', ref, ' --out ', out1)
  sent_run2 = paste0('cat ', out1, '.log',  ' | grep "Observed scale h2:" > ', out1, '.txt')
  cat(sent_run1, "\n")
  cat(sent_run2, "\n")
}
sink()
system(paste0("bash ", '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/get_h2_1kg.sh'))



library(readr)
####### generate run_sh file for PRSBridge using 1kg small ref #######
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
      run_sh = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_small/run_sh/chr", chr, "/run_percent",percent,"_alpha",alpha,".sh")
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sent_dir <- paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_small/result/chr", chr, "/percent",percent,"/alpha",alpha)
      output_dir = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_small/result/chr", chr, "/percent",percent,"/alpha",alpha)
      ref = '/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_1kg/chr'
      sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/hm3_sumdat.txt')
      hfile = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/chr', chr, '_1kg_h2.txt')
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
  
  run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_small/run_sh/chr",chr,"/run_all.sh")
  sink(run_all)
  mem = 2
  if(chr<21) mem = 5
  for (run_sh_i in 1:length(run_sh_list)) {
    run_sh = run_sh_list[run_sh_i]
    submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
    cat(submit_sh, "\n")
  }
  sink()
  system(paste0("cd /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_small/run_sh/chr",chr))
  system(paste0("bash ", run_all))
}


############### generate run_sh file for PRSBridge using 1kg large################
for (chr in c(1:22)) {
  #alpha_list = c(0.5, 0.25, 0.125)
  alpha_list = c("auto")
  #percent_list = c(0)
  percent_list = c(0.2, 0.4, 0.6, 0.8)
  run_sh_list = c()
  for (alpha_i in 1:length(alpha_list)) {
    for (percent_i in 1:length(percent_list)) {
      alpha = alpha_list[alpha_i]
      percent = percent_list[percent_i]
      run_sh = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_large/run_sh/chr", chr, "/run_percent",percent,"_alpha",alpha,".sh")
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sent_dir <- paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_large/result/chr", chr, "/percent",percent,"/alpha",alpha)
      output_dir = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_large/result/chr", chr, "/percent",percent,"/alpha",alpha)
      ref = '/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_1kg_new/chr'
      sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/hm3_sumdat.txt')
      hfile = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/chr', chr, '_1kg_h2.txt')
      h2_txt = readLines(hfile)
      h2_txt = strsplit(h2_txt, ': ')[[1]][2]
      h2 = as.numeric(strsplit(h2_txt, split = " " )[[1]][1])
      h2_se = parse_number(strsplit(h2_txt, split = " " )[[1]][2])
      if (h2 < 0) {
        h2 = 0
      }
      sent_run <- paste0("python /dcs04/nilanjan/data/ydun/PRS_Bridge/PRS_Bayes/PRSBridge.py --percent ", percent,
                         " --chr ", chr, " --alpha ", alpha,  " --ref ", ref,  " --sumdat ", sumdat,  " --h2 ", h2,  " --h2_se ", h2_se
                         ,  " --method ", 'cg',  " --output ", output_dir)
      
      cat(sent_dir, "\n")
      cat("module load anaconda \n")
      cat(sent_run, "\n")
      sink()
    }
  }
  
  run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_large/run_sh/chr",chr,"/run_all.sh")
  sink(run_all)
  mem = 5
  if(chr<21) mem = 12
  for (run_sh_i in 1:length(run_sh_list)) {
    run_sh = run_sh_list[run_sh_i]
    submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
    cat(submit_sh, "\n")
  }
  sink()
  system(paste0("cd /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/Bridge_large/run_sh/chr",chr))
  system(paste0("bash ", run_all))
}


####### generate run_sh file for PRSBridge using ukbb small ref #######
for (chr in c(1:22)) {
  #alpha_list = c(0.5, 0.25, 0.125)
  alpha_list = c("auto")
  percent_list = c(0)
  #percent_list = c(0.2, 0.4, 0.6, 0.8)
  run_sh_list = c()
  for (alpha_i in 1:length(alpha_list)) {
    for (percent_i in 1:length(percent_list)) {
      alpha = alpha_list[alpha_i]
      percent = percent_list[percent_i]
      run_sh = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_small/run_sh/chr", chr, "/run_percent",percent,"_alpha",alpha,".sh")
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sent_dir <- paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_small/result/chr", chr, "/percent",percent,"/alpha",alpha)
      output_dir = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_small/result/chr", chr, "/percent",percent,"/alpha",alpha)
      ref = '/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/chr'
      sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/hm3_sumdat.txt')
      hfile = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/chr', chr, '_h2.txt')
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
  
  run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_small/run_sh/chr",chr,"/run_all.sh")
  sink(run_all)
  mem = 2
  if(chr<21) mem = 5
  for (run_sh_i in 1:length(run_sh_list)) {
    run_sh = run_sh_list[run_sh_i]
    submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
    cat(submit_sh, "\n")
  }
  sink()
  system(paste0("cd /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_small/run_sh/chr",chr))
  system(paste0("bash ", run_all))
}


############### generate run_sh file for PRSBridge using ukbb large################
for (chr in c(1:22)) {
  #alpha_list = c(0.5, 0.25, 0.125)
  alpha_list = c("auto")
  #percent_list = c(0)
  percent_list = c(0.2, 0.4, 0.6, 0.8)
  run_sh_list = c()
  for (alpha_i in 1:length(alpha_list)) {
    for (percent_i in 1:length(percent_list)) {
      alpha = alpha_list[alpha_i]
      percent = percent_list[percent_i]
      run_sh = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_large/run_sh/chr", chr, "/run_percent",percent,"_alpha",alpha,".sh")
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sent_dir <- paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_large/result/chr", chr, "/percent",percent,"/alpha",alpha)
      output_dir = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_large/result/chr", chr, "/percent",percent,"/alpha",alpha)
      ref = '/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb_new/chr'
      sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/hm3_sumdat.txt')
      hfile = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/chr', chr, '_h2.txt')
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
  
  run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_large/run_sh/chr",chr,"/run_all.sh")
  sink(run_all)
  mem = 5
  if(chr<21) mem = 12
  for (run_sh_i in 1:length(run_sh_list)) {
    run_sh = run_sh_list[run_sh_i]
    submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
    cat(submit_sh, "\n")
  }
  sink()
  system(paste0("cd /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/Bridge_large/run_sh/chr",chr))
  system(paste0("bash ", run_all))
}



############### generate run_sh file for PRScs using ukbb ref################
a_list = c(0.5, 1.5)
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
for (chr in c(1:22)) {
  for (phi_i in 1:length(phi_list)) {
    for (a in a_list) {
      phi = phi_list[phi_i]
      name = name_list[phi_i]
      run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/ukbb/PRScs/run_sh/run_chr', chr, name, '_a', a,'.sh')
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/PRScs_sumdat.txt')
      ldsc = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc-hm3.txt')
      sumdat_tmp = bigreadr::fread2(ldsc)
      N = as.integer(floor(mean(sumdat_tmp$N)))
      output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/ukbb/PRScs/result/chr', chr)
      
      sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/ukbb/ldblk_ukbb_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                         output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi)
      if(phi == 0){
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/ukbb/ldblk_ukbb_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr)
      }
      cat("module load anaconda \n")
      cat(sent_run, "\n")
      bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr)
      name1 = name_list_prs[phi_i]
      result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/ukbb/PRScs/result/chr', chr,
                          '_pst_eff_a', a, '_b0.5_phi', name1, '_chr', chr)
      prs_output = paste0(result_dir, '_prs.txt')
      prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                      paste0('--score ', result_dir, '.txt'),
                      paste0(' 2 5 6 sum --bfile ', bfile),
                      " --threads 1",
                      paste0(' --out ', prs_output))
      #cat(prscode, "\n")
      sink()
    }
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/PRScs/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch  --mem=', mem,'G ', run_sh)
  cat(submit_sh, "\n")
}
sink()
setwd(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/PRScs/run_sh/"))
system(paste0("bash ", run_all))


############### generate run_sh file for PRScs using 1kg################
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
for (chr in c(1:22)) {
  for (phi_i in 1:length(phi_list)) {
    for (a in a_list) {
      phi = phi_list[phi_i]
      name = name_list[phi_i]
      run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/1kg/PRScs/run_sh/run_chr', chr, name, '_a', a,'.sh')
      run_sh_list = c(run_sh_list, run_sh)
      sink(run_sh)
      cat('#!/bin/bash', "\n")
      sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/PRScs_sumdat.txt')
      ldsc = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc-hm3.txt')
      sumdat_tmp = bigreadr::fread2(ldsc)
      N = as.integer(floor(mean(sumdat_tmp$N)))
      output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/1kg/PRScs/result/chr', chr)
      
      sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/1kg/ldblk_1kg_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                         output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi)
      if(phi == 0){
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/1kg/ldblk_1kg_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr)
      }
      cat("module load anaconda \n")
      cat(sent_run, "\n")
      bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr)
      name1 = name_list_prs[phi_i]
      result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/1kg/PRScs/result/chr', chr,
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

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/PRScs/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch  --mem=', mem,'G ', run_sh)
  cat(submit_sh, "\n")
}
sink()
setwd(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/PRScs/run_sh/"))
system(paste0("bash ", run_all))


############### generate run_sh file for ldpred2 using 1kg################
race = 1
run_sh_list = c()
for (ref_N in c('1kg')) {
  for (chr in 1:22) {
    run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/ldpred2/run_sh/chr', chr,'_ref_N', ref_N,'.sh')
    run_sh_list = c(run_sh_list, run_sh)
    sink(run_sh)
    sent_run = paste0('Rscript /dcs04/nilanjan/data/ydun/SinglePRS/code/ldpred2_trait_1kg.R ', race, ' ', chr, ' ', ref_N, ' ', trait)
    cat("module load conda_R \n")
    cat(sent_run, "\n")
    sink()
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ldpred2/run_sh/run_all_1kg.sh")
sink(run_all)
for (chr in 1:22) {
  mem = 5
  if(chr>13) mem = 2.5
  for (ref_N in c('1kg')) {
    run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/ldpred2/run_sh/chr', chr,'_ref_N', ref_N,'.sh')
    submit_sh = paste0('qsub -cwd -pe local 4 -l mem_free=', mem,'G,h_vmem=', mem,"G,h='!(compute-053|compute-054|compute-066)' ", run_sh)
    cat(submit_sh, "\n")
  }
}
sink()
system(paste0("bash ", run_all))


############### generate run_sh file for ldpred2 using ukbb################
race = 1
run_sh_list = c()
for (ref_N in c(5000)) {
  for (chr in 1:22) {
    run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/ldpred2/run_sh/chr', chr,'_ref_N', ref_N,'.sh')
    run_sh_list = c(run_sh_list, run_sh)
    sink(run_sh)
    sent_run = paste0('Rscript /dcs04/nilanjan/data/ydun/SinglePRS/code/ldpred2_trait.R ', race, ' ', chr, ' ', ref_N, ' ', trait)
    cat("module load conda_R \n")
    cat(sent_run, "\n")
    sink()
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ldpred2/run_sh/run_all_ukbb.sh")
sink(run_all)
for (chr in 1:22) {
  mem = 5
  if(chr>13) mem = 2.5
  for (ref_N in c(5000)) {
    run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/ldpred2/run_sh/chr', chr,'_ref_N', ref_N,'.sh')
    submit_sh = paste0('qsub -cwd -pe local 4 -l mem_free=', mem,'G,h_vmem=', mem,"G,h='!(compute-053|compute-054|compute-066)' ", run_sh)
    cat(submit_sh, "\n")
  }
}
sink()
system(paste0("cd /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ldpred2/run_sh/"))
system(paste0("bash ", run_all))


############### generate run_sh file for ldpred2 using blked################
race = 1
run_sh_list = c()
for (ref_N in c('blk_provided')) {
  for (chr in 1:22) {
    run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/ldpred2/run_sh/chr', chr,'_ref_N', ref_N,'.sh')
    run_sh_list = c(run_sh_list, run_sh)
    sink(run_sh)
    sent_run = paste0('Rscript /dcs04/nilanjan/data/ydun/SinglePRS/code/ldpred2_trait_blk.R ', race, ' ', chr, ' ', ref_N, ' ', trait)
    cat("module load conda_R \n")
    cat(sent_run, "\n")
    sink()
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ldpred2/run_sh/run_all_blk.sh")
sink(run_all)
for (chr in 1:22) {
  mem = 5
  if(chr>13) mem = 2.5
  for (ref_N in c('blk_provided')) {
    run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/ldpred2/run_sh/chr', chr,'_ref_N', ref_N,'.sh')
    submit_sh = paste0('qsub -cwd -pe local 4 -l mem_free=', mem,'G,h_vmem=', mem,"G,h='!(compute-053|compute-054|compute-066)' ", run_sh)
    cat(submit_sh, "\n")
  }
}
sink()
system(paste0("cd /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ldpred2/run_sh/"))
system(paste0("bash ", run_all))



############### generate run_sh file for PRScs_proj################
a_list = c("1")
phi_list = c(1, 0.01, 0.0001, 0.000001)
percent_list = c(0.2, 0.4, 0.6, 0.8)
name_list = c('_1', '_e2', '_e4', '_e6')
name_list_prs = c('1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
ref = '1kg'
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj/result"))
for (chr in c(1:22)) {
  for (phi_i in 1:length(phi_list)) {
    for (percent_i in percent_list) {
      for (a in a_list) {
        phi = phi_list[phi_i]
        name = name_list[phi_i]
        run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_proj/run_sh/run_chr', chr, name, '_percent', percent_i, '_a', a, '.sh')
        run_sh_list = c(run_sh_list, run_sh)
        sink(run_sh)
        cat('#!/bin/bash', "\n")
        cat(paste0('mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/', ref, '/PRScs_proj/result/percent', percent_i), "\n")
        sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/PRScs_sumdat.txt')
        ldsc = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc-hm3.txt')
        sumdat_tmp = bigreadr::fread2(ldsc)
        N = as.integer(floor(mean(sumdat_tmp$N)))
        output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/', ref, '/PRScs_proj/result/percent', percent_i,'/chr', chr, '')
        
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --eigenval_rm=', percent_i)
        if(phi == 0){
          sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
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
                        " --threads 1",
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
  submit_sh = paste0('sbatch  --mem=', mem,'G ', run_sh)
  cat(submit_sh, "\n")
}
sink()
setwd(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj/run_sh/"))
system("rm *.out")
system(paste0("bash ", run_all))

############### generate run_sh file for PRScs_proj_diff################
ref = "1kg"
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diff/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diff/result"))
#a_list = c("0.5", "1.0", "1.5")
a_list = c("1.0")
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
percent_list = c(0, 0.2, 0.4, 0.6, 0.8)
name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
for (chr in c(1:22)) {
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
        sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/PRScs_sumdat.txt')
        ldsc = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc-hm3.txt')
        sumdat_tmp = bigreadr::fread2(ldsc)
        N = as.integer(floor(mean(sumdat_tmp$N)))
        output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/', ref, '/PRScs_proj_diff/result/percent', percent_i,'/chr', chr, '')
        
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj_diff.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --eigenval_rm=', percent_i)
        if(phi == 0){
          sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj_diff.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur  --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
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
                        paste0(' 2 5 6 sum --bfile ', bfile),
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


############### generate run_sh file for PRScs_proj_diag################
ref = "1kg"
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diag/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diag/result"))
#a_list = c("0.5", "1.0", "1.5")
a_list = c("1.0")
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
#percent_list = c("N", "0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1")
percent_list = c("1")
name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
for (chr in c(1:22)) {
  for (phi_i in 1:length(phi_list)) {
    for (percent_i in percent_list) {
      for (a in a_list) {
        phi = phi_list[phi_i]
        name = name_list[phi_i]
        run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_proj_diag/run_sh/run_chr', chr, name, '_percent', percent_i, '_a', a, '.sh')
        run_sh_list = c(run_sh_list, run_sh)
        sink(run_sh)
        cat('#!/bin/bash', "\n")
        cat(paste0('mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/', ref, '/PRScs_proj_diag/result/percent', percent_i), "\n")
        sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/PRScs_sumdat.txt')
        ldsc = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc-hm3.txt')
        sumdat_tmp = bigreadr::fread2(ldsc)
        N = as.integer(floor(mean(sumdat_tmp$N)))
        output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/', ref, '/PRScs_proj_diag/result/percent', percent_i,'/chr', chr, '')
        
        percent_i_input = percent_i
        if (percent_i == "N") percent_i_input = 1/sqrt(N)
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj_diag.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --eigenval_rm=', percent_i_input)
        if(phi == 0){
          sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/PRScs_proj_diag.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/', ref, '/ldblk_', ref, '_eur  --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                             output_dir, ' --a=', a, ' --chrom=', chr, ' --eigenval_rm=', percent_i_input)
        }
        cat("module load anaconda \n")
        cat(sent_run, "\n")
        bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr)
        name1 = name_list_prs[phi_i]
        result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/', ref, '/PRScs_proj_diag/result/percent', percent_i, '/chr', chr,
                            '_pst_eff_a', a, '_b0.5_phi', name1, '_chr', chr)
        prs_output = paste0(result_dir, '_prs.txt')
        prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                        paste0('--score ', result_dir, '.txt'),
                        paste0(' 2 5 6 sum --bfile ', bfile),
                        " --threads 1",
                        paste0(' --out ', prs_output))
        #cat(prscode, "\n")
        sink()
      }
    }
  }
}

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diag/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch --mem=', mem,'G --exclude=compute-053,compute-054,compute-057,compute-058,compute-069,compute-089,compute-096 ', run_sh)
  cat(submit_sh, "\n")
}
sink()
setwd(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diag/run_sh/"))
system(paste0("bash ", run_all))

############### generate run_sh file for PRScs_threshold using ukbb ref################
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/PRScs_threshold/run_sh"))
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/PRScs_threshold/run_sh"))
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
threshold_list = c(0.01, 100)
name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
for (chr in c(1:22)) {
  for (phi_i in 1:length(phi_list)) {
    for (threshold_i in threshold_list) {
      for (a in c("1")) {
        phi = phi_list[phi_i]
        name = name_list[phi_i]
        run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/ukbb/PRScs_threshold/run_sh/run_chr', chr, name, '_threshold', threshold_i, '_a', a, '.sh')
        run_sh_list = c(run_sh_list, run_sh)
        sink(run_sh)
        cat('#!/bin/bash', "\n")
        cat(paste0('mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/ukbb/PRScs_threshold/result/threshold', threshold_i), "\n")
        sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/PRScs_sumdat.txt')
        ldsc = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc-hm3.txt')
        sumdat_tmp = bigreadr::fread2(ldsc)
        N = as.integer(floor(mean(sumdat_tmp$N)))
        output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/ukbb/PRScs_threshold/result/threshold', threshold_i,'/chr', chr, '')
        
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/non_convergence/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/ukbb/ldblk_ukbb_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --threshold=', threshold_i)
        if(phi == 0){
          sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/non_convergence/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/ukbb/ldblk_ukbb_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                             output_dir, ' --a=', a, ' --chrom=', chr, ' --threshold=', threshold_i)
        }
        cat("module load anaconda \n")
        cat(sent_run, "\n")
        bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr)
        name1 = name_list_prs[phi_i]
        result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/ukbb/PRScs_threshold/result/threshold', threshold_i, '/chr', chr,
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

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/PRScs_threshold/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch  --mem=', mem,'G ', run_sh)
  cat(submit_sh, "\n")
}
sink()
setwd(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/ukbb/PRScs_threshold/run_sh/"))
system(paste0("bash ", run_all))

############### generate run_sh file for PRScs_threshold using 1kg ref################
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
threshold_list = c(0.01, 100)
name_list = c('_auto', '_1', '_e2', '_e4', '_e6')
name_list_prs = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
run_sh_list = c()
for (chr in c(1:22)) {
  for (phi_i in 1:length(phi_list)) {
    for (threshold_i in threshold_list) {
      for (a in c("1")) {
        phi = phi_list[phi_i]
        name = name_list[phi_i]
        run_sh = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/1kg/PRScs_threshold/run_sh/run_chr', chr, name, '_threshold', threshold_i, '_a', a, '.sh')
        run_sh_list = c(run_sh_list, run_sh)
        sink(run_sh)
        cat('#!/bin/bash', "\n")
        cat(paste0('mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/1kg/PRScs_threshold/result/threshold', threshold_i), "\n")
        sumdat = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/PRScs_sumdat.txt')
        ldsc = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/data/chr', chr, '/ldsc-hm3.txt')
        sumdat_tmp = bigreadr::fread2(ldsc)
        N = as.integer(floor(mean(sumdat_tmp$N)))
        output_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait,'/1kg/PRScs_threshold/result/threshold', threshold_i,'/chr', chr, '')
        
        sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/non_convergence/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/1kg/ldblk_1kg_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N, ' --out_dir=',
                           output_dir, ' --a=', a, ' --chrom=', chr, ' --phi=', phi, ' --threshold=', threshold_i)
        if(phi == 0){
          sent_run <- paste0('python /dcs04/nilanjan/data/ydun/PRScs/non_convergence/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/ydun/PRScs/ref_data/1kg/ldblk_1kg_eur --bim_prefix=/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/chr', chr,' --sst_file=', sumdat, ' --n_gwas=', N ,' --out_dir=',
                             output_dir, ' --a=', a, ' --chrom=', chr, ' --threshold=', threshold_i)
        }
        cat("module load anaconda \n")
        cat(sent_run, "\n")
        bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning/chr', chr)
        name1 = name_list_prs[phi_i]
        result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/1kg/PRScs_threshold/result/threshold', threshold_i, '/chr', chr,
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

run_all = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/PRScs_threshold/run_sh/run_all.sh")
sink(run_all)
mem = 3
for (run_sh_i in 1:length(run_sh_list)) {
  run_sh = run_sh_list[run_sh_i]
  submit_sh = paste0('sbatch  --mem=', mem,'G ', run_sh)
  cat(submit_sh, "\n")
}
sink()
setwd(paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/1kg/PRScs_threshold/run_sh/"))
system(paste0("bash ", run_all))


