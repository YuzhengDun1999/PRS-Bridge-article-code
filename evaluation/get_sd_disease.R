############## PRS-Bridge ###########
library(RISCA)
cal_AUC = function(dat){
  roc_obj = roc.binary(status = "y", variable = "PRS", confounders = ~age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                       data = dat, precision = seq(0.05,0.95, by=0.01))
  return(2*(qnorm(roc_obj$auc)**2))
}
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
method = args[2]
ref = args[3]
#trait = 'BC'
#method = 'small'
#ref = 'ukbb'
system(paste0("mkdir -p /dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/"))
if (method == "small") { result_dir = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "/", ref, "/result_projected_prior/") }
if (method == "large") { result_dir = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "_large/", ref, "/result_projected_prior/") }
alpha_list = c(0.5, 0.25, 0.125)
percent_list = c(0, 0.2, 0.4, 0.6, 0.8)
for (chr in c(1:22)) {
  bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr)
  for (alpha_i in 1:length(alpha_list)) {
    for (percent_i in 1:length(percent_list)) {
      alpha = alpha_list[alpha_i]; percent = percent_list[percent_i]
      output_dir = paste0(result_dir, "chr", chr, "/percent",percent,"/alpha",alpha)
      prs_output = paste0(output_dir, '/prsresult_validation.txt')
      prsscore = paste0(output_dir, 'prsfile.txt')
      prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                      paste0('--score ', prsscore),
                      paste0(' 1 2 3 sum --bfile ',bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
    }
  }
}



library(dplyr)
cov = rbind(bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/tuning_cov.txt')),
            bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/validation_cov.txt')))
cov = cov %>% select(-FID)
names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
result = data.frame(matrix(ncol = 5,nrow = 0))
names(result) = c("seed", "alpha", "percent", "tuning_R2", "validation_R2"); result_i = 1
for (alpha in alpha_list) {
  for (percent in percent_list) {
    cov$PRS = 0
    for (chr in 1:22) {
      PRS = rbind(bigreadr::fread2(paste0(result_dir, "chr", chr, "/percent", percent, "/alpha", alpha, "/prsresult.txt.profile")),
                  bigreadr::fread2(paste0(result_dir, "chr", chr, "/percent", percent, "/alpha", alpha, "/prsresult_validation.txt.profile")))
      cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
    }
    for (seed in 1:10) {
      set.seed(seed)
      idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
      tuning = cov[idx, ]; validation = cov[-idx, ]
      result[result_i,] = c(seed, alpha, percent, cal_AUC(tuning), cal_AUC(validation)); result_i = result_i + 1
    }
  }
}
readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/Bridge_", method, '_', ref,'_all_result.txt'))

result_validation = data.frame(matrix(ncol = 4, nrow = 0))
names(result_validation) = c("seed", "alpha", "percent", "validation_R2"); result_validation_i = 1

for (seed_i in 1:10) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$alpha, result_tmp$percent, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
alpha_best = result_validation %>% group_by(alpha) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(alpha_best = alpha_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/Bridge_", method, '_', ref, '_validation_result_std.txt'))


############## PRS-Bridge-auto ###########
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
method_list = c("large", "small")
ref_list = c("ukbb", "1kg")
percent_list = c(0.2, 0.4, 0.6, 0.8)
alpha_list = c("auto")
for (j in 1:2) {
  method = method_list[j]; ref = ref_list[j]
  if (method == "small") { result_dir = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "/", ref, "/result_projected_prior/") }
  if (method == "large") { result_dir = paste0("/dcs04/nilanjan/data/ydun/PRS_Bridge/", trait, "_large/", ref, "/result_projected_prior/") }
  for (chr in c(1:22)) {
    for (alpha_i in 1:length(alpha_list)) {
      for (percent_i in 1:length(percent_list)) {
        bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/chr', chr)
        alpha = alpha_list[alpha_i]; percent = percent_list[percent_i]
        output_dir = paste0(result_dir, "chr", chr, "/percent", percent, "/alpha", alpha)
        prs_output = paste0(output_dir, '/prsresult.txt')
        prsscore = paste0(output_dir, '/coef.txt')
        prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                        paste0('--score ', prsscore),
                        paste0(' 1 3 4 sum --bfile ', bfile),
                        paste0(' --out ', prs_output))
        system(prscode)
        bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr)
        prs_output = paste0(output_dir, '/prsresult_validation.txt')
        prsscore = paste0(output_dir, '/coef.txt')
        prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                        paste0('--score ', prsscore),
                        paste0(' 1 3 4 sum --bfile ', bfile),
                        paste0(' --out ', prs_output))
        system(prscode)
      }
    }
  }
  
  library(dplyr)
  cov = rbind(bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/tuning_cov.txt')),
              bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/validation_cov.txt')))
  cov = cov %>% select(-FID)
  names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
  result = data.frame(matrix(ncol = 5,nrow = 0))
  names(result) = c("seed", "alpha", "percent", "tuning_R2", "validation_R2"); result_i = 1
  for (alpha in alpha_list) {
    for (percent in percent_list) {
      cov$PRS = 0
      for (chr in 1:22) {
        PRS = rbind(bigreadr::fread2(paste0(result_dir, "chr", chr, "/percent", percent, "/alpha", alpha, "/prsresult.txt.profile")),
                    bigreadr::fread2(paste0(result_dir, "chr", chr, "/percent", percent, "/alpha", alpha, "/prsresult_validation.txt.profile")))
        cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
      }
      for (seed in 1:10) {
        set.seed(seed)
        idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
        tuning = cov[idx, ]; validation = cov[-idx, ]
        result[result_i,] = c(seed, alpha, percent, cal_AUC(tuning), cal_AUC(validation)); result_i = result_i + 1
      }
    }
  }
  readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref,'_all_result.txt'))
  
  result = result %>% mutate(percent = as.numeric(percent), tuning_R2 = as.numeric(tuning_R2), validation_R2 = as.numeric(validation_R2))
  result_validation = data.frame(matrix(ncol = 3, nrow = 0))
  names(result_validation) = c("percent", "R2", "std"); result_validation_i = 1
  for (percent_i in percent_list) {
    result_tmp = result %>% filter(percent == percent_i)
    result_validation[result_validation_i,] = c(percent_i, mean(result_tmp$validation_R2), sqrt(var(result_tmp$validation_R2))); result_validation_i = result_validation_i + 1
  }
  readr::write_tsv(result_validation, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/Bridge_", method, "_auto", '_', ref, '_validation_result_std.txt'))
}



############## PRS-CS and PRS-CS-auto ###########
library(dplyr)
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
Ref = c('ukbb', '1kg')
#ref = '1kg'
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
name_list = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
method = "PRScs"
for (ref in Ref) {
  result_dir = paste0("/dcs04/nilanjan/data/ydun/PRScs/", trait, "/result_", ref, "/")
  ### for tuning calculate prs
  for (chr in c(1:22)) {
    for (phi_i in 1:length(phi_list)) {
      bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr)
      phi = phi_list[phi_i]; name = name_list[phi_i]
      prs_output = paste0(result_dir, '/chr', chr, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, '_validation_prs.txt')
      prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                      paste0('--score ', result_dir, '/chr', chr, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, '.txt'),
                      paste0(' 2 5 6 sum --bfile ', bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
      bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/chr', chr)
      prs_output = paste0(result_dir, '/chr', chr, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, '_prs.txt')
      prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                      paste0('--score ', result_dir, '/chr', chr, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, '.txt'),
                      paste0(' 2 5 6 sum --bfile ', bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
    }
  }
  
  library(dplyr)
  cov = rbind(bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/tuning_cov.txt')),
              bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/validation_cov.txt')))
  cov = cov %>% select(-FID)
  names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
  result = data.frame(matrix(ncol = 4, nrow = 0))
  names(result) = c("seed", "phi", "tuning_R2", "validation_R2"); result_i = 1
  for (phi_i in 1:length(phi_list)){  
    phi = phi_list[phi_i]; name = name_list[phi_i]
    cov$PRS = 0
    for (chr in 1:22) {
      PRS = rbind(bigreadr::fread2(paste0(result_dir, '/chr', chr, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, "_prs.txt.profile")),
                  bigreadr::fread2(paste0(result_dir, '/chr', chr, '_pst_eff_a1_b0.5_phi', name, '_chr', chr, "_validation_prs.txt.profile")))
      cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
    }
    for (seed in 1:10) {
      set.seed(seed)
      idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
      tuning = cov[idx, ]; validation = cov[-idx, ]
      result[result_i,] = c(seed, phi, cal_AUC(tuning), cal_AUC(validation)); result_i = result_i + 1
    }
  }
  readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref,'_all_result.txt'))
  result_auto = result %>% filter(phi == 0)
  readr::write_tsv(data.frame(phi_best = "auto", R2 = mean(result_auto$validation_R2), std = sqrt(var(result_auto$validation_R2))), 
                   paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_auto_', ref, '_validation_result_std.txt'))
  
  result_validation = data.frame(matrix(ncol = 3, nrow = 0))
  names(result_validation) = c("seed", "phi", "validation_R2"); result_validation_i = 1
  
  for (seed_i in 1:10) {
    result_tmp = result %>% filter(phi > 0) %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
    result_validation[result_validation_i,] = c(seed_i, result_tmp$phi, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
  }
  phi_best = result_validation %>% group_by(phi) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))
  
  readr::write_tsv(data.frame(phi_best = phi_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                   paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref, '_validation_result_std.txt'))
}


############## PRS-CS-proj and PRS-CS-auto-proj ###########
library(dplyr)
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
Ref = c('ukbb', '1kg')
#ref = '1kg'
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
name_list = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
percent_list = c(0.2, 0.4, 0.6, 0.8)
a_list = c("0.5", "1.0", "1.5")
method = "PRScs_proj"
for (ref in Ref) {
  for (a in a_list) {
    ### for tuning calculate prs
    for (chr in c(1:22)) {
      for (phi_i in 1:length(phi_list)) {
        for (percent in percent_list) {
          result_dir = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj/result/")
          bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr)
          phi = phi_list[phi_i]; name = name_list[phi_i]
          prs_output = paste0(result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '_validation_prs.txt')
          prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                          paste0('--score ', result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '.txt'),
                          paste0(' 2 5 6 sum --bfile ', bfile),
                          paste0(' --out ', prs_output))
          system(prscode)
          bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/chr', chr)
          prs_output = paste0(result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '_prs.txt')
          prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                          paste0('--score ', result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '.txt'),
                          paste0(' 2 5 6 sum --bfile ', bfile),
                          paste0(' --out ', prs_output))
          system(prscode)
        }
      }
    }
    
    library(dplyr)
    cov = rbind(bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/tuning_cov.txt')),
                bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/validation_cov.txt')))
    cov = cov %>% select(-FID)
    names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
    result = data.frame(matrix(ncol = 5, nrow = 0))
    names(result) = c("seed", "phi", "percent", "tuning_R2", "validation_R2"); result_i = 1
    for (phi_i in 1:length(phi_list)){ 
      for (percent in percent_list) {
        phi = phi_list[phi_i]; name = name_list[phi_i]
        cov$PRS = 0
        for (chr in 1:22) {
          PRS = rbind(bigreadr::fread2(paste0(result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, "_prs.txt.profile")),
                      bigreadr::fread2(paste0(result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, "_validation_prs.txt.profile")))
          cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
        }
        for (seed in 1:10) {
          set.seed(seed)
          idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
          tuning = cov[idx, ]; validation = cov[-idx, ]
          result[result_i,] = c(seed, phi, percent, cal_AUC(tuning), cal_AUC(validation)); result_i = result_i + 1
        }
      }
    }
    readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref, '_a', a, '_all_result.txt'))
    result_auto = result %>% filter(phi == 0)
    result_auto_out = data.frame(matrix(nrow = 0, ncol = 4)); result_auto_out_i = 1
    for (percent_i in percent_list) {
      result_auto_tmp = result_auto %>% filter(percent == percent_i)
      result_auto_out[result_auto_out_i,] = c("auto", percent_i, mean(result_auto_tmp$validation_R2), sqrt(var(result_auto_tmp$validation_R2)))
    }
    names(result_auto_out) = c("phi_best", "percent", "R2", "std")
    readr::write_tsv(result_auto_out, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_auto_', ref, '_a', a, '_validation_result_std.txt'))
    
    result_validation = data.frame(matrix(ncol = 3, nrow = 0))
    names(result_validation) = c("seed", "phi", "validation_R2"); result_validation_i = 1
    
    for (seed_i in 1:10) {
      result_tmp = result %>% filter(phi > 0) %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
      result_validation[result_validation_i,] = c(seed_i, result_tmp$phi, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
    }
    phi_best = result_validation %>% group_by(phi) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))
    
    readr::write_tsv(data.frame(phi_best = phi_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                     paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref, '_a', a, '_validation_result_std.txt'))
  }
}


############## PRS-CS-proj-diff and PRS-CS-auto-proj-diff ###########
library(dplyr)
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
Ref = c('ukbb', '1kg')
#ref = '1kg'
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
name_list = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
percent_list = c(0, 0.2, 0.4, 0.6, 0.8)
a_list = c("0.5", "1.0", "1.5")
method = "PRScs_proj_diff"
for (ref in Ref) {
  for (a in a_list) {
    ### for tuning calculate prs
    for (chr in c(1:22)) {
      for (phi_i in 1:length(phi_list)) {
        for (percent in percent_list) {
          result_dir = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_proj_diff/result/")
          bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr)
          phi = phi_list[phi_i]; name = name_list[phi_i]
          prs_output = paste0(result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '_validation_prs.txt')
          prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                          paste0('--score ', result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '.txt'),
                          paste0(' 2 5 6 sum --bfile ', bfile),
                          paste0(' --out ', prs_output))
          system(prscode)
          bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/chr', chr)
          prs_output = paste0(result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '_prs.txt')
          prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                          paste0('--score ', result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '.txt'),
                          paste0(' 2 5 6 sum --bfile ', bfile),
                          paste0(' --out ', prs_output))
          system(prscode)
        }
      }
    }
    
    library(dplyr)
    cov = rbind(bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/tuning_cov.txt')),
                bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/validation_cov.txt')))
    cov = cov %>% select(-FID)
    names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
    result = data.frame(matrix(ncol = 5, nrow = 0))
    names(result) = c("seed", "phi", "percent", "tuning_R2", "validation_R2"); result_i = 1
    for (phi_i in 1:length(phi_list)){ 
      for (percent in percent_list) {
        phi = phi_list[phi_i]; name = name_list[phi_i]
        cov$PRS = 0
        for (chr in 1:22) {
          PRS = rbind(bigreadr::fread2(paste0(result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, "_prs.txt.profile")),
                      bigreadr::fread2(paste0(result_dir, "/percent", percent, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, "_validation_prs.txt.profile")))
          cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
        }
        for (seed in 1:10) {
          set.seed(seed)
          idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
          tuning = cov[idx, ]; validation = cov[-idx, ]
          result[result_i,] = c(seed, phi, percent, cal_AUC(tuning), cal_AUC(validation)); result_i = result_i + 1
        }
      }
    }
    readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref, '_a', a, '_all_result.txt'))
    result_auto = result %>% filter(phi == 0)
    result_auto_out = data.frame(matrix(nrow = 0, ncol = 4)); result_auto_out_i = 1
    for (percent_i in percent_list) {
      result_auto_tmp = result_auto %>% filter(percent == percent_i)
      result_auto_out[result_auto_out_i,] = c("auto", percent_i, mean(result_auto_tmp$validation_R2), sqrt(var(result_auto_tmp$validation_R2)))
    }
    names(result_auto_out) = c("phi_best", "percent", "R2", "std")
    readr::write_tsv(result_auto_out, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_auto_', ref, '_a', a, '_validation_result_std.txt'))
    
    result_validation = data.frame(matrix(ncol = 3, nrow = 0))
    names(result_validation) = c("seed", "phi", "validation_R2"); result_validation_i = 1
    
    for (seed_i in 1:10) {
      result_tmp = result %>% filter(phi > 0) %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
      result_validation[result_validation_i,] = c(seed_i, result_tmp$phi, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
    }
    phi_best = result_validation %>% group_by(phi) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))
    
    readr::write_tsv(data.frame(phi_best = phi_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                     paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref, '_a', a, '_validation_result_std.txt'))
  }
}


############## PRS-CS-threshold ###########
library(dplyr)
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
Ref = c('ukbb', '1kg')
phi_list = c(0, 1, 0.01, 0.0001, 0.000001)
name_list = c('auto', '1e+00', '1e-02', '1e-04', '1e-06')
threshold_list = c("0.01", "0.1", "1", "10", "100")
a_list = c("1.0")
method = "PRScs_threshold"
for (ref in Ref) {
  for (a in a_list) {
    ### for tuning calculate prs
    for (chr in c(1:22)) {
      for (phi_i in 1:length(phi_list)) {
        for (threshold in threshold_list) {
          result_dir = paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/", ref, "/PRScs_threshold/result/")
          bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr)
          phi = phi_list[phi_i]; name = name_list[phi_i]
          prs_output = paste0(result_dir, "/threshold", threshold, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '_validation_prs.txt')
          prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                          paste0('--score ', result_dir, "/threshold", threshold, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '.txt'),
                          paste0(' 2 5 6 sum --bfile ', bfile),
                          paste0(' --out ', prs_output))
          system(prscode)
          bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/chr', chr)
          prs_output = paste0(result_dir, "/threshold", threshold, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '_prs.txt')
          prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                          paste0('--score ', result_dir, "/threshold", threshold, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '.txt'),
                          paste0(' 2 5 6 sum --bfile ', bfile),
                          paste0(' --out ', prs_output))
          system(prscode)
        }
      }
    }
    
    library(dplyr)
    cov = rbind(bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/tuning_cov.txt')),
                bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/validation_cov.txt')))
    cov = cov %>% select(-FID)
    names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
    result = data.frame(matrix(ncol = 5, nrow = 0))
    names(result) = c("seed", "phi", "threshold", "tuning_R2", "validation_R2"); result_i = 1
    for (phi_i in 1:length(phi_list)){ 
      for (threshold in threshold_list) {
        phi = phi_list[phi_i]; name = name_list[phi_i]
        cov$PRS = 0
        for (chr in 1:22) {
          PRS = rbind(bigreadr::fread2(paste0(result_dir, "/threshold", threshold, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, "_prs.txt.profile")),
                      bigreadr::fread2(paste0(result_dir, "/threshold", threshold, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, "_validation_prs.txt.profile")))
          cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
        }
        for (seed in 1:10) {
          set.seed(seed)
          idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
          tuning = cov[idx, ]; validation = cov[-idx, ]
          result[result_i,] = c(seed, phi, threshold, cal_AUC(tuning), cal_AUC(validation)); result_i = result_i + 1
        }
      }
    }
    readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref, '_a', a, '_all_result.txt'))
    result$threshold = as.numeric(result$threshold)
    result_auto = result %>% filter(phi == 0)
    result_auto_out = data.frame(matrix(nrow = 0, ncol = 4)); result_auto_out_i = 1
    for (threshold_i in threshold_list) {
      result_auto_tmp = result_auto %>% filter(threshold == as.numeric(threshold_i))
      result_auto_out[result_auto_out_i,] = c("auto", threshold_i, mean(result_auto_tmp$validation_R2), sqrt(var(result_auto_tmp$validation_R2))); result_auto_out_i = result_auto_out_i + 1
    }
    names(result_auto_out) = c("phi_best", "threshold", "R2", "std")
    readr::write_tsv(result_auto_out, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_auto_', ref, '_a', a, '_validation_result_std.txt'))
    
    phi_best_threshold = data.frame(matrix(ncol = 4, nrow = 0)); phi_best_threshold_i = 1
    names(phi_best_threshold) = c("phi_best", "threshold", "R2", "std")
    for (threshold_i in threshold_list) {
      result_validation = data.frame(matrix(ncol = 3, nrow = 0))
      names(result_validation) = c("seed", "phi_best", "validation_R2"); result_validation_i = 1
      for (seed_i in 1:10) {
        result_tmp = result %>% filter(phi > 0) %>% filter(threshold == as.numeric(threshold_i)) %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
        result_validation[result_validation_i,] = c(seed_i, result_tmp$phi, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
      }
      phi_best = result_validation %>% group_by(phi_best) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))
      phi_best_threshold[phi_best_threshold_i,] = c(phi_best[1,1], threshold_i, mean(result_validation$validation_R2), sqrt(var(result_validation$validation_R2))); phi_best_threshold_i = phi_best_threshold_i + 1
    }
    readr::write_tsv(phi_best_threshold, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref, '_a', a, '_validation_result_std.txt'))
      #readr::write_tsv(data.frame(phi_best = phi_best[1,1], threshold = threshold_i, R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
      #                 paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', ref, '_a', a, '_validation_result_std.txt'))
  }
}


############## LDpred2 ##############
library('dplyr')
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
N = c(args[2])
#N_list = c('N5000', '1kg', 'Nblk_provided')
#N_list = c('N5000')
method = 'ldpred2'
result_dir = paste0('/dcs04/nilanjan/data/ydun/ldpred2/', trait, "/")
### for tuning calculate prs
for (chr in c(1:22)) {
  bfile = paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr)
  for (x in 1:63) {
    prs_output = paste0(result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', N, '_validation_prs_', x,'.txt')
    if(file.exists(paste0(result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', N, '_x_', x, '.txt'))){
      prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink'),
                      paste0('--score ', result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', N, '_x_', x, '.txt'),
                      paste0(' 1 3 4',' sum --bfile ',bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
    }
  }
}

cov = rbind(bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/tuning_cov.txt')),
            bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/validation_cov.txt')))
cov = cov %>% select(-FID)
names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
result = data.frame(matrix(ncol = 4,nrow = 0))
names(result) = c("seed", "tune_i", "tuning_R2", "validation_R2"); result_i = 1
for (x in 1:63) {
  cov$PRS = 0
  for (chr in 1:22) {
    PRS = rbind(bigreadr::fread2(paste0(result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', N, '_prs_', x,".txt.profile")),
                bigreadr::fread2(paste0(result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', N, '_validation_prs_', x,".txt.profile")))
    cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:10) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, x, cal_AUC(tuning), cal_AUC(validation)); result_i = result_i + 1
  }
}
readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', N,'_all_result.txt'))

result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "tune_i", "validation_R2"); result_validation_i = 1

for (seed_i in 1:10) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$tune_i, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
tune_i_best = result_validation %>% group_by(tune_i) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(tune_i_best = tune_i_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', N, '_validation_result_std.txt'))



############## Lassosum ##############
library('dplyr')
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
N = c(args[2])
#N_list = c('N5000', 'N1kg')
#N_list = c('N5000')
method = 'Lassosum'
result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/' ,trait, '/LASSOSUM2/')

### for tuning calculate prs
for (chr in c(1:22)) {
  prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                  paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-EUR-ref_', N, '.txt'),
                  paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-123 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr),
                  paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/ref_", N, "validation"))
  #system(prscode)
  prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                  paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-EUR-ref_', N, '.txt'),
                  paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-123 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/chr', chr),
                  paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/ref_", N, "tuning"))
  #system(prscode)
}

cov = rbind(bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/tuning_cov.txt')),
            bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/validation_cov.txt')))
cov = cov %>% select(-FID)
names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
result = data.frame(matrix(ncol = 4, nrow = 0))
names(result) = c("seed", "tune_i", "tuning_R2", "validation_R2"); result_i = 1
PRS_list_tuning = lapply(1:22, function(chr) {
  file_path <- bigreadr::fread2(paste0(result_dir, '/chr', chr, "/ref_", N, "tuning.sscore"))
})
PRS_list_validation = lapply(1:22, function(chr) {
  file_path <- bigreadr::fread2(paste0(result_dir, '/chr', chr, "/ref_", N, "validation.sscore"))
})

for (x in 1:80) {
  cov$PRS = 0
  for (chr in 1:22) {
    PRS = rbind( PRS_list_tuning[[chr]], PRS_list_validation[[chr]])
    PRS_tmp = PRS[, c(2, x + 4)]; names(PRS_tmp) = c("IID", "SCORESUM")
    cov = merge(cov, PRS_tmp, by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:10) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, x, cal_AUC(tuning), cal_AUC(validation)); result_i = result_i + 1
  }
}
readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', N,'_all_result.txt'))

result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "tune_i", "validation_R2"); result_validation_i = 1

for (seed_i in 1:10) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$tune_i, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
tune_i_best = result_validation %>% group_by(tune_i) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(tune_i_best = tune_i_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', N, '_validation_result_std.txt'))




############## Lassosum2 ##############
library('dplyr')
args = commandArgs(trailingOnly = TRUE)
trait = args[1]
N = c(args[2])
#N_list = c('N5000', 'N1kg', 'Nblk_provided')
#N_list = c('N5000')
method = 'Lassosum2'
result_dir = paste0('/dcs04/nilanjan/data/ydun/SinglePRS/' ,trait, '/LASSOSUM2/')

### for tuning calculate prs
for (chr in c(1:22)) {
  prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                  paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-EUR-ref_', N, '.txt'),
                  paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-123 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/chr', chr),
                  paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/ref_", N, "validation"))
  system(prscode)
  prscode = paste(paste0('/dcs04/nilanjan/data/ydun/tools/plink2'),
                  paste0('--score ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, '/LASSOSUM2effect-hm3-EUR-ref_', N, '.txt'),
                  paste0(' 1 3 cols=+scoresums,-scoreavgs --score-col-nums 4-123 --bfile /dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/chr', chr),
                  paste0(' --out ', '/dcs04/nilanjan/data/ydun/SinglePRS/', trait, '/LASSOSUM2/chr', chr, "/ref_", N, "tuning"))
  system(prscode)
}

cov = rbind(bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/tuning_', trait, '/tuning_cov.txt')),
            bigreadr::fread2(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/validation/validation_', trait, '/validation_cov.txt')))
cov = cov %>% select(-FID)
names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
result = data.frame(matrix(ncol = 4,nrow = 0))
names(result) = c("seed", "tune_i", "tuning_R2", "validation_R2"); result_i = 1
PRS_list_tuning = lapply(1:22, function(chr) {
  file_path <- bigreadr::fread2(paste0(result_dir, '/chr', chr, "/ref_", N, "tuning.sscore"))
})
PRS_list_validation = lapply(1:22, function(chr) {
  file_path <- bigreadr::fread2(paste0(result_dir, '/chr', chr, "/ref_", N, "validation.sscore"))
})

for (x in 1:120) {
  cov$PRS = 0
  for (chr in 1:22) {
    PRS = rbind( PRS_list_tuning[[chr]], PRS_list_validation[[chr]])
    PRS_tmp = PRS[, c(2, x + 4)]; names(PRS_tmp) = c("IID", "SCORESUM")
    cov = merge(cov, PRS_tmp, by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:10) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, x, cal_AUC(tuning), cal_AUC(validation)); result_i = result_i + 1
  }
}
readr::write_tsv(result, paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', N,'_all_result.txt'))

result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "tune_i", "validation_R2"); result_validation_i = 1

for (seed_i in 1:10) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$tune_i, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
tune_i_best = result_validation %>% group_by(tune_i) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(tune_i_best = tune_i_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0("/dcs04/nilanjan/data/ydun/SinglePRS/", trait, "/result/", method, '_', N, '_validation_result_std.txt'))

