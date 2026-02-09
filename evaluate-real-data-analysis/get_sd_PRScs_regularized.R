####### This script is used for evaluation of PRS-CS-Regularized ###########
outcome = "continuous" # change to continuous if your outcome is continuous 
source("evaluation/evaluation_metrics.R")

######### Change the parameters here if you use different set of parameters in run-methods/PRScs.R
a_list = c(1)
phi_list = c(1, 0.01, 0.0001, 0.000001)
name_list = c('1e+00', '1e-02', '1e-04', '1e-06')
regularized_list = c("1")
Ref = c('ukbb') # LD reference data you use in PRScs.R, can be 'ukbb' and '1kg'
trait = "Example"

library(dplyr)
method = "PRScs_proj"
for (ref in Ref) {
  for (a in a_list) {
    ### for tuning calculate prs
    for (chr in c(1:22)) {
      for (phi_i in 1:length(phi_list)) {
        for (regularized in regularized_list) {
          result_dir = paste0(trait, "/", ref, "/", method, "/result/")
          bfile = paste0('validation/validation/chr', chr)
          phi = phi_list[phi_i]; name = name_list[phi_i]
          prs_output = paste0(result_dir, "/regularized", regularized, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '_validation_prs.txt')
          prscode = paste(paste0('plink'),
                          paste0('--score ', result_dir, "/regularized", regularized, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '.txt'),
                          paste0(' 2 5 6 sum --bfile ', bfile),
                          paste0(' --out ', prs_output))
          system(prscode)
          bfile = paste0('validation/tuning/chr', chr)
          prs_output = paste0(result_dir, "/regularized", regularized, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '_prs.txt')
          prscode = paste(paste0('plink'),
                          paste0('--score ', result_dir, "/regularized", regularized, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, '.txt'),
                          paste0(' 2 5 6 sum --bfile ', bfile),
                          paste0(' --out ', prs_output))
          system(prscode)
        }
      }
    }
    
    library(dplyr)
    cov = rbind(bigreadr::fread2(paste0('validation/tuning/', trait, '_cov.txt')),
                bigreadr::fread2(paste0('validation/validation/', trait, '_cov.txt')))
    cov = cov %>% select(-FID)
    names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
    result = data.frame(matrix(ncol = 5, nrow = 0))
    names(result) = c("seed", "phi", "regularized", "tuning_R2", "validation_R2"); result_i = 1
    for (phi_i in 1:length(phi_list)){ 
      for (regularized in regularized_list) {
        phi = phi_list[phi_i]; name = name_list[phi_i]
        cov$PRS = 0
        for (chr in 1:22) {
          PRS = rbind(bigreadr::fread2(paste0(result_dir, "/regularized", regularized, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, "_prs.txt.profile")),
                      bigreadr::fread2(paste0(result_dir, "/regularized", regularized, '/chr', chr, '_pst_eff_a', a, '_b0.5_phi', name, '_chr', chr, "_validation_prs.txt.profile")))
          cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
        }
        for (seed in 1:100) {
          set.seed(seed)
          idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
          tuning = cov[idx, ]; validation = cov[-idx, ]
          result[result_i,] = c(seed, phi, regularized, evaluation(tuning), evaluation(validation)); result_i = result_i + 1
        }
      }
    }
    readr::write_tsv(result, paste0(trait, "/result/", method, '_', ref, '_a', a, '_all_result.txt'))
    
    result_validation = data.frame(matrix(ncol = 3, nrow = 0))
    names(result_validation) = c("seed", "phi", "validation_R2"); result_validation_i = 1
    
    for (seed_i in 1:100) {
      result_tmp = result %>% filter(phi > 0) %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
      result_validation[result_validation_i,] = c(seed_i, result_tmp$phi, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
    }
    phi_best = result_validation %>% group_by(phi) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))
    
    readr::write_tsv(data.frame(phi_best = phi_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                     paste0(trait, "/result/", method, '_', ref, '_a', a, '_validation_result_std.txt'))
  }
}
