####### This script is used for evaluation of PRS-Brdge ###########

outcome = "continuous" # change to continuous if your outcome is continuous 
source("evaluation/evaluation_metrics.R")

#### Example: To evaluate using 1kg small ref in run-methods/PRSBridge.R for BMI
trait = "Example" # input your trait used in run-methods/PRSBridge.R
LD_size = "small" # input what LD_size you use in run-methods/PRSBridge.R, can be "large", and "small"
ref = "1kg" # input what LD reference data you use in run-methods/PRSBridge.R

######### change your alpha_list and percent_list here if you use different parameters in run-methods/PRSBridge.R
alpha_list = c(0.5, 0.25, 0.125)
percent_list = c(0.2, 0.4, 0.6, 0.8)


system(paste0("mkdir -p ", trait, "/result/"))
method = paste0("Bridge_", LD_size)
result_dir = paste0(trait, "/", ref, "/", method, "/result/")
for (chr in c(1:22)) {
  for (alpha_i in 1:length(alpha_list)) {
    for (percent_i in 1:length(percent_list)) {
      bfile = paste0('tuning/chr', chr) #### path to your tuning bfile
      alpha = alpha_list[alpha_i]; percent = percent_list[percent_i]
      output_dir = paste0(result_dir, "chr", chr, "/percent",percent,"/alpha",alpha)
      prs_output = paste0(output_dir, '/prsresult.txt')
      prsscore = paste0(output_dir, '/coef.txt')
      prscode = paste(paste0('plink'),
                      paste0('--score ', prsscore),
                      paste0(' 1 3 4 sum --bfile ', bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
      bfile = paste0('validation/chr', chr) #### path to your validation bfile
      prs_output = paste0(output_dir, '/prsresult_validation.txt')
      prsscore = paste0(output_dir, '/coef.txt')
      prscode = paste(paste0('plink'),
                      paste0('--score ', prsscore),
                      paste0(' 1 3 4 sum --bfile ', bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
    }
  }
}

library(dplyr)
cov = rbind(bigreadr::fread2(paste0('tuning/', trait, '_cov.txt')), #### path to your covariate file for tuning
            bigreadr::fread2(paste0('validation/', trait, '_cov.txt'))) #### path to your covariate file for validation
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
    for (seed in 1:100) {
      set.seed(seed)
      idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
      tuning = cov[idx, ]; validation = cov[-idx, ]
      result[result_i,] = c(seed, alpha, percent, evaluation(tuning), evaluation(validation)); result_i = result_i + 1
    }
  }
}
readr::write_tsv(result, paste0(trait, "/result/", method, '_', ref, '_all_result.txt'))

result_validation = data.frame(matrix(ncol = 4, nrow = 0))
names(result_validation) = c("seed", "alpha", "percent", "validation_R2"); result_validation_i = 1

for (seed_i in 1:100) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$alpha, result_tmp$percent, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
alpha_best = result_validation %>% group_by(alpha) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(alpha_best = alpha_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0(trait, "/result/", method, '_', ref, '_validation_result_std.txt'))



############## PRS-Bridge-auto ###########
alpha_list = c("auto")
for (chr in c(1:22)) {
  for (alpha_i in 1:length(alpha_list)) {
    for (percent_i in 1:length(percent_list)) {
      bfile = paste0('validation/tuning/chr', chr)
      alpha = alpha_list[alpha_i]; percent = percent_list[percent_i]
      output_dir = paste0(result_dir, "chr", chr, "/percent",percent,"/alpha",alpha)
      prs_output = paste0(output_dir, '/prsresult.txt')
      prsscore = paste0(output_dir, '/coef.txt')
      prscode = paste(paste0('plink'),
                      paste0('--score ', prsscore),
                      paste0(' 1 3 4 sum --bfile ', bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
      bfile = paste0('validation/validation/chr', chr)
      prs_output = paste0(output_dir, '/prsresult_validation.txt')
      prsscore = paste0(output_dir, '/coef.txt')
      prscode = paste(paste0('plink'),
                      paste0('--score ', prsscore),
                      paste0(' 1 3 4 sum --bfile ', bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
    }
  }
}
  
library(dplyr)
cov = rbind(bigreadr::fread2(paste0('tuning/', trait, '_cov.txt')), #### path to your covariate file for tuning
            bigreadr::fread2(paste0('validation/', trait, '_cov.txt'))) #### path to your covariate file for validation
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
      result[result_i,] = c(seed, alpha, percent, evaluation(tuning), evaluation(validation)); result_i = result_i + 1
    }
  }
}
readr::write_tsv(result, paste0(trait, "/result/", method, '_', ref,'_all_result.txt'))

result = result %>% mutate(percent = as.numeric(percent), tuning_R2 = as.numeric(tuning_R2), validation_R2 = as.numeric(validation_R2))
result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("percent", "R2", "std"); result_validation_i = 1
for (percent_i in percent_list) {
  result_tmp = result %>% filter(percent == percent_i)
  result_validation[result_validation_i,] = c(percent_i, mean(result_tmp$validation_R2), sqrt(var(result_tmp$validation_R2))); result_validation_i = result_validation_i + 1
}
readr::write_tsv(result_validation, paste0( trait, "/result/", method, "_auto", '_', ref, '_validation_result_std.txt'))
