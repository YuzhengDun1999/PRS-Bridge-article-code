####### This script is used for evaluation of LDpred2 ###########

outcome = "continuous" # change to continuous if your outcome is continuous 
source("evaluation/evaluation_metrics.R")

library('dplyr')
trait = "Example"
ref_N = "1kg" # Same as in run-methods/ldpred2.R

method = 'ldpred2'
result_dir = paste0(trait, '/ldpred2/')
### for tuning calculate prs
for (chr in c(1:22)) {
  bfile = paste0('validation/tuning/chr', chr)
  for (x in 1:63) {
    prs_output = paste0(result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', ref_N, '_prs_', x,'.txt')
    if(file.exists(paste0(result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', ref_N, '_x_', x, '.txt'))){
      prscode = paste(paste0('plink'),
                      paste0('--score ', result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', ref_N, '_x_', x, '.txt'),
                      paste0(' 1 3 4',' sum --bfile ',bfile),
                      paste0(' --out ', prs_output))
      system(prscode)
    }
  }
}

cov = rbind(bigreadr::fread2(paste0('validation/tuning/', trait, '_cov.txt')),
            bigreadr::fread2(paste0('validation/validation/', trait, '_cov.txt')))
cov = cov %>% select(-FID)
names(cov) = c("IID", paste0("PC", 1:10), "age", "sex", "y")
result = data.frame(matrix(ncol = 4, nrow = 0))
names(result) = c("seed", "tune_i", "tuning_R2", "validation_R2"); result_i = 1
for (x in 1:63) {
  cov$PRS = 0
  for (chr in 1:22) {
    PRS = rbind(bigreadr::fread2(paste0(result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', ref_N, '_prs_', x,".txt.profile")),
                bigreadr::fread2(paste0(result_dir, 'chr', chr, '/ldpred2effect-hm3-EUR-ref_', ref_N, '_validation_prs_', x,".txt.profile")))
    cov = merge(cov, PRS %>% select(IID, SCORESUM), by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:10) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, x, evaluation(tuning), evaluation(validation)); result_i = result_i + 1
  }
}
readr::write_tsv(result, paste0(trait, "/result/", method, '_', ref_N,'_all_result.txt'))

result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "tune_i", "validation_R2"); result_validation_i = 1

for (seed_i in 1:10) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$tune_i, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
tune_i_best = result_validation %>% group_by(tune_i) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(tune_i_best = tune_i_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0(trait, "/result/", method, '_', N, '_validation_result_std.txt'))


