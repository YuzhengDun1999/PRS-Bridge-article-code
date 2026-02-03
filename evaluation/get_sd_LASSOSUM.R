####### This script is used for evaluation of ALSSOSUM ###########

outcome = "continuous" # change to continuous if your outcome is continuous 
source("evaluation/evaluation_metrics.R")

############## Lassosum ##############
library('dplyr')
trait = "Example"
ref_N = "1kg" #Same as ref_N in run-methods/LASSOSUM.R


method = 'Lassosum'
result_dir = paste0(trait, '/LASSOSUM/')

cov = rbind(bigreadr::fread2(paste0('validation/tuning/', trait, '_cov.txt')),
            bigreadr::fread2(paste0('validation/validation/', trait, '_cov.txt')))
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
    PRS_tmp = PRS[, c(2, x + 5)]; names(PRS_tmp) = c("IID", "SCORESUM")
    cov = merge(cov, PRS_tmp, by = "IID") %>% mutate(PRS = PRS + SCORESUM) %>% select(-SCORESUM)
  }
  for (seed in 1:10) {
    set.seed(seed)
    idx = sample(1:nrow(cov), size = round(nrow(cov) / 2))
    tuning = cov[idx, ]; validation = cov[-idx, ]
    result[result_i,] = c(seed, x, evaluation(tuning), evaluation(validation)); result_i = result_i + 1
  }
}
readr::write_tsv(result, paste0(trait, "/result/", method, '_', N,'_all_result.txt'))

result_validation = data.frame(matrix(ncol = 3, nrow = 0))
names(result_validation) = c("seed", "tune_i", "validation_R2"); result_validation_i = 1

for (seed_i in 1:10) {
  result_tmp = result %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
  result_validation[result_validation_i,] = c(seed_i, result_tmp$tune_i, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
}
tune_i_best = result_validation %>% group_by(tune_i) %>% summarise(total_count = n()) %>% as.data.frame() %>% filter(total_count == max(total_count))

readr::write_tsv(data.frame(tune_i_best = tune_i_best[1,1], R2 = mean(result_validation$validation_R2), std = sqrt(var(result_validation$validation_R2))), 
                 paste0( trait, "/result/", method, '_', N, '_validation_result_std.txt'))


