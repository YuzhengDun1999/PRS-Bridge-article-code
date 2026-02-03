trait_list = c("Example") # a list of trait to get relative efficiency, baseline is LASSOSUM.

get_validation = function(dat){
  result_validation = data.frame(matrix(ncol = 2, nrow = 0))
  names(result_validation) = c("seed", "validation_R2"); result_validation_i = 1
  for (seed_i in 1:100) {
    result_tmp = dat %>% filter(seed == seed_i) %>% filter(tuning_R2 == max(tuning_R2))
    result_validation[result_validation_i,] = c(seed_i, result_tmp$validation_R2); result_validation_i = result_validation_i + 1
  }
  return(result_validation)
}
get_RR = function(dat1, dat2){
  RR = dat2 / dat1
  return(c(mean(RR), sqrt(var(RR))))
}
library(dplyr)
result = data.frame(matrix(ncol = 5, nrow = 0))
names(result) = c("method", "trait", "ref", "RR", "std"); result_i = 1
for (trait in trait_list) {
  Lassosum = get_validation(bigreadr::fread2(paste0(trait, "/result/Lassosum_N5000_all_result.txt")))
  ldpred2_blk = get_validation(bigreadr::fread2(paste0(trait, "/result/ldpred2_Nblk_provided_all_result.txt")))
  ldpred2_N5000 = get_validation(bigreadr::fread2(paste0(trait, "/result/ldpred2_N5000_all_result.txt")))
  PRScs = get_validation(bigreadr::fread2(paste0(trait, "/result/PRScs_ukbb_all_result.txt")))
  Bridge_small = get_validation(bigreadr::fread2(paste0(trait, "/result/Bridge_small_ukbb_all_result.txt")))
  Bridge_large = get_validation(bigreadr::fread2(paste0(trait, "/result/Bridge_large_ukbb_all_result.txt")))
  result[result_i, ] = c("ldpred2_N5000", trait, "ukbb", get_RR(Lassosum$validation_R2, ldpred2_N5000$validation_R2)); result_i = result_i + 1
  result[result_i, ] = c("ldpred2_blk", trait, "ukbb", get_RR(Lassosum$validation_R2, ldpred2_blk$validation_R2)); result_i = result_i + 1
  result[result_i, ] = c("PRScs", trait, "ukbb", get_RR(Lassosum$validation_R2, PRScs$validation_R2)); result_i = result_i + 1
  result[result_i, ] = c("Bridge_small", trait, "ukbb", get_RR(Lassosum$validation_R2, Bridge_small$validation_R2)); result_i = result_i + 1
  result[result_i, ] = c("Bridge_large", trait, "ukbb", get_RR(Lassosum$validation_R2, Bridge_large$validation_R2)); result_i = result_i + 1
  Lassosum = get_validation(bigreadr::fread2(paste0(trait, "/result/Lassosum_N1kg_all_result.txt")))
  ldpred2_N5000 = get_validation(bigreadr::fread2(paste0(trait, "/result/ldpred2_1kg_all_result.txt")))
  PRScs = get_validation(bigreadr::fread2(paste0(trait, "/result/PRScs_1kg_all_result.txt")))
  Bridge_small = get_validation(bigreadr::fread2(paste0(trait, "/result/Bridge_small_1kg_all_result.txt")))
  Bridge_large = get_validation(bigreadr::fread2(paste0(trait, "/result/Bridge_large_1kg_all_result.txt")))
  result[result_i, ] = c("ldpred2_N5000", trait, "1kg", get_RR(Lassosum$validation_R2, ldpred2_N5000$validation_R2)); result_i = result_i + 1
  result[result_i, ] = c("PRScs", trait, "1kg", get_RR(Lassosum$validation_R2, PRScs$validation_R2)); result_i = result_i + 1
  result[result_i, ] = c("Bridge_small", trait, "1kg", get_RR(Lassosum$validation_R2, Bridge_small$validation_R2)); result_i = result_i + 1
  result[result_i, ] = c("Bridge_large", trait, "1kg", get_RR(Lassosum$validation_R2, Bridge_large$validation_R2)); result_i = result_i + 1
}
readr::write_tsv(result, "RR.txt")
