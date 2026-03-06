###### provide example code to run pipeline under the section Real Data Results

TRAIT="example"
PATH_TO_1kg="example/1kg"
PATH_TO_UKBB="example/UKBB"

Rscript run-real-data-analysis/sumdat_QC.R ${TRAIT}
Rscript run-real-data-analysis/LASSOSUM.R ${TRAIT} ${PATH_TO_1kg} 1kg
Rscript run-real-data-analysis/LASSOSUM.R ${TRAIT} ${PATH_TO_UKBB} ukbb

Rscript run-real-data-analysis/ldpred2.R ${TRAIT} ${PATH_TO_1kg} 1kg
Rscript run-real-data-analysis/ldpred2.R ${TRAIT} ${PATH_TO_UKBB} ukbb
Rscript run-real-data-analysis/ldpred2.R ${TRAIT} ${PATH_TO_BLK} NBlk

PATH_TO_GENO="example/UKBB"
sbatch run-real-data-analysis/PRScs.sh ${TRAIT} ${PATH_TO_GENO}
sbatch run-real-data-analysis/PRScs_proj.sh ${TRAIT} ${PATH_TO_GENO}
sbatch run-real-data-analysis/PRScs_regularized.sh ${TRAIT} ${PATH_TO_GENO}
sbatch run-real-data-analysis/PRScs_threshold.sh ${TRAIT} ${PATH_TO_GENO}
sbatch run-real-data-analysis/PRScs_threshold.sh ${TRAIT} ${PATH_TO_GENO}

sbatch run-real-data-analysis/PRSBridge.sh ${TRAIT}

