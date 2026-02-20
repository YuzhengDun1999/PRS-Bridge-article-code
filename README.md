
# PRS-Bridge-article-code

This repository contains instructions and scripts to reproduce all the synthetic and real data analysis results in the PRS-Bridge article, benchmarking the following methods: PRS-Bridge, LDpred2, LASSOSUM, and PRS-CS, PRS-CS-Projection, PRS-CS-threshold, and PRS-CS-Regularized.
All the (bash) commands to reproduce the results should be run from the repository's root directory.


## Prerequisites
-   **External software:**  
Download methods in the root directory. 
<br>_"Download methods in the root directory" makes no sense. 
Did you mean something more like "Place downloaded software executable files in the root directory"? (Though this probably isn't quite right b/c many of the software packages don't come with executables.) Be more clear and explicit._ <br>
The necessary LD-related files for the methods should be placed in the 'data/' directory.
	- PRS-Bridge: [https://github.com/YuzhengDun1999/PRSBridge](https://github.com/YuzhengDun1999/PRSBridge) 
    &mdash; download all LD reference data for EUR population by following the instructions provided in the tutorial in the GitHub repo.
	- PRS-CS [https://github.com/getian107/PRScs](https://github.com/getian107/PRScs) 
    &mdash; download LD reference data for EUR population by following method's tutorial in github.
	-  LASSOSUM: [https://github.com/tshmak/lassosum](https://github.com/tshmak/lassosum). 
	- LDpred2: [https://privefl.github.io/bigsnpr/articles/LDpred2.html](https://privefl.github.io/bigsnpr/articles/LDpred2.html) 
    &mdash; download [HapMap3 variants with independent LD blocks](https://doi.org/10.6084/m9.figshare.19213299), this link can also be found in their tutorial.
	- Plink 1.9 [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/) and 2.0 [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/2.0/).
-   **Data access application:** Access to UK Biobank for real trait analyses. 
    (The synthetic data and the summary statistics for disease traits are public available.)

## Data Download and Pre-process
-   **Synthetic data:**  
    Download from  [https://doi.org/10.7910/DVN/COXHAP](https://doi.org/10.7910/DVN/COXHAP).
-   **Real data (disease traits):**
Download the GWAS summary statistics for each trait and get access to the Individual-level genotype and phenotype data UK Biobank (application required). The individual-level data are used to evaluate each method. Run `Rscript process_continuous_dat.R TRAIT PATH_TO_PHENO PATH_TO_PSAM` to generate the corresponding GWAS summary statistics for each disease traits, as well as the disease status and covariate information for individuals. `TRAIT` is the name of the trait (defined in our code below). It must be the same in this whole pipeline. `PATH_TO_PHENO` should be the path to your UK Biobank phenotype file. The file must contain the following field IDs: 'f.eid', 'f.31.0.0', 'f.53.0.0', 'f.21022.0.0', 'f.22009', 'f.22020', 'f.21000.0.0', as well as the fields starting with 'f.20001', 'f.20002', 'f.41270', and 'f.41280'. `PATH_TO_PSAM` should be the path to your UK Biobank '.psam' file. The '.psam' file from any chromosome is acceptable. The resulting GWAS summary statistics will be stored in the file '{TRAIT}/sumdat_Rcov.txt'. The phenotype and covariate information used for tuning and validation will be stored in 'tuning/{TRAIT}_cov.txt' and 'validation/{TRAIT}_cov.txt'.
    -   Breast Cancer: Use the link in [https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-results](https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-results) to download 'oncoarray_bcac_public_release_oct17.txt.gz'. After decompressing the file, place 'oncoarray_bcac_public_release_oct17.txt' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R BC` to get the processed summary statistics stored in `BC/sumdat_Rcov.txt`
    -    Coronary Artery disease: Download 'CARDIoGRAMplusC4D 1000 Genomes-based GWAS – Additive' from [https://cardiogramplusc4d.org/data-downloads/](https://cardiogramplusc4d.org/data-downloads/). Place 'cad.add.160614.website.txt' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R CAD` to get the processed summary statistics stored in `CAD/sumdat_Rcov.txt`
    -    Depression: Download from [https://doi.org/10.6084/m9.figshare.21655784](https://doi.org/10.6084/m9.figshare.21655784). Place 'daner_pgc_mdd_meta_w2_no23andMe_rmUKBB' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R Depression` to get the processed summary statistics stored in `Depression/sumdat_Rcov.txt`.
    -    Inflammatory Bowel Disease: Download 'Latest combined GWAS and Immunochip trans-ancestry' from [https://www.ibdgenetics.org/](https://www.ibdgenetics.org/). Place 'EUR.IBD.gwas_info03_filtered.assoc' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R IBD` to get the processed summary statistics stored in `IBD/sumdat_Rcov.txt`.
     -    Rheumatoid Arthritis: Download 'Eurpean RA GWAS meta-analysis' from [http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz](http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz). Place 'RA_GWASmeta_European_v2.txt.gz' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R RA` to get the processed summary statistics stored in `RA/sumdat_Rcov.txt`.
    <br>_(Afaik, the ACC form isn't going to be published. And you don't provide the download links in the manuscript. Please provide them here explicitly, making it clear which link corresponds to which phenotype. Also describe what file formats they come in.)_ <br>
    _Yuzheng: I am not sure if it's necessary to describe the file formats. They are all text files and process_disease_dat.R read the files with the same name._ <br> 
    _Aki: The thing is, there is no way people would know that files "daner_pgc_mdd_meta_w2_no23andMe_rmUKBB" and "EUR.IBD.gwas_info03_filtered.assoc" are plain text files unless you explicitly say so. Pay attention to the fact how the reproducibility editor got confused about .psam, .bed, .bim and .fam files._ 

    _Yuzheng: For CAD, I have no idea why the summary statistics cannot be downloaded now. It can be directly downloaded several years ago. I find this in their website and I do not know if this is the reason: Please note that CARDIoGRAMplusC4D is currently undertaking the analyses listed under [Ongoing Projects](https://cardiogramplusc4d.org/ongoing-projects/), please contact the relevant PI if you wish to discuss collaborating. Can I just revise the ACC form and move CAD summary statistics to data needed request? This issue is same for inflammatory bowel disease. In their original website, it will appear 'This link has been deleted' if you click the original link._ <br>
    _Aki: If the original datasets were publicly available, can you just re-post your copy of CAD and IBD in a public data repository? (You could even post all the dataset there.)_

-   **Real data (continuous traits):**
	First apply for the access to the individual-level genotype and phenotype data from UK Biobank. 
    <br>_I guess the phenotype data needs to be placed in the folder at `${PATH_TO_PSAM}`? Explain it._ <br>
    Then run `Rscript process_continuous_dat.R ${TRAIT} ${FIELD_ID} ${PATH_TO_PHENO} ${PATH_TO_PSAM}` to generate the corresponding GWAS summary statistics for each continuous trait. 
    <br>_Denote all the variables as I've done above and make sure to introduce them as variables defining specific pieces of info._  <br>
    `TRAIT` is a variable defining the user-specified name of the trait. 
    <br>_What does this mean? The code below seems to require you name the traits specifically as "BMI", "RHR", etc._ <br>
    It must be the same in this whole pipeline. 
    <br>_It is likely easier to follow if you specify the names and use it consistently throughout the instructions._ <br>
    `FIELD_ID` is the corresponding field ID of the trait in the UK Biobank phenotype file. 
    <br>_Unclear what "the corresponding filed ID" corresponds to._  <br>
    `PATH_TO_PHENO` is a variable specifying the path to your UK Biobank phenotype file. The file must contain the following field IDs: 'f.eid', 'f.31.0.0', 'f.21022.0.0', 'f.22009', 'f.22020', and 'f.21000.0.0', as well as the field specified by the second input parameter. `PATH_TO_PSAM` should be the path to your UK Biobank '.psam' file. The '.psam' file from any chromosome is acceptable. The resulting GWAS summary statistics will be stored in the file `TRAIT/sumdat_Rcov.txt`.
     -    BMI: Run `Rscript run-real-data-analysis/process_continuous_dat.R BMI f.21001.0.0 PATH_TO_PHENO PATH_TO_PSAM`. 
     <br>_Note it isn't entirely clear here if `BMI` is a variable or a literal string "BMI"._ <br>
     -    Resting Heart Rate: Run `Rscript process_continuous_dat.R RHR f.102.0.0 PATH_TO_PHENO PATH_TO_PSAM`. 
     -    High-density lipoprotein: Run `Rscript run-real-data-analysis/process_continuous_dat.R HDL f.30760.0.0 PATH_TO_PHENO PATH_TO_PSAM`.
     -    Low-density lipoprotein: Run `Rscript run-real-data-analysis/process_continuous_dat.R LDL f.30780.0.0 PATH_TO_PHENO PATH_TO_PSAM`.
     -    Apolipoprotein A: Run `Rscript run-real-data-analysis/process_continuous_dat.R APOEA f.30630.0.0 PATH_TO_PHENO PATH_TO_PSAM`.
     -    Apolipoprotein B: Run `Rscript run-real-data-analysis/process_continuous_dat.R APOEB f.30640.0.0 PATH_TO_PHENO PATH_TO_PSAM`.

To reproduce the numerical results via the steps outlined below, all the downloaded files should be place in the root directory of this GitHub directory.


## Steps to Reproduce Numerical Results

### PRS-CS's non-convergence behavior in the absence of the ad-hoc constraint on prior variance (Figure 1)

_Make the casing of subsection titles consistent; i.e. use either Title Case or Sentence case consistently._

1. Run the version of PRS-CS without the ad-hoc constraint via the command `python PRS-CS-proj/PRScs_noconverge.py --ref_dir=data/ldblk_1kg_eur --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_noconverge --chrom=22 --write_pst=True --phi=1e-04 --thin=1 --n_iter=800 --n_burnin=1 --seed=1`.
The resulting MCMC output will be stored in the file `data/coef_noconverge_pst_eff_a1_b0.5_phi1e-04_chr22.txt`.
2. Run the version of PRS-CS that uses our projection approach in place of the ad-hoc constraint via
`python PRS-CS-proj/PRScs_proj.py --ref_dir=data/ldblk_1kg_eur --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_proj --chrom=22 --write_pst=True --phi=1e-04 --thin=1 --n_iter=800 --n_burnin=1 --seed=1 --eigenval_rm=0.2`.
The resulting MCMC output will be stored in the file `data/coef_proj_pst_eff_a1_b0.5_phi1e-04_chr22.txt`.
3. Plot the result via `Rscript plot-results/plot_figure1.R`.

### Synthetic Data Results
1.  Run `Rscript run-synthetic-data-analysis/process_synthetic_dat.R ${RHO} ${GA}` to preprocess the GWAS summary statistics under the settings as specified by the variables `RHO` and `GA`. `RHO` controls the proportion of causal SNPs. It takes values 1, 2, or 3, corresponding to causal SNP proportions of 0.01, 0.001, and 0.0005, respectively. 
<br>_It would have been better to call it something like `causal_snp_prop_setting`, but fine._ <br>
`GA` controls the genetic architecture. 
It takes values 1, 4, or 5, corresponding to strong negative selection, no negative selection, and mild negative selection, respectively. 
<br>_This will naturally raises the question of what the value of `2` and `3` would do._ <br>
These parameter choices are taken from the simulation settings of Zhang et al. (2023). Further information on the synthetic data can be found at: https://doi.org/10.7910/DVN/COXHAP. The processed GWAS summary statistics used as input for each method are stored in the directory: `rho${RHO}GA${GA}/sumdat/`.
2. Run `Rscript run-synthetic-data-analysis/PRSBridge.R RHO GA`, `Rscript run-synthetic-data-analysis/PRScs.R RHO GA`, and `Rscript run-synthetic-data-analysis/ldpred2.R RHO GA PATH_TO_1kg`; the input parameters are the same as those in `process_synthetic_dat.R` and `PATH_TO_1kg` is the path to the 1kg individual-level genotype data used to construct the LD reference panel.
3. Run `run-synthetic-data-analysis/evaluation.R RHO GA`; input parameters are the same as `process_synthetic_dat.R`.
4. Plot the results via `Rscript run-synthetic-data-analysis/plot.R`.


### Real Data Results
The computational demands make it unrealistic to reproduce the full results on a local machine, especially for PRS-Bridge and PRS-CS; therefore, the provided R scripts instead generate a set of .sh files that can be submitted to a high-performance computing cluster (HPC) running on SLURM.  The commands for running each method are specified within the generated .sh files. For LASSOSUM and ldpred2, we also provide scripts for directly running the methods, though we recommend similarly executing them on an HPC system.
1. Run the quality control step via `Rscript run-real-data-analysis/sumdat_QC.R TRAIT`, where `TRAIT` should match the input parameter specified in the **Data Download and Pre-process** section. 
<br>_It is a bit weird to use a bold font here. Just write 'in the section "Data Download and Pre-process".'_ <br>
This will generate files in the directory `{TRAIT}/data/`.
2.  Run the PRS methods on the preprocessed data:
     -   Lassosum: Run `Rscript run-real-data-analysis/LASSOSUM.R TRAIT PATH_TO_1kg 1kg` and `Rscript run-real-data-analysis/LASSOSUM.R TRAIT PATH_TO_UKBB ukbb`, where  `TRAIT` should match the input parameter specified in the previous step; `PATH_TO_1kg` and `PATH_TO_UKBB` are the paths to the individual-level genotype data used to construct the corresponding LD reference panels. This will generate coefficients in the directory `{TRAIT}/LASSOSUM/`.
     -   LDpred2: Run `Rscript run-real-data-analysis/ldpred2.R TRAIT PATH_TO_1kg 1kg`, `Rscript run-real-data-analysis/ldpred2.R TRAIT PATH_TO_UKBB ukbb`, and `Rscript run-real-data-analysis/ldpred2.R TRAIT PATH_TO_BLK NBlk`, where  `TRAIT` should match the input parameter specified in the previous step; `PATH_TO_1kg` and `PATH_TO_UKBB` are the paths to the individual-level genotype data used to construct the corresponding LD reference panels; `PATH_TO_BLK` is the path to the downloaded LD reference data with independent block. Please see LDpred2 tutorial for how to process it. This will generate coefficients in the directory `{TRAIT}/chr{1..22}/`.
     -   PRS-CS and its extensions: Run `Rscript run-real-data-analysis/{method}.R TRAIT`, where  `TRAIT` should match the input parameter specified in the previous step, and `{method}` can be chosen from `PRScs` (PRS-CS), `PRScs_proj` (PRS-CS-Projection), `PRScs_regularized` (PRS-CS-Regularized), `PRScs_threshold` (PRS-CS-threshold). This will generate a brunch of script files to run PRS-CS and its extensions in the directory `{TRAIT}/1kg/{method}/run_sh` and `{TRAIT}/ukbb/{method}/run_sh`. 
     <br>_The reason you use an Rscript to generate bash scripts is for the same reasons as you described below for PRS-CS needing the GWAS sample size info?_ <br>
     The jobs will then be automatically submitted via SLURM. 
     <br>_How could jobs be "automatically submitted" (without a user doing anything)? This is confusing._ <br>
     The coefficients will be stored in the directory `{TRAIT}/1kg/{method}/result` and `{TRAIT}/ukbb/{method}/result`.
     -   PRS-Bridge: Run `Rscript run-real-data-analysis/PRSBridge.R TRAIT`, where  `TRAIT` should match the input parameter specified in the previous step. This will generate a brunch of script files to run PRS-Bridge in the directory `{TRAIT}/1kg/Bridge_small/run_sh`,  `{TRAIT}/1kg/Bridge_large/run_sh`, `{TRAIT}/ukbb/Bridge_small/run_sh` and `{TRAIT}/ukbb/Bridge_large/run_sh`. The jobs will then be automatically submitted via SLURM. The coefficients will be stored in the directory `{TRAIT}/1kg/Bridge_small`,  `{TRAIT}/1kg/Bridge_large`, `{TRAIT}/ukbb/Bridge_small` and `{TRAIT}/ukbb/Bridge_large`.
 
**Note:** 
Under the `data` folder, we have provided a preprocessed sample dataset which can be directly used to run the version of PRS-CS without the ad hoc constraint on Chromosome 22 and generate the MCMC samples as plotted in Figure 1 of the manuscript.
<br>_Place this note under the section "PRS-CS's non-convergence behavior in the absence of the ad-hoc constraint on prior variance (Figure 1)"?_ <br>

### Evaluate and Plot Real Data Analysis Results
Having run all the methods on all the datasets, evaluate the performance of each methods via `Rscript evaluate-real-data-analysis/get_sd_{method_name}.R TRAIT OUTCOME` where `method_name` can be chosen from `PRSBridge`, `ldpred2`, `PRScs`, `PRScs_proj`, `PRScs_regularized`, `PRScs_threshold`, and `ldpred2`; `TRAIT` should match the input parameter specified in the previous step; `OUTCOME` should be set to 'continuous' for continuous traits and 'disease' for disease traits. The performance for each method will be stored in `{TRAIT}/result`. After evaluating all methods on all traits, the methods' performances relative to Lassosum can be summarized by running `Rscript evaluate-real-data-analysis/RE_lassosum.R`. The resulting file will be stored in root directory.

Finally, the results can be graphically summarized via `Rscript plot-result/plot.R`, which reproduces Figure 4 and 5.
The other figures in the main manuscript and the supplement can similarly be reproduced by running appropriate scripts under the `plot-result/` folder.

 
## Directory Contents
This is an optional read for those interested in learning more about the specific files in the above steps to reproduce all the results.

### data
This folder contains the following files:
-   `PRScs_sumdat.txt`,  `chr22.bim`: used to generate posterior samples of the regression coefficients used to create Figure 1. Specifically, `PRScs_sumdat.txt` provides the input GWAS summary statistics for original and extension of PRS-CS, and `chr22.bim` contains the individual-level genotype variant information (PLINK `.bim` file) used as the input for original and extension of PRS-CS.
-   `sumdat_Rcov.txt`: Example input summary statistics used by the script in the run-real-data-analysis/sumdat_QC.R to demonstrate how input data should be formatted before preparing data for running all methods.

### PRS-CS-proj

This folder provides our modification of the PRS-CS code, the original version of which is found at [https://github.com/getian107/PRScs](https://github.com/getian107/PRScs). 
The files included are:

-   `PRScs_noconverge.py`: Implements a modified version of PRS-CS that removes the ad hoc constraint on the prior variance. 
See **PRS-CS's non-convergence behavior in the absence of the ad-hoc constraint (Figure 1)** for example usage.
<br>_Again, weird to use a bold font here._ <br>
- `PRScs_threshold.py`: Implements versions of PRS-CS modified to use different values of the upper threshold on the prior variance. 
Example usage: `python PRS-CS-proj/PRScs_threshold.py --ref_dir=data/ldblk_1kg_eur --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_threshold --chrom=22 --phi=1e-04 --seed=1 --threshold=0.1`.
<br>_(I don't understand why you need to have a user run this R scripts to generate the bash scripts.
Why can't you just directly provide the bash scripts generated by this R script and place them in a folder `slurm_scripts/`, for example?
You can then just provide a bash command, e.g. `sbatch ...`, to submit the simulations to SLURM.
That way, it is far more transparent what steps are involved in generating the results.)
Yuzheng: PRS-CS needs the GWAS sample size as input (a scalar). Different traits and different SNPs have different sample sizes (missing rate for different SNPs are different). Instead of manually setting the sample size. I use R script to automatically calculate the mean sample size across SNPs and then use that as input. If the current pipeline will be acceptable by the reviewer, I will keep the current._ 
-   `PRScs_proj.py`: Implements versions of PRS-CS modified to use projected summary statistics. 
See **PRS-CS's non-convergence behavior in the absence of the ad-hoc constraint (Figure 1)** for example usage.
-   `PRScs_Regularized.py`: Implements versions of PRS-CS modified to use Regularized LD matrix. 
Example usage: `python PRS-CS-proj/PRScs_threshold.py --ref_dir=data/ldblk_1kg_eur --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_threshold --chrom=22 --phi=1e-04 --seed=1 --regularized=1`.