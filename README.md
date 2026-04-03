


# PRS-Bridge-article-code

This repository contains instructions and scripts to reproduce all the synthetic and real data analysis results in the PRS-Bridge article, benchmarking the following methods: PRS-Bridge, LDpred2, LASSOSUM, and PRS-CS, PRS-CS-Projection, PRS-CS-threshold, and PRS-CS-Regularized.
All the (bash) commands to reproduce the results should be run from the repository's root directory.


## Prerequisites
-   **External software:**  
Download PRS-Bridge, PRS-CS, PLINK 1.9, and PLINK 2.0 from the links below and place the packages and executable files in the root directory. 
Install LASSOSUM and LDpred2 in R according to the instructions in their respective tutorials.
The necessary LD-related files for the methods should be placed in the 'data/' directory.
	- PRS-Bridge: [https://github.com/YuzhengDun1999/PRSBridge](https://github.com/YuzhengDun1999/PRSBridge).
	<br> Additionally, download pre-estimated LD matrices for EUR population required by PRS-Bridge using the following links; the downloaded files can be uncompressed via a command `tar -zxvf`: 
		- [1000 Genomes Small Block](https://zenodo.org/records/18706275)
		- [UK Biobank Large Block](https://zenodo.org/records/18706275)
		- [1000 Genomes Large Block](https://zenodo.org/records/18673493)
		- [UK Biobank Small Block](https://zenodo.org/records/18673493 "AFR reference")
	- PRS-CS [https://github.com/getian107/PRScs](https://github.com/getian107/PRScs).
		<br> Additionally, download pre-estimated LD matrices for EUR population required by PRS-CS using the following links; the downloaded files can be uncompressed via a command `tar -zxvf`: 
		- [UK Biobank ](https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=0 "EUR reference")
		- [1000 Genomes](https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0 "EUR reference")
	-  LASSOSUM: [https://github.com/tshmak/lassosum](https://github.com/tshmak/lassosum). 
	- LDpred2: [https://privefl.github.io/bigsnpr/articles/LDpred2.html](https://privefl.github.io/bigsnpr/articles/LDpred2.html).
    <br> Additionally, download pre-estimated LD matrix with block structure required by LDpred2 using the following link: [HapMap3 variants with independent LD blocks](https://doi.org/10.6084/m9.figshare.19213299).
	- Plink 1.9 [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/) and 2.0 [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/2.0/).
-   **Data access application:** Access to UK Biobank for real trait analyses. 
    (The synthetic data and the summary statistics for disease traits are public available.)

**Note:** 
For the pre-estimated LD matrices required by LDpred2 and PRS-CS, the download links are not maintained by us and are only confirmed to be current as of March, 2026 ---
if necessary, please refer to the corresponding methods’ websites for the latest download links.


## Data Download and Preprocess
-   **Synthetic data:**  
Download the PLINK files (`.bim`, `.bed`, `.fam`) with the prefix `EUR_`, the phenotype files named `EUR_pheno_rho*_size_4_GA_{1,4,5}`, and the GWAS summary statistics files named `EUR_summary_rho*_size_4_GA_{1,4,5}` from  [https://doi.org/10.7910/DVN/COXHAP](https://doi.org/10.7910/DVN/COXHAP).

-   **Real data (disease traits):**
Download the GWAS summary statistics for each trait and get access to the Individual-level genotype and phenotype data UK Biobank. The individual-level data are used to evaluate each method. Run `Rscript run-real-data-analysis/process_disease_dat.R ${TRAIT} ${PATH_TO_PHENO} ${PATH_TO_GENO}` to generate the corresponding GWAS summary statistics for each disease traits, as well as the disease status and covariate information for individuals. 
 `TRAIT` is a variable specifying the trait name, to be set as one of the following: `"BC"`, `"CAD"`, `"Depression"`, `"IBD"`, and `"RA"`. 
 `PATH_TO_PHENO` is a variable specifying the path to your UK Biobank phenotype file in the .rds format. The file must contain the following field IDs:  'f.eid', 'f.31.0.0', 'f.53.0.0', 'f.21022.0.0', 'f.22009', 'f.22020', 'f.21000.0.0', as well as the fields starting with 'f.20001', 'f.20002', 'f.41270', and 'f.41280'. 
    `PATH_TO_GENO` is a variable specifying the basename of the PLINK files for your UK Biobank data, i.e. the files should be named `${PATH_TO_GENO}.{bim,bed,fam}`. 
    PLINK files from any chromosome are acceptable. 
    The resulting GWAS summary statistics will be stored in the file `${TRAIT}/sumdat.txt`. The phenotype and covariate information used for tuning and validation will be stored in `tuning/${TRAIT}_cov.txt` and `validation/${TRAIT}_cov.txt`.
    -   Breast Cancer: Use the link in [https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-results](https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-results) and download 'oncoarray_bcac_public_release_oct17.txt.gz'. After decompressing the file, place the text file 'oncoarray_bcac_public_release_oct17.txt' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R BC ${PATH_TO_PHENO} ${PATH_TO_GENO}` to get the processed summary statistics stored in `BC/sumdat.txt`
    -    Coronary Artery disease: Download 'CARDIoGRAMplusC4D 1000 Genomes-based GWAS – Additive' from [https://cardiogramplusc4d.org/data-downloads/](https://cardiogramplusc4d.org/data-downloads/). Place the text file  'cad.add.160614.website.txt' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R CAD ${PATH_TO_PHENO} ${PATH_TO_GENO}` to get the processed summary statistics stored in `CAD/sumdat.txt`
    -    Depression: Download from [https://doi.org/10.6084/m9.figshare.21655784](https://doi.org/10.6084/m9.figshare.21655784). Place the text file  'daner_pgc_mdd_meta_w2_no23andMe_rmUKBB' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R Depression ${PATH_TO_PHENO} ${PATH_TO_GENO}` to get the processed summary statistics stored in `Depression/sumdat.txt`.
    -    Inflammatory Bowel Disease: Download 'Latest combined GWAS and Immunochip trans-ancestry' from [https://www.ibdgenetics.org/](https://www.ibdgenetics.org/). Place the text file 'EUR.IBD.gwas_info03_filtered.assoc' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R IBD ${PATH_TO_PHENO} ${PATH_TO_GENO}` to get the processed summary statistics stored in `IBD/sumdat.txt`.
     -    Rheumatoid Arthritis: Download 'Eurpean RA GWAS meta-analysis' from [http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz](http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz). Place the text file 'RA_GWASmeta_European_v2.txt.gz' in the root directory. Then run `Rscript run-real-data-analysis/process_disease_dat.R RA ${PATH_TO_PHENO} ${PATH_TO_GENO}` to get the processed summary statistics stored in `RA/sumdat.txt`.

**Note:** 
For coronary artery disease and inflammatory bowel disease, we obtained the original GWAS datasets in 2022. These datasets are no longer publicly available as of February 2026. The corresponding references are provided in Table 1 of the manuscript. Interested researchers may contact the corresponding authors of the original studies to request access.

-   **Real data (continuous traits):**
	First apply for the access to the individual-level genotype and phenotype data from UK Biobank. 
    Then run `Rscript run-real-data-analysis/process_continuous_dat.R ${TRAIT} ${FIELD_ID} ${PATH_TO_PHENO} ${PATH_TO_GENO}` to generate the corresponding GWAS summary statistics for each continuous trait. 
    `TRAIT` is a variable specifying the trait name, to be set as one of the followings: `"BMI"`, `"RHR"`, `"HDL"`, `"LDL"`, `"APOEA"`, and `"APOEB"`. 
    `FIELD_ID` is a variable denoting the UK Biobank field ID corresponding to the specified `${TRAIT}` in the phenotype file.
    `PATH_TO_PHENO` is a variable specifying the path to your UK Biobank phenotype file stores in the .rds format. The file must contain the following field IDs: 'f.eid', 'f.31.0.0', 'f.21022.0.0', 'f.22009', 'f.22020', and 'f.21000.0.0', as well as the field ID specified by `${FIELD_ID}`.
    `PATH_TO_GENO` is a variable specifying the basename of the PLINK files for your UK Biobank data, i.e. the files should be named `${PATH_TO_GENO}.{bim,bed,fam}`. 
    The resulting GWAS summary statistics will be stored in the file `${TRAIT}/sumdat.txt`. The phenotype and covariate information used for tuning and validation will be stored in `tuning/${TRAIT}_cov.txt` and `validation/${TRAIT}_cov.txt`.
     -    BMI: Run `Rscript run-real-data-analysis/process_continuous_dat.R BMI f.21001.0.0 ${PATH_TO_PHENO} ${PATH_TO_GENO}`. 
     -    Resting Heart Rate: Run `Rscript process_continuous_dat.R RHR f.102.0.0 ${PATH_TO_PHENO} ${PATH_TO_GENO}`. 
     -    High-density lipoprotein: Run `Rscript run-real-data-analysis/process_continuous_dat.R HDL f.30760.0.0 ${PATH_TO_PHENO} ${PATH_TO_GENO}`.
     -    Low-density lipoprotein: Run `Rscript run-real-data-analysis/process_continuous_dat.R LDL f.30780.0.0 ${PATH_TO_PHENO} ${PATH_TO_GENO}`.
     -    Apolipoprotein A: Run `Rscript run-real-data-analysis/process_continuous_dat.R APOEA f.30630.0.0 ${PATH_TO_PHENO} ${PATH_TO_GENO}`.
     -    Apolipoprotein B: Run `Rscript run-real-data-analysis/process_continuous_dat.R APOEB f.30640.0.0 ${PATH_TO_PHENO} ${PATH_TO_GENO}`.

To reproduce the numerical results via the steps outlined below, all the downloaded files should be place in the root directory of this GitHub directory.


## Steps to Reproduce Numerical Results

### PRS-CS's Non-convergence Behavior in the Absence of the Ad-hoc Constraint on Prior Variance (Figure 1)

1. Run the version of PRS-CS without the ad-hoc constraint via the command `python PRS-CS-proj/PRScs_noconverge.py --ref_dir=data/ldblk_1kg_eur --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_noconverge --chrom=22 --write_pst=True --phi=1e-04 --thin=1 --n_iter=800 --n_burnin=1 --seed=1`.
The resulting MCMC output will be stored in the file `data/coef_noconverge_pst_eff_a1_b0.5_phi1e-04_chr22.txt`.
2. Run the version of PRS-CS that uses our projection approach in place of the ad-hoc constraint via
`python PRS-CS-proj/PRScs_proj.py --ref_dir=data/ldblk_1kg_eur --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_proj --chrom=22 --write_pst=True --phi=1e-04 --thin=1 --n_iter=800 --n_burnin=1 --seed=1 --eigenval_rm=0.2`.
The resulting MCMC output will be stored in the file `data/coef_proj_pst_eff_a1_b0.5_phi1e-04_chr22.txt`.
3. Plot the result via `Rscript plot-result/plot_fig1.R`.

**Note:** 
Under the `data` folder, we have provided a preprocessed sample dataset which can be directly used to reproduce Figure 1 following the above steps.



### Synthetic Data Results

Due to the computational efforts involved, reproducing the full results in a reasonable amount of time requires running the tasks in parallel on a computing cluster.
For PRS-Bridge and PRS-CS, therefore, we provide bash scripts to run the methods by submitting them via `sbatch` to a cluster operating under Simple Linux Utility for Resource Management (SLURM). 
For LASSOSUM and ldpred2, we provide R scripts for directly running the methods, though we recommend similarly executing them on a computing cluster.

1.  Run `Rscript run-synthetic-data-analysis/process_synthetic_dat.R ${RHO} ${GA}` to preprocess the GWAS summary statistics under the settings as specified by the variables `RHO` and `GA`. 
`RHO` is a variable specifying the proportion of causal SNPs. It takes values 1, 2, or 3, corresponding to causal SNP proportions of 0.01, 0.001, and 0.0005, respectively. 
`GA` is a variable specifying the genetic architecture. It takes values 1, 4, or 5, corresponding to strong negative selection, no negative selection, and mild negative selection, respectively. 
The values of 2 and 3 correspond to multi-ancestry settings, which are not considered in our article. 
These parameter choices are taken from the simulation settings of Zhang et al. (2023). Further information on the synthetic data can be found at: https://doi.org/10.7910/DVN/COXHAP. The processed GWAS summary statistics used as input for each method are stored in the directory: `rho${RHO}GA${GA}/sumdat/`.
3. Run `sbatch run-synthetic-data-analysis/PRSBridge.sh ${RHO} ${GA}`, `sbatch run-synthetic-data-analysis/PRScs.sh ${RHO} ${GA}`, and `Rscript run-synthetic-data-analysis/ldpred2.R ${RHO} ${GA} ${PATH_TO_1kg}`; the input parameters are the same as those in `process_synthetic_dat.R` and `PATH_TO_1kg` is a variable specifying the path to the 1kg individual-level genotype data used to construct the LD reference panel.
4. Run `Rscript run-synthetic-data-analysis/evaluation.R ${RHO} ${GA}`; input parameters are the same as `process_synthetic_dat.R`.
5. Plot the results via `Rscript run-synthetic-data-analysis/plot.R`.


### Real Data Results

As with the synthetic data results, reproducing the full results here requires significant computational efforts.
Therefore, for PRS-Bridge, PRS-CS, and its modifications, we provide bash scripts to submit the tasks to a computing cluster via `sbatch`. 

1. Run the quality control step via `Rscript run-real-data-analysis/sumdat_QC.R ${TRAIT}`, where the variable `TRAIT` is as described in the "Data Download and Preprocess" section.
The resulting files will be stored in the directory `${TRAIT}/data/`.
2.  Run the PRS methods on the preprocessed data:
     -   Lassosum: Run `Rscript run-real-data-analysis/LASSOSUM.R ${TRAIT} ${PATH_TO_1kg} 1kg` and `Rscript run-real-data-analysis/LASSOSUM.R ${TRAIT} ${PATH_TO_UKBB} ukbb`, where `PATH_TO_1kg` and `PATH_TO_UKBB` are variables specifying the paths to the individual-level genotype data used to construct the corresponding LD reference panels. The estimated regression coefficients will be stored in the directory `${TRAIT}/LASSOSUM/`.
     -   LDpred2: Run `Rscript run-real-data-analysis/ldpred2.R ${TRAIT} ${PATH_TO_1kg} 1kg`, `Rscript run-real-data-analysis/ldpred2.R ${TRAIT} ${PATH_TO_UKBB} ukbb`, and `Rscript run-real-data-analysis/ldpred2.R ${TRAIT} ${PATH_TO_BLK} NBlk`, where `PATH_TO_1kg` and `PATH_TO_UKBB` are variables specifying the paths to the individual-level genotype data used to construct the corresponding LD reference panels, and `PATH_TO_BLK` is a variable specifying the path to the downloaded LD reference data with independent block. The estimated regression coefficients will be stored in the directory `${TRAIT}/chr{1..22}/`.
     -   PRS-CS and its modifications: Run `sbatch run-real-data-analysis/${method}.sh ${TRAIT} ${PATH_TO_GENO}`, where `${method}` should be chosen from `PRScs`, `PRScs_proj`, `PRScs_regularized`, and `PRScs_threshold`. `PATH_TO_GENO` specifies the basename of the PLINK files for your UK Biobank data, i.e. the files should be named `${PATH_TO_GENO}.{bim,bed,fam}`. This command submits SLURM array jobs to run PRS-CS and its modifications in parallel. The estimated regression coefficients will be stored in the directory `${TRAIT}/1kg/${method}/result` and `${TRAIT}/ukbb/${method}/result`.
     -   PRS-Bridge: Run `sbatch run-real-data-analysis/PRSBridge.sh ${TRAIT}`. This command submits SLURM array jobs to run PRS-Bridge in parallel. The estimated regression coefficients will be stored in the directory `${TRAIT}/1kg/Bridge_small`,  `${TRAIT}/1kg/Bridge_large`, `${TRAIT}/ukbb/Bridge_small`, and `${TRAIT}/ukbb/Bridge_large`.
 
**Note:** 
The `example` folder contains "mock" datasets that have the same structure as the real datasets after you download and preprocess them following the instructions in the "Data Download and Preprocess" section.
These mock datasets are provided to illustrate the real-data analysis workflow:
they can be used to run the script `example-run-real-data-analysis.sh` &mdash; like the real datasets can be used to run the scripts under `run-real-data-analysis/` &mdash; to generate SNP effect size estimates from all the methods.
More precisely, `example/sumdat.txt` provides mock input GWAS summary statistics, and `example/UKBB` and `example/1kg` provide mock individual-level genotype data used to construct LD reference data. 
Running `bash example-run-real-data-analysis.sh` generates the SNP effect size estimates based on the mock datasets.
 
### Evaluate and Plot Real Data Analysis Results
Having run all the methods on all the datasets, evaluate their performances via `Rscript evaluate-real-data-analysis/get_sd_${method_name}.R ${TRAIT} ${OUTCOME} ${PATH_TO_GENO}` 
where `${method_name}` should be chosen from `PRSBridge`, `ldpred2`, `PRScs`, `PRScs_proj`, `PRScs_regularized`, `PRScs_threshold`, and `ldpred2`; 
`TRAIT` is as described in the earlier sections; 
`OUTCOME` should be set to 'continuous' for continuous traits and to 'disease' for binary disease traits; 
`PATH_TO_GENO` specifies the basename of the PLINK files for your UK Biobank data, i.e. the files should be named `${PATH_TO_GENO}.{bim,bed,fam}`. 
The summary of each method's performance will be stored in `${TRAIT}/${method}/result`. After evaluating all the methods' performances, the methods' performances relative to Lassosum can be summarized by running `Rscript evaluate-real-data-analysis/RE_lassosum.R`. The resulting file will be stored in the directory `${trait}/result/`.

Finally, the results are graphically summarized via `Rscript plot-result/plot.R`, which reproduces Figure 4 and 5.
The other figures in the main manuscript and the supplement can similarly be reproduced by running appropriate scripts under the `plot-result/` folder.

**Note:** 
The tuning and validation data necessary for this evaluation step are generated under the folders `tuning` and `validation` when you run `run-real-data-analysis/process_disease_dat.R` and `run-real-data-analysis/process_continuous_dat.R` as described in the "Data Download and Preprocess" section.
The `tuning` and `validation` folders as provided only contain "mock" individual-level covariate data for illustration purposes;
these data are used for evaluating the mock data analysis results produced by the script `example-evaluate-real-data-analysis.sh` in the previous step.
Running `bash example-evaluate-real-data-analysis.sh` evaluates and visualizes the mock data analysis results.
 
 
## Directory Contents
This is an optional read for those interested in learning more about the specific files in the above steps to reproduce all the results.

### data
This folder contains the following files:
-   `PRScs_sumdat.txt`,  `chr22.bim`: used to generate posterior samples of the regression coefficients used to create Figure 1. Specifically, `PRScs_sumdat.txt` provides the input GWAS summary statistics for the original and modified versions of PRS-CS, and `chr22.bim` provides the individual-level genotype variant information used as an input for the original and modified versions of PRS-CS.


### PRS-CS-proj

This folder provides our modification of the PRS-CS code, the original version of which is found at [https://github.com/getian107/PRScs](https://github.com/getian107/PRScs). 
The files included are:

-   `PRScs_noconverge.py`: Implements a modified version of PRS-CS that removes the ad hoc constraint on the prior variance. 
See the section "PRS-CS's Non-convergence Behavior in the Absence of the Ad-hoc Constraint on Prior Variance (Figure 1)" for example usage.
- `PRScs_threshold.py`: Implements versions of PRS-CS modified to use different values of the upper threshold on the prior variance. 
Example usage: `python PRS-CS-proj/PRScs_threshold.py --ref_dir=data/ldblk_1kg_eur --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_threshold --chrom=22 --phi=1e-04 --seed=1 --threshold=0.1`.
-   `PRScs_proj.py`: Implements versions of PRS-CS modified to use projected summary statistics. 
See the section "PRS-CS's Non-convergence Behavior in the Absence of the Ad-hoc Constraint on Prior Variance (Figure 1)" for example usage.
-   `PRScs_Regularized.py`: Implements versions of PRS-CS modified to use a regularized LD matrix. 
Example usage: `python PRS-CS-proj/PRScs_threshold.py --ref_dir=data/ldblk_1kg_eur --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_threshold --chrom=22 --phi=1e-04 --seed=1 --regularized=1`.