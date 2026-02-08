# PRS-Bridge-article-code

This repository contains instruction and scripts to reproduce all the synthetic and real data analysis results in the PRS-Bridge article, benchmarking the following methods: PRS-Bridge, LDpred2, LASSOSUM, and PRS-CS, PRS-CS-Projection, PRS-CS-threshold, and PRS-CS-Regularized.
All the (bash) commands to reproduce the results should be run from the repository's root directory.


## Prerequisites
-   **External software:**  PRS-Bridge [https://github.com/YuzhengDun1999/PRSBridge](https://github.com/YuzhengDun1999/PRSBridge); PRS-CS [https://github.com/getian107/PRScs](https://github.com/getian107/PRScs); Plink 1.9 [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/) and 2.0 [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/2.0/); LDpred2 [https://privefl.github.io/bigsnpr/articles/LDpred2.html](https://privefl.github.io/bigsnpr/articles/LDpred2.html); LASSOSUM [https://github.com/tshmak/lassosum](https://github.com/tshmak/lassosum).
-   **Data:**  UK Biobank access for real trait analyses; synthetic data and summary statistics for disease traits public available.

## Data Download
-   **Synthetic data:**  
    Download from  [https://doi.org/10.7910/DVN/COXHAP](https://doi.org/10.7910/DVN/COXHAP).
-   **Real data:**
    -   Disease traits: check our manuscript and acc_form for reference and download link.
    <br>_(Afaik, the ACC form isn't going to be published. And you don't provide the download links in the manuscript. Please provide them here explicitly, making it clear which link corresponds to which phenotype. Also describe what file formats they come in.)_
    -   Continuous traits: Individual-level genotype and phenotype data UK Biobank (application required). Follow the instruction in run-real-data-analysis/process_continuous_dat.R to get summary statistics.

To reproduce the numerical results via the steps outlined below, all the downloaded files should be place in the root directory of this GitHub directory.
<br>_(This will create a mess, but it looks like you hardcoded the relative paths this way. It makes more sense to place all the downloaded files in a folder like `raw_data` and the processed ones in `processed_data`, but it is also acceptable to keep things as they are.)_


## Steps to Reproduce Numerical Results

### PRS-CS's non-convergence behavior in the absence of the ad-hoc constraint (Figure 1)

_(It is confusing that your current file `plot_fig1.R` actually runs a Python script to generate the necessary outputs, rather than simply making the plot out of existing outputs. Adopt the workflow as blow and modify the R script appropriately.)_

1. Run the version of PRS-CS without the ad-hoc constraint via the command `python PRS-CS-proj/PRScs_noconverge.py --ref_dir=LD_REFERNCE_DIR --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_noconverge --chrom=22 --write_pst=True --phi=1e-04 --thin=1 --n_iter=800 --n_burnin=1 --seed=1` where `LD_REFERNCE_DIR` should be set as .... 
The resulting MCMC output will be stored it in the file `...`.
_(Can you just `--ref_dir="data/...` instead?)_
2. Run the version of PRS-CS without the ad-hoc constraint but with our projection approach via
`python PRS-CS-proj/PRScs_proj.py --ref_dir=LD_REFERNCE_DIR--bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_proj --chrom=22 --write_pst=True --phi=1e-04 --thin=1 --n_iter=800 --n_burnin=1 --seed=1 --eigenval_rm=0.2`.
The resulting MCMC output will be stored it in the file `...`.
3. Plot the results via `Rscript plot-results/plot_figure1.R`.

### Synthetic Data Results
_Fill in_

### Real Data Results on Continuous Diseases
_Fill in_

### Real Data Results on Binary Diseases
_(I realize the step-by-step instructions I wrote up below probably aren't exactly correct, but the point is to give you an example of an actually reproducible workflow. 
In particular, users shouldn't have to read the comments inside your R scripts and modify/run different pieces of codes manually; that wouldn't be very reproducible.)_

1. Run `Rscript run-real-data-analysis/process_disease_dat.R`. 
This will generate a file  `data/sumdat_Rcov.txt` containing processed data, to be fed into the PRS pipelines.<br>
_(It looks like all the analyses use the same file name:
https://github.com/YuzhengDun1999/PRS-Bridge-article-code/blob/main/run-real-data-analysis/process_disease_dat.R#L12?
So you have to manually run a specific part of the script for each trait? 
This is such a bad design. 
Ideally, you should modify the code so that the datasets have different names.
If that takes too much work, at least add an input argument, like `run-real-data-analysis/process_disease_dat.R RA`, so that you can still run the whole pipeline more programmatically.)_
2. Run the quality control via `Rscript run-real-data-analysis/sumdat_QC.R`.
This will generate files `.../...`. 
3. Run the PRS methods on the preprocessed data via `Rscript run-real-data-analysis/${method_name}.R` where `method_name` can be chosen from `PRSBridge`, `PRScs`, `PRScs_proj`, `PRScs_regularized`, `PRScs_threshold`, and `ldpred2`.

**Note:** 
Under the `data` folder, we have provided a preprocessed sample dataset which can be used, without having to complete Step 1 and 2above, to run the version of PRS-CS without the ad hoc constraint on Chromosome 22 and generate the MCMC samples as plotted in Figure 1 of the manuscript.
<br>_(Note how much more details I put here compared to your "reproduce Figure 1." Your explanations tend to be way too incomplete.)_

### Plot results
Having run all the methods on all the datasets, Figure 4 and 5 can be reproduced via `Rscript plot-result/plot.R`.
Other figures in the main manuscript and the supplement can be reproduced by running other appropriate scripts under `plot-result/`.

 
## Directory Contents
This is an optional read for those interested in learning more about the specific files in the above steps to reproduce all the results.

### data
This folder contains the data used to generate posterior samples of the regression coefficients used to create Figure 1.
<br>_(Be more precise about what you mean by "data"; I guess `PRScs_sumdat.txt` contains the summary-level data and `chr22.bim` contains the individual-level data? These are important pieces of information.)_<br>
Also contained is an example summary-statistics file provided for illustration only.
<br>_(What does this mean? Illustration of what?)_<br>
Specifically, the files included are:
-   `PRScs_sumdat.txt`,  `chr22.bim`: For Figure 1. _(Be more precise.)_
-   `sumdat_Rcov.txt`: Example input summary statistics used by the preprocessing scripts in the run-real-data-analysis/ directory to demonstrate how input data should be formatted before preparing data for running all methods.

### PRS-CS-proj

This folder provides our modification of the PRS-CS code, the original version of which is found at [https://github.com/getian107/PRScs](https://github.com/getian107/PRScs). 
The files included are:

-   `PRScs_noconverge.py`: Implements a modified version of PRS-CS that removes the ad hoc constraint on the prior variance. 
See `plot-results/plot_fig1.R` for example usage.
<br>_(There doesn't seem to be any good reason here to bury the example usage inside the script. 
You can be more explicit and reader-friendly by providing the usage here; e.g. "To obtain the MCMC output shown in Figure 1, you can run `python PRS-CS-proj/PRScs_noconverge.py --ref_dir=LD_REFERNCE_DIR --bim_prefix=data/chr22  --sst_file=data/PRScs_sumdat.txt --n_gwas=284389 --out_dir=data/coef_noconverge --chrom=22 --write_pst=True --phi=1e-04 --thin=1 --n_iter=800 --n_burnin=1 --seed=1`.")_
- `PRScs_threshold.py`: Implements versions of PRS-CS modified to use different values of the upper threshold on the prior variance. 
See the scripts files generated by `run-real-data-analysis/PRScs_threshold.R` for example usage.
<br>_(I don't understand why you need to have a user run this R scripts to generate the bash scripts.
Why can't you just directly provide the bash scripts generated by this R script and place them in a folder `slurm_scripts/`, for example?
You can then just provide a bash command, e.g. `sbatch ...`, to submit the simulations to SLURM.
That way, it is far more transparent what steps are involved in generating the results.)_
-   `PRScs_proj.py`: Implements versions of PRS-CS modified to use projected summary statistics. 
See the scripts files generated by `run-real-data-analysis/PRScs_proj.R` for example usage.
-   `PRScs_Regularized.py`: Implements versions of PRS-CS modified to use Regularized LD matrix. 
See the scripts files generated by `run-real-data-analysis/PRScs_Regularized.R` for example usage.

### run-synthetic-data-analysis

This folder contains scripts to reproduce the numerical results of Section 4.1 based on the synthetic data.
The steps to reproduce the results are as follows:

1.  **Download synthetic data**  (see above)
2.  **Preprocess summary statistics:** Run `process_synthetic_dat.R`; detailed explanations of the input parameters are provided within the script.
3.  **Run PRS methods:**  Run `PRSBridge.R`, `PRScs.R`, and `ldpred2.R`; the input parameters are the same as those in `process_synthetic_dat.R`. Please see each method's tutorial on how to estimate or download LD reference panel.
4.  **Evaluate results:** Run `evaluation.R`; input parameters are the same as `process_synthetic_dat.R`.
5.  **Visualize:**  Run `plot.R` to produce Figure 3.

### run-real-data-analysis

This folder contains scripts to run each method on the real-data summary statistics. 
The computational demands make it unrealistic to reproduce the full results on a local machine, especially for PRS-Bridge and PRS-CS;
therefore, the provided R scripts instead generate a set of .sh files that can be submitted to a high-performance computing cluster (HPC) running on SLURM. 
The commands for running each method are specified within the generated .sh files. 
For LASSOSUM and ldpred2, we also provide scripts for directly running the methods, though we recommend similarly executing them on an HPC system.

**Step-by-step workflow:**

1.  **Preprocess GWAS summary statistics:**
    -   Disease traits: run  `process_disease_dat.R`
    -   Continuous traits (UK Biobank): run  `process_continuous_dat.R` to get summary statistics
2.  **Quality control:**  run  `sumdat_QC.R`
3.  **Run each PRS method:**  
    Use method-specific scripts (e.g.,  `PRScs.R`,  `PRSBridge.R`, etc.), ensuring the correct LD reference panel is specified (see each method's tutorial on how to estimate or download LD reference panel).

### evaluate-real-data-analysis
This folder contains scripts for evaluating the methods' performances:

-   `RE_Lassosum.R`: Computes relative efficiency for each method.
-   `get_sd_*.R`: Computes $R^2$ (or transformed AUC) and standard errors for each method:
    -   `get_sd_PRSBridge.R`
    -   `get_sd_PRScs.R`
    - ...

### plot-results

This folder contains scripts for visualization of results and reproducing figures:

-   `plot.R`: Generates results for Figures 4 and 5.
-   `plot_fig1.R`: Reproduces Figure 1 (posterior coefficients).
-   All plot scripts are annotated to specify which figure they reproduce and what input parameters to use in scripts in run-real-data-analysis section.
