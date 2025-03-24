#!/usr/bin/env bash
#SBATCH -c 10
#SBATCH --mem=20GB
#SBATCH -N 1
#SBATCH --job-name network
#SBATCH --output=logs/NG-%a.%A.out # STDOUT file, %j for job id, %N for hostname
#SBATCH --error=logs/NG-%a.%A.err  # STDERR fil
#SBATCH --partition=medium
#SBATCH --time=24:00:00          # total run time limit in HH:MM:SS


### To activate the conda environment

#Before running this script activate your conda environment using : "conda activate R_env"

### Run your Job

fastspar --yes --otu_table df_otu.tsv --correlation median_correlation.tsv --covariance median_covariance.tsv
fastspar_bootstrap --otu_table df_otu.tsv --number 1000 --prefix bootstrap_counts/df_otu
parallel fastspar --yes --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
fastspar_pvalues --otu_table df_otu.tsv --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_df_otu --permutations 1000 --outfile pvalues.tsv
Rscript Network_fastspar_Rscript.R
