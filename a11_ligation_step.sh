#!/bin/bash
#SBATCH -J ligation_step
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH -A ACF-UTK0171
#SBATCH --qos=condo
#SBATCH --partition=condo-trowan1
#SBATCH -t 00:45:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mhennige@vols.utk.edu 

singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools concat \
-n \
-f imputation/imputed_phased_snp50_file_names.txt \
-Ob \
-o imputation/ligated/ligated_imputed_phased_snp50_testing_animals_seq.chr25.bcf