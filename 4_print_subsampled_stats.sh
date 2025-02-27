#!/bin/bash
#SBATCH -J stats_subsampled_file
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -A ACF-UTK0171
#SBATCH --qos=condo
#SBATCH --partition=condo-trowan1
#SBATCH -t 00:05:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mhennige@vols.utk.edu 

singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools stats \
samples/snp50_testing_animals_seq.chr25.bcf.gz >> stats/snp50_testing_animals_seq.chr25.stats.out
