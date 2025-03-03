#!/bin/bash
#SBATCH -J index_subsampled_samples
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH -A ACF-UTK0171
#SBATCH --qos=condo
#SBATCH --partition=condo-trowan1
#SBATCH -t 00:05:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mhennige@vols.utk.edu 

singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools index \
--threads 2 \
--o samples/snp50_testing_animals_seq.chr25.bcf.gz.csi \
samples/snp50_testing_animals_seq.chr25.bcf.gz
