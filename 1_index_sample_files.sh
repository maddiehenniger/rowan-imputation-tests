#!/bin/bash
#SBATCH -J index_sample_files
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -A ACF-UTK0171
#SBATCH --qos=condo
#SBATCH --partition=condo-trowan1
#SBATCH -t 00:05:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mhennige@vols.utk.edu 

singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools index \
--threads 8 \
--o samples/Testing_animals_seq.chr25.vcf.gz.csi \
samples/Testing_animals_seq.chr25.vcf.gz