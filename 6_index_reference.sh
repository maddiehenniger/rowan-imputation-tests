#!/bin/bash
#SBATCH -J index_reference_db
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
--o databases/MU_HD_only.chr25.bcf.csi \
databases/MU_HD_only.chr25.bcf