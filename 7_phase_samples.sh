#!/bin/bash
#SBATCH -J phase_samples
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH -A ACF-UTK0171
#SBATCH --qos=condo
#SBATCH --partition=condo-trowan1
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mhennige@vols.utk.edu 

singularity exec https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--hb60d31d_0 SHAPEIT5_phase_common \
--input samples/snp50_testing_animals_seq.chr25.bcf.gz \
--reference databases/MU_HD_only.chr25.bcf \
--thread 24 \
--region 25 \
--output phasing/phased_snp50_testing_animals_seq.chr25.bcf
