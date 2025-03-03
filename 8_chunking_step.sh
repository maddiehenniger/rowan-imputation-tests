#!/bin/bash
#SBATCH -J chunk_step
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -A ACF-UTK0171
#SBATCH --qos=condo
#SBATCH --partition=condo-trowan1
#SBATCH -t 00:45:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mhennige@vols.utk.edu 

../../z_rowan_imp_pipeline_tests/3_imputation/impute5/impute5_v1.2.0/imp5Chunker_v1.2.0_static \
--h databases/MU_HD_only.chr25.bcf \
--g phasing/phased_snp50_testing_animals_seq.chr25.bcf \
--r 25 \
--l stats/phased_chunking.log \
--o imputation/phased_snp50_testing_animals_seq_chunked_coords.txt