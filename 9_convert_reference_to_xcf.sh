#!/bin/bash
#SBATCH -J convert_ref_xcf
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -A ACF-UTK0171
#SBATCH --qos=condo
#SBATCH --partition=condo-trowan1
#SBATCH -t 00:45:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mhennige@vols.utk.edu 

../../z_rowan_imp_pipeline_tests/3_imputation/impute5/impute5_v1.2.0/xcftools_static view \
--input databases/MU_HD_only.chr25.bcf \
--output databases/MU_HD_only.chr25_xcf.bcf \
--format sh \
--region 25 \
--thread 8 \
--maf 0.03125 \
--log stats/xcf_log.out