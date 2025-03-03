#!/bin/bash
#SBATCH -J run_impute5
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -A ACF-UTK0171
#SBATCH --qos=condo
#SBATCH --partition=condo-trowan1
#SBATCH -t 00:45:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mhennige@vols.utk.edu 

while IFS= read -r line; do
    chr=$(echo "$line" | awk '{print $2}')
    region=$(echo "$line" | awk '{print $4}')
    buffer=$(echo "$line" | awk '{print $3}')
    count=$(echo "$line" | awk '{print $1}')
    out_file="imputation/imputed_phased_snp50_testing_animals_seq.${chr}_${count}.bcf"
    log_file="stats/imputated_chunks_testing_animals_seq.${chr}_${count}_${region}.out"
    ../../z_rowan_imp_pipeline_tests/3_imputation/impute5/impute5_v1.2.0/impute5_v1.2.0_static \
    --h databases/MU_HD_only.chr25_xcf.bcf \
    --g phasing/phased_snp50_testing_animals_seq.chr25.bcf \
    --r ${region} \
    --buffer-region ${buffer} \
    --o ${out_file} \
    --l ${log_file}
done < imputation/phased_snp50_testing_animals_seq_chunked_coords.txt