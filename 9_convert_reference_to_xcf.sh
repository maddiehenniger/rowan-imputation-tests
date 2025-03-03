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

