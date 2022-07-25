#!/bin/bash
#SBATCH --time=11:50:59
#SBATCH --output=16_Run_iFibroblasts_COPD.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100GB
#SBATCH --gpus-per-node=1

set -e

module load R/4.1.0-gnu9.1

cd /fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Codes
Rscript 16_iFibroblasts_COPD_comparison.R
