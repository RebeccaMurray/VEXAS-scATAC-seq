#!/bin/bash
#SBATCH --job-name=chromvar
#SBATCH --mail-type=all
#SBATCH --mail-user=rmurray@nygenome.org
#SBATCH --partition=dev
#SBATCH --mem=100g
#SBATCH --cpus-per-task=10
#SBATCH --output=logs/%x_%j.log

module load gdal
module load gsl/2.6
module load zlib/1.2.11 
module load macs2
module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23_seurat_5

/gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23_seurat_5/bin/Rscript \
/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/R/chromVAR/chromvar_analysis_with_recalling_peaks.R $1
