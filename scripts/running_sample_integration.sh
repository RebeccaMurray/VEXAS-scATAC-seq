#!/bin/bash
#SBATCH --job-name=signac_merge_samples
#SBATCH --mail-type=all
#SBATCH --mail-user=rmurray@nygenome.org
#SBATCH --partition=pe2
#SBATCH --mem=150g
#SBATCH --output=logs/%x_%j.log
#SBATCH --cpus-per-task=20

module load gdal
module load gsl/2.6
module load zlib/1.2.11 
module load macs2
module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23_seurat_5

/gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23_seurat_5/bin/Rscript /gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/R/sample_integration/integrating_samples.R
