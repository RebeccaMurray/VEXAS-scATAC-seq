#!/bin/bash
#SBATCH --job-name=gene_scores
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmurray@nygenome.org
#SBATCH --partition=dev
#SBATCH --mem=150g
#SBATCH --output=logs/%x_%j.log
#SBATCH --cpus-per-task=5

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23_seurat_5

/gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23_seurat_5/bin/Rscript /gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/R/differential_gene_scores/VEXAS_MUT_vs_WT/running_lmm_for_gene_accessibility.R
