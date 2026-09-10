#!/usr/bin/env bash
#SBATCH --job-name=genexpr_chunk0
#SBATCH --no-requeue
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000M
#SBATCH --gres=gpu:h100:1
#SBATCH --partition=gpu_diasfrazer
#SBATCH --output=slurm-genexpr-chunk0-%j.out
#SBATCH --error=slurm-genexpr-chunk0-%j.err

set -eo pipefail

cd /users/diasfrazer/manglada/projects/wigttsoip
source /users/diasfrazer/manglada/miniforge3/etc/profile.d/conda.sh
conda activate wigttsoip

SAMPLES=$(cat data/prep/alphagenome-genexpr/full_chunks/chunk_0_samples.txt)

python workflows/02-preprocess_data/scripts/alphagenome_genexpr.py \
    --vcf data/raw/1000genomes/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.recalibrated_variants.vcf.gz \
    --fasta data/raw/GENCODE/release_46/GRCh38.primary_assembly.genome.fa.gz \
    --gtf data/raw/GENCODE/release_46/gencode.v46.annotation.gtf.parquet \
    --weights-dir data/raw/articles/Avsec2026/alphagenome_jax/all-folds \
    --gene-ids-file data/prep/alphagenome-atlas/hbb_window/gene_ids.txt \
    --samples $SAMPLES \
    --chromosome chr11 \
    --tss 5229394 \
    --window-size 1048576 \
    --batch-size 2 \
    --output results/02-preprocess_data/alphagenome/hbb_window/genexpr_personalized/full/chunk_0

echo "Done!"
