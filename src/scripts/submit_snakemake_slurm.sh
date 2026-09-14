#!/usr/bin/env bash
#SBATCH --job-name=snakepipe
#SBATCH --no-requeue
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6000M
#SBATCH --gres=none
#SBATCH --partition=genoa64
#SBATCH --qos=marathon
#SBATCH --output=snakepipe.%x_%A.out
#SBATCH --error=snakepipe.%x_%A.err

set -eo pipefail

# unpack arguments
SNAKEMAKE_COMMAND=$1

# run command
eval $SNAKEMAKE_COMMAND

echo "Done!"
