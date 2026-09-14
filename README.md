# "Hey AlphaGenome, is the whole greater than the sum of its parts?"

Code to reproduce our blogpost testing whether AlphaGenome's prediction for an individual's
whole combination of HBB-window variants differs from the naive sum of each variant's own
single-variant Atlas effect.

## Installation

```shell
conda env create --file environment.yaml
conda activate wigttsoip
```

Also requires:
- An AlphaGenome Atlas API key, exported as `ALPHAGENOME_API_KEY` (used by `01-obtain_data`'s Atlas download rules).
- Local AlphaGenome JAX weights and GPU access for `02-preprocess_data`'s personalized-sequence forward passes (via [`alphagenome_research`](https://github.com/google-deepmind/alphagenome_research)).

## Structure

Three Snakemake stages, each its own `config/config.yaml`-driven pipeline:

1. **`workflows/01-obtain_data`** - downloads the 1000 Genomes 30x high-coverage chr11 VCF and
   the AlphaGenome JAX weights, and queries the AlphaGenome Atlas for single-variant RNA_SEQ
   effect scores across the HBB window (GTEx whole blood + all genes in the window).
2. **`workflows/02-preprocess_data`** - for each 1000 Genomes individual, builds their
   personalized HBB-window sequence (every SNV they carry, unphased) and scores the predicted
   RNA_SEQ change per gene/track with AlphaGenome's own gene-mask scorer, so results are
   directly comparable to the Atlas single-variant scores from step 1. Also runs one
   reference-only forward pass to get each gene's baseline expression level.
3. **`workflows/03-analysis`** - joins each individual's combined-variant effect against the
   sum of their own variants' single-variant effects, and computes per-gene top-k
   (strongest-variants-first) correlation curves. `notebooks/hbb.ipynb` (rendered via
   `figures.smk`) turns this into the blogpost's figures: variant burden, combined-vs-summed
   effect distributions and correlation (all variants and top-3), correlation-vs-k for HBB and
   for every gene in the window, and the same curves faceted by each gene's baseline expression.

## Running

Dev-scale rules (small sample/window) are always in `rule all`; long-running, GPU- or
quota-bound full-scale rules are left out and must be triggered explicitly by requesting their
output path (see the NOTE comment at the bottom of each stage's `Snakefile`).

```shell
sbatch src/scripts/submit_snakemake_slurm.sh 'snakemake -s workflows/01-obtain_data/Snakefile --use-conda --cluster "sbatch --parsable --job-name=snakejob --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.runtime} --partition={resources.partition} --qos={resources.qos} --gres={resources.gres} --output=slurm-%j.out --error=slurm-%j.err" --jobs 5 --latency-wait 30'

# same pattern for 02-preprocess_data and 03-analysis

# full-scale: trigger an out-of-rule-all target explicitly, e.g.
sbatch src/scripts/submit_snakemake_slurm.sh 'snakemake -s workflows/01-obtain_data/Snakefile --use-conda --cluster "..." --jobs 5 --latency-wait 30 data/raw/alphagenome-atlas/hbb_window/chunks/.done_full'
```