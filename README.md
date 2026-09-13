# "Hey AlphaGenome, is the whole greater than the sum of its parts?"

This is the code to reproduce our blogpost on whether AlphaGenome makes predictions for
combinations of variants that are different than the sum of the effects of
different variants predicted alone.

## Installation

```shell
conda env create --file environment.yaml
conda activate wigttsoip
```

Also requires:
- An AlphaGenome Atlas API key, exported as `ALPHAGENOME_API_KEY` (used by `01-obtain_data`'s
  Atlas download rules).
- Local AlphaGenome JAX weights and GPU access for `02-preprocess_data`'s personalized-sequence
  forward passes (via [`alphagenome_research`](https://github.com/google-deepmind/alphagenome_research),
  submitted to this cluster's `gpu_diasfrazer`/`gpu` SLURM partitions - not the remote Atlas API).

## Structure

### 1. Obtain data (`workflows/01-obtain_data`)
- Download the 1000 Genomes 30x high-coverage chr11 VCF (`20201028_3202_raw_GT_with_annot`),
  scoped to the HBB window.
- Download the AlphaGenome JAX model weights (for the local GPU forward passes in step 2).
- Query the AlphaGenome Atlas API for single-variant RNA_SEQ effect scores of every SNV in
  the HBB window (dev: TSS +/- 5kb; full: the full 1,048,576bp window), and extract an
  HBB/whole blood (GTEx) long-format parquet from the cached chunks.

### 2. Preprocess data (`workflows/02-preprocess_data`)
- For each 1000 Genomes individual, build their personalized HBB-window sequence by applying
  every SNV they carry (unphased - no maternal/paternal split, indels/SVs skipped).
- Run both the reference and personalized sequences through the local AlphaGenome model
  (GPU), and score the predicted RNA_SEQ change per gene per track using AlphaGenome's own
  official gene-mask scorer (`GeneMaskLFCScorer`), so results are directly comparable to the
  Atlas single-variant scores from step 1.
- Subset the personalized-effect batches down to HBB/whole blood (GTEx) rows into one parquet.

### 3. Analysis (`workflows/03-analysis`)
- Get the distinct SNVs actually carried by the 1000 Genomes samples in scope (independent
  of any track/score data), then restrict Atlas's single-variant scores to just that set -
  Atlas scores every possible alt allele at every position, but only a fraction is ever
  polymorphic in a given sample set, so the join below runs against a much smaller table.
- Join, per sample and track (`track_name` alone - see note below), each individual's
  combination-of-variants effect (step 2, `whole_score`) against the sum of their variants'
  own single-variant effects (step 1, `sum_parts_score`) into a single long-format parquet -
  the additive/no-epistasis null against which to test for synergy or saturation.
- `notebooks/hbb.ipynb` (`figures.smk`, executed in place via `jupyter nbconvert`): a minimal
  four-panel story - (A) whole vs. sum of independent effects, (B) compression vs. summed
  effect magnitude, (C) effective dimensionality (strongest single variant / top-k curve),
  (D) residual context-dependence after removing compression - plus a weak-tail
  accumulation control and a supplementary section (GMM clustering, extreme individuals).

**Track-matching note**: the join matches on `track_name` alone, not also `track_strand`, as
a defensive default. For the GTEx whole blood polyA track used here, both the
personalized-sequence model (step 2) and Atlas's singles table (step 1) agree on
`strand='.'`, so this isn't currently working around any mismatch - kept anyway since
Atlas's singles table never has more than one `track_strand` value per `track_name`, so
matching on `track_name` alone can't create duplicate-match fan-out.

## Running

Each workflow stage is its own Snakefile with its own `config/config.yaml`-driven `rule all`.
Dev-scale rules (small sample/window sizes, fast to iterate on) are always included in
`rule all`; long-running, quota- or GPU-time-bound full-scale rules are deliberately left out
of `rule all` and must be triggered explicitly by requesting their output path (see the NOTE
comment at the bottom of each stage's `Snakefile`).

```shell
# dev-scale, any stage - iterate quickly (small sample list, small window)
sbatch src/scripts/submit_snakemake_slurm.sh 'snakemake -s workflows/01-obtain_data/Snakefile --use-conda --cluster "sbatch --parsable --job-name=snakejob --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.runtime} --partition={resources.partition} --qos={resources.qos} --gres={resources.gres} --output=slurm-%j.out --error=slurm-%j.err" --jobs 5 --latency-wait 30'

# same pattern for the other two stages
sbatch src/scripts/submit_snakemake_slurm.sh 'snakemake -s workflows/02-preprocess_data/Snakefile --use-conda --cluster "..." --jobs 5 --latency-wait 30'
sbatch src/scripts/submit_snakemake_slurm.sh 'snakemake -s workflows/03-analysis/Snakefile --use-conda --cluster "..." --jobs 5 --latency-wait 30'

# full-scale: trigger an out-of-rule-all target explicitly, e.g.
sbatch src/scripts/submit_snakemake_slurm.sh 'snakemake -s workflows/01-obtain_data/Snakefile --use-conda --cluster "..." --jobs 5 --latency-wait 30 data/raw/alphagenome-atlas/hbb_window/chunks/.done_full'
```

`src/scripts/submit_snakemake_slurm.sh` submits an `eval` of whatever Snakemake command you
pass as `sbatch --parsable ...` cluster jobs; **always submit via `sbatch`/`salloc`, never run
Snakemake or any heavy script directly on the login node.**

### The full 1,048,576bp / 3,202-sample personalized run

`predict_alphagenome_genexpr_full` (step 2, full scale) is the single most expensive rule in
the pipeline (~50s/sample on an H100 for the full-window forward pass + scoring, ~44.5h for
all 3,202 samples). The default, reproducible way to run it is the same as any other
explicit target - `sbatch src/scripts/submit_snakemake_slurm.sh` wrapping a
`snakemake --cluster ...` call targeting its output path (see "Running" above).

For this analysis, it was instead split across two GPUs by hand (via two now-removed,
cluster-specific `sbatch` scripts calling `alphagenome_genexpr.py` directly against a
disjoint half of the sample list each, writing to their own `full/chunk_{0,1}/`
subdirectory) to roughly halve the wall-clock time, then stopped early once 1,842 of the
3,202 samples were done rather than waiting for the rest. Both chunks' output was reconciled
into the single flat `full/` directory `predict_alphagenome_genexpr_full` itself would have
produced:

```shell
python src/scripts/merge_genexpr_full_chunks.py \
    --output results/02-preprocess_data/alphagenome/hbb_window/genexpr_personalized/full \
    --chunk-dirs results/02-preprocess_data/alphagenome/hbb_window/genexpr_personalized/full/chunk_0 \
                 results/02-preprocess_data/alphagenome/hbb_window/genexpr_personalized/full/chunk_1 \
    --vcf data/raw/1000genomes/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.recalibrated_variants.vcf.gz
```

**This stays reproducible via plain Snakemake**: nothing about the resulting `full/`
directory is special-cased - it's the exact same `batch_*.parquet` shape
`predict_alphagenome_genexpr_full` produces on its own, so every downstream rule
(`subset_hbb_wholeblood_genexpr_personalized_full` onward) runs against it unmodified, and running
`predict_alphagenome_genexpr_full` itself later would simply compute whatever samples are
still missing. The one gotcha: since Snakemake didn't run the job that wrote into `full/`
itself, it doesn't know that output is "complete" and raises `IncompleteFilesException` on
the next real Snakemake invocation - after independently confirming the directory's content
is valid (e.g. every `batch_*.parquet` reads cleanly and the sample counts match what you
expect), clear the stale `.snakemake/incomplete/` marker so Snakemake trusts it again:

```shell
snakemake -s workflows/02-preprocess_data/Snakefile --cleanup-metadata \
    results/02-preprocess_data/alphagenome/hbb_window/genexpr_personalized/full
```

(if that reports "metadata was not present", the real fix is deleting the matching
base64-named file directly under `.snakemake/incomplete/` - decode each with
`echo <name> | base64 -d` to confirm which path it corresponds to before removing it).
