# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

See `README.md` for the research question, data flow, and full stage-by-stage description.
In short: a three-stage Snakemake pipeline (`01-obtain_data`, `02-preprocess_data`,
`03-analysis`) testing whether AlphaGenome predicts a 1000 Genomes individual's combination
of HBB-window variants to differ from the naive sum of each variant's own single-variant
Atlas effect.

## Running the workflows

**Never run Snakemake, `salloc`, or any pipeline/model code directly on the login node** -
always submit through `sbatch`. The standard pattern, matching every rule's own per-job
resource spec via `--cluster`:

```bash
sbatch src/scripts/submit_snakemake_slurm.sh 'snakemake -s workflows/<N>-<name>/Snakefile --use-conda --cluster "sbatch --parsable --job-name=snakejob --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.runtime} --partition={resources.partition} --qos={resources.qos} --gres={resources.gres} --output=slurm-%j.out --error=slurm-%j.err" --jobs 5 --latency-wait 30 [target]'
```

`submit_snakemake_slurm.sh` submits the Snakemake *driver* itself as one long-running
`marathon`-qos job; the driver then submits each rule as its own real `sbatch` job sized to
that rule's own `resources:` block. Without `--cluster`, `snakemake -j N` runs every rule
locally inside whatever single allocation you happened to request - fine for quick dev-scale
iteration inside a short `salloc`, but it silently ignores each rule's own `resources.mem_mb`
(all rules share the one outer allocation), which is exactly how a real full-scale run got
OOM-killed on a step whose own rule asked for far more memory than the ad hoc outer job had.

`snakemake` itself lives in the `base` conda env, not this project's `wigttsoip` env (only
individual rules, via `conda: "wigttsoip"` + `--use-conda`, run inside that env) - if you
need to invoke it directly rather than through `submit_snakemake_slurm.sh`, use the full
path (`/users/diasfrazer/manglada/miniforge3/bin/snakemake`) rather than activating
`wigttsoip` first, or you'll get `snakemake: command not found`.

## Cluster operational notes (lessons learned the hard way)

- **GPU gres alias quirk**: a `sbatch: HTTP/1.1 400 Bad Request: ... "Unknown GPU type: 'h100'"` message on submission is benign - the job still submits and runs correctly. Some
  monitoring/accounting hook on this cluster doesn't recognize the short GPU-type alias
  (`h100`) that `slurm.conf`'s own GresTypes actually uses; it doesn't block scheduling.
- **Never restart per-sample batch numbering at 0 in a resumable script.** `alphagenome_genexpr.py`'s `already_processed_samples()` correctly scans existing
  `batch_*.parquet` files to decide which samples still need scoring - but on every fresh
  invocation it used to number new output files starting at `batch_000000` again regardless
  of what already existed, silently overwriting (destroying) earlier completed batches
  instead of adding to them. Confirmed via cross-referencing job logs: **219 batches
  destroyed and rewritten across a handful of resubmissions** (~438 sample-scorings' worth
  of wasted GPU compute), with the final on-disk sample count barely changing between runs
  as a symptom. Fixed by always continuing numbering from `next_batch_index()` (highest
  existing index + 1) - any similarly "resumable, writes numbered output files" script must
  do the same, mirroring `merge_genexpr_full_chunks.py`'s own `next_batch_index()` pattern.
- **duckdb: never carry a large/duplicated string column through an `UNNEST` explode.**
  `summarize_whole_vs_sum_parts.py`'s join used to `SELECT` the *original* semicolon-joined
  `variant_id` string (up to ~3,800 SNVs per sample in the full 1,048,576bp window - tens of
  KB of text) alongside an `UNNEST(string_split(variant_id, ';'))` column. DuckDB duplicates
  every other selected column across each row `UNNEST` produces, so that whole string got
  copied into every one of a sample's thousands of exploded rows - tens of GB of pure
  redundant duplication, not the join itself, is what OOM-killed the job even at 32GB.
  Fixed by exploding through a lightweight surrogate `ROW_NUMBER()` key and reattaching the
  original string once at the very end via a final join back to the small source table.
  **Diagnose this class of bug by asking "does a wide/large column ride along through the
  fan-out?" before reaching for more memory.**
- **duckdb `SET memory_limit` under SLURM needs real headroom, not an exact match.** Without
  an explicit `memory_limit`, duckdb sizes its buffer pool off the node's total physical RAM,
  not the job's actual cgroup `--mem` allocation - it over-allocates and gets OS OOM-killed
  rather than spilling to disk gracefully. Setting `memory_limit` to *exactly* the job's
  `--mem` isn't safe either (confirmed: still OOM-killed at 32000MB with `memory_limit =
  '32000MB'`, MaxRSS landing right at the cgroup boundary) - duckdb's tracked buffer pool
  doesn't account for the Python process itself or duckdb's own overhead outside that pool.
  `summarize_whole_vs_sum_parts.py`'s `--memory-limit-mb` is wired to match `resources.mem_mb`
  as a convention (so both numbers come from one place in the rule), but the real fix for an
  OOM is almost always removing the underlying memory blowup (see above), not just raising
  the number.
- **Stale `.snakemake/incomplete/` markers after writing a Snakemake-declared output via an
  out-of-DAG script.** Reconciling manually-split GPU work (see `README.md`'s "full 1MB /
  3,202-sample personalized run" section) writes into a directory Snakemake itself didn't
  create, so the next real Snakemake run raises `IncompleteFilesException` even though the
  content is fine. `snakemake --cleanup-metadata <path>` is the documented fix, but if it
  reports "metadata was not present" the real marker lives directly under
  `.snakemake/incomplete/` as a base64-encoded filename - decode each with `echo <name> |
  base64 -d` to confirm which path it corresponds to before deleting it directly.
- **Stale directory-wide `.snakemake/locks`.** All Snakefiles sharing one project root share
  one `.snakemake/` lock - a driver job that's still alive (even stuck forever polling a
  cluster sub-job that already died) blocks *every* Snakemake invocation against this repo,
  not just its own Snakefile. `ps aux`/`sacct` to confirm nothing genuine is running, cancel
  the stuck driver, then `snakemake --unlock`.

### Snakemake gotcha: never make an in-place-executed file the rule's declared `output:`

Snakemake deletes a rule's declared `output:` files right before running its shell command
(so it can detect failed jobs reliably). `render_hbb_notebook` (`workflows/03-analysis/rules/figures.smk`) executes `notebooks/hbb.ipynb` in place via `jupyter nbconvert --inplace`
- when the notebook itself was the declared `output:`, Snakemake deleted it moments before
nbconvert tried to open it (`pattern 'notebooks/hbb.ipynb' matched no files`, and the file
had to be reconstructed from scratch). Fixed the same way `alphagenome_finetuning_rna`'s
CLAUDE.md documents for a resumable training checkpoint: declare a separate completion
marker (`output: touch("notebooks/.hbb_rendered")`) as the real Snakemake output, and list
the actual file the shell command reads/writes in place as an `input:` instead - its content
changing still triggers a rerun, but Snakemake never deletes it first.

## Configuration

All paths/params are centralized in `config/config.yaml` (`gencode`, `alphagenome_atlas`,
`thousand_genomes`, `alphagenome_jax`, `alphagenome_genexpr`, `analysis`, `gnomad`). See
`README.md` for what each stage actually does with them.

## Conda environments

Single project env, `environment.yaml` -> `wigttsoip` (pandas, pyarrow, pyranges, anndata,
pysam, duckdb, the vendored/patched `alphagenome` package, `jax[cuda12]` +
`alphagenome_research` for local GPU forward passes, `ipykernel`/`nbconvert` for
`notebooks/hbb.ipynb`'s headless execution via `figures.smk`). `snakemake` itself is in the
`base` env - see "Running the workflows" above.
