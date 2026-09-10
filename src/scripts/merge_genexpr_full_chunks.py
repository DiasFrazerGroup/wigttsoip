"""One-off reconciliation after manually splitting predict_alphagenome_genexpr_full's
3,202 samples across two GPUs (gpu_diasfrazer/h100 + gpu/7g.80gb) via custom sbatch
scripts instead of Snakemake, to cut wall-clock time roughly in half.

Moves every batch_*.parquet out of full/chunk_0/ and full/chunk_1/ up into full/
itself, renumbering to continue after whatever batches were already there (from
before the chunked jobs started), then removes the now-empty chunk directories.
Afterwards, full/ is a flat directory of batch_*.parquet exactly as if
predict_alphagenome_genexpr_full had produced it directly - Snakemake's own
`output: directory(...)` check for that rule is then satisfied with no
Snakemake-visible difference from a normal single-job run.
"""

import argparse
import shutil
from pathlib import Path

import pandas as pd
import pysam


def next_batch_index(output_dir):
    existing = sorted(output_dir.glob("batch_*.parquet"))
    if not existing:
        return 0
    last = existing[-1].stem  # "batch_000005"
    return int(last.split("_")[1]) + 1


def merge(output_dir, chunk_dirs):
    output_dir = Path(output_dir)
    idx = next_batch_index(output_dir)
    n_moved = 0
    for chunk_dir in chunk_dirs:
        chunk_dir = Path(chunk_dir)
        for f in sorted(chunk_dir.glob("batch_*.parquet")):
            dest = output_dir / f"batch_{idx:06d}.parquet"
            shutil.move(str(f), str(dest))
            idx += 1
            n_moved += 1
        remaining = list(chunk_dir.glob("*"))
        if remaining:
            print(f"Warning: {chunk_dir} still has files, not removing: {remaining}")
        else:
            chunk_dir.rmdir()
    print(f"Moved {n_moved} batch files into {output_dir}, next index would be {idx}")


def verify(output_dir, vcf_path):
    with pysam.VariantFile(vcf_path) as f:
        all_samples = set(f.header.samples)
    done = set()
    for f in Path(output_dir).glob("batch_*.parquet"):
        done.update(pd.read_parquet(f, columns=["sample"])["sample"].unique())
    missing = all_samples - done
    print(f"{len(done)}/{len(all_samples)} samples present in {output_dir}")
    if missing:
        print(f"Missing {len(missing)} samples: {sorted(missing)[:10]}...")
    else:
        print("All samples accounted for.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Canonical flat output directory (e.g. .../full)")
    parser.add_argument("--chunk-dirs", nargs="+", required=True)
    parser.add_argument("--vcf", required=True, help="For sample-completeness verification")
    args = parser.parse_args()

    merge(args.output, args.chunk_dirs)
    verify(args.output, args.vcf)
