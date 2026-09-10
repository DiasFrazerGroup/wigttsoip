"""Read AlphaGenome Atlas variant effect scores for a gene/region/track subset from a chunk cache.

Companion to download_alphagenome_atlas_genexpr.py's --interval mode: that script never
merges an interval's chunks into one big file (for a wide interval, that single file would
be hundreds of GB). Instead it leaves the chunks as individually queryable per-region h5ad
files under --chunk-cache-dir, one file per grid chunk. This script reads only the chunk(s)
overlapping the requested region, and within each chunk, reads only the row/column slice of
the underlying X/quantiles matrices that actually match the requested genes/tracks - via
h5py directly, never loading a whole chunk's (obs x tracks) matrix into memory.
"""

import argparse
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from download_alphagenome_atlas_genexpr import DEFAULT_CHUNK_SIZE, grid_chunks, load_gencode_gene_ids, parse_interval


def _decode(arr):
    """h5py returns fixed/variable-length strings as bytes - decode to str."""
    return np.array([v.decode() if isinstance(v, bytes) else v for v in arr])


def _read_categorical(group, col):
    cat_group = group[col]
    categories = _decode(cat_group["categories"][:])
    codes = cat_group["codes"][:]
    return categories, codes


def _read_column(group, col):
    """Read a var/obs column as a plain string array, regardless of whether this
    particular chunk stored it as a pandas-categorical group (categories + codes) or a
    plain string dataset - anndata's write_h5ad picks one or the other per-chunk based
    on the column's actual cardinality in that chunk, so the same named column (e.g.
    'name', 'strand') can be either, inconsistently across chunks."""
    node = group[col]
    if isinstance(node, h5py.Group):
        categories, codes = _read_categorical(group, col)
        return categories[codes]
    return _decode(node[:])


def chunks_overlapping(chunk_cache_dir, interval, chunk_size):
    """Cached chunk files (that actually exist) overlapping `interval`."""
    paths = []
    for chunk in grid_chunks(interval, chunk_size):
        path = chunk_cache_dir / f"chunk_{chunk.chromosome}_{chunk.start}_{chunk.end}.h5ad"
        if path.exists():
            paths.append(path)
    return paths


def extract_chunk_long(path, interval, versioned_gene_ids=None, biosample_names=None):
    """Read only the matching row/column slice of one chunk's h5ad, as a long-format DataFrame.

    Returns None if nothing in this chunk matches the filters.
    """
    with h5py.File(path, "r") as f:
        gene_id_all = _read_column(f["obs"], "gene_id")
        if versioned_gene_ids is not None:
            gene_mask = np.isin(gene_id_all, list(versioned_gene_ids))
        else:
            gene_mask = np.ones(len(gene_id_all), dtype=bool)

        variant_all = _read_column(f["obs"], "variant")
        # Parse positions only once per *unique* variant string, not once per row - in
        # this gene-dense region each variant repeats once per overlapping gene, so
        # row count can be far larger than the number of distinct variants.
        unique_variants, inverse = np.unique(variant_all, return_inverse=True)
        positions_per_unique = np.array([int(v.split(":")[1]) for v in unique_variants])
        position_mask = ((positions_per_unique >= interval.start) & (positions_per_unique < interval.end))[inverse]

        row_mask = gene_mask & position_mask
        row_idx = np.where(row_mask)[0]
        if len(row_idx) == 0:
            return None

        biosample_name_all = _read_column(f["var"], "biosample_name")
        if biosample_names is not None:
            col_mask = np.isin(biosample_name_all, list(biosample_names))
        else:
            col_mask = np.ones(len(biosample_name_all), dtype=bool)
        col_idx = np.where(col_mask)[0]
        if len(col_idx) == 0:
            return None

        # Fancy-index one axis at a time (h5py doesn't support both at once): rows first
        # (the expensive axis to skip - only ~1/n_genes of n_obs), then columns in memory.
        scores = f["X"][row_idx, :][:, col_idx]
        quantiles = f["layers/quantiles"][row_idx, :][:, col_idx]

        gene_name_all = _read_column(f["obs"], "gene_name")
        obs = pd.DataFrame({
            "variant": variant_all[row_idx],
            "gene_id": gene_id_all[row_idx],
            "gene_name": gene_name_all[row_idx],
        })

        # "name" + "strand" are the library's own unique-track key (see
        # alphagenome.models.variant_scorers.tidy_anndata's docstring: "track_name and
        # track_strand" - the combined-effect side of the analysis is keyed the same
        # way). assay_title alone is NOT sufficient: e.g. K562 can have two distinct
        # "polyA plus RNA-seq" tracks (one strand='.', one strand='-') that are
        # different underlying tracks.
        var = pd.DataFrame({
            "track_id": _decode(f["var/_index"][col_idx]),
            "track_name": _read_column(f["var"], "name")[col_idx],
            "track_strand": _read_column(f["var"], "strand")[col_idx],
            "biosample_name": biosample_name_all[col_idx],
            "assay_title": _read_column(f["var"], "Assay title")[col_idx],
            "data_source": _read_column(f["var"], "data_source")[col_idx],
            "gtex_tissue": _read_column(f["var"], "gtex_tissue")[col_idx],
        })

    n_rows, n_cols = scores.shape
    long = pd.concat(
        [obs.loc[obs.index.repeat(n_cols)].reset_index(drop=True)]
        + [pd.concat([var] * n_rows, ignore_index=True)],
        axis=1,
    )
    long["score"] = scores.reshape(-1)
    long["quantile"] = quantiles.reshape(-1)
    return long


def _extract_chunk_long_star(args):
    path, interval, versioned_gene_ids, biosample_names = args
    return extract_chunk_long(path, interval, versioned_gene_ids=versioned_gene_ids, biosample_names=biosample_names)


def extract_to_batches(
    chunk_cache_dir,
    interval,
    output_dir,
    chunk_size=DEFAULT_CHUNK_SIZE,
    gene_ids=None,
    gene_annotation=None,
    biosample_names=None,
    max_workers=1,
    chunks_per_batch=50,
):
    """Write (variant, gene, track, score, quantile) rows for `interval` as a directory of
    zstd-compressed batch_*.parquet files, reading only the chunks and row/column slices
    that match `gene_ids`/`biosample_names`.

    Each source chunk is a separate h5ad file with its own fixed per-file overhead (open,
    decode obs/var string columns, np.unique on variant strings) - for a wide interval (e.g.
    the full 1,048,576bp window at chunk_size_bp=512, ~2048 chunk files) that per-chunk
    overhead dominates wall time when done serially with zero progress output, so this is
    parallelized across processes (h5py decoding is CPU-bound, not I/O-bound enough for
    threads to help) and prints progress periodically.

    Batches (not a single concatenated DataFrame) so peak memory stays bounded by
    `chunks_per_batch` chunks' worth of rows rather than the whole interval's - a returned
    per-chunk DataFrame is flushed to disk and dropped as soon as `chunks_per_batch` of them
    have accumulated, instead of being held (and IPC-pickled back from worker processes) for
    every one of ~2048 chunks before one final concat.
    """
    chunk_cache_dir = Path(chunk_cache_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    chunk_files = chunks_overlapping(chunk_cache_dir, interval, chunk_size)
    if not chunk_files:
        raise ValueError(f"No cached chunks overlap {interval} in {chunk_cache_dir}")

    versioned_gene_ids = None
    if gene_ids is not None:
        if gene_annotation is None:
            raise ValueError("gene_annotation is required to resolve gene_ids")
        versioned_gene_ids, _ = load_gencode_gene_ids(gene_annotation, gene_ids)

    print(f"Extracting {len(chunk_files)} chunk files with {max_workers} worker(s)...", flush=True)
    jobs = [(f, interval, versioned_gene_ids, biosample_names) for f in chunk_files]

    def write_batch(buf, batch_idx):
        if not buf:
            return
        batch_df = pd.concat(buf, ignore_index=True)
        tmp = output_dir / f"batch_{batch_idx:06d}.parquet.tmp"
        final = output_dir / f"batch_{batch_idx:06d}.parquet"
        batch_df.to_parquet(tmp, compression="zstd", index=False)
        tmp.rename(final)
        print(f"Wrote {final} ({len(batch_df)} rows)", flush=True)

    buf = []
    batch_idx = 0
    n_total_rows = 0
    t0 = time.time()

    def handle_part(part):
        nonlocal buf, batch_idx, n_total_rows
        if part is not None:
            buf.append(part)
            n_total_rows += len(part)
        if len(buf) >= chunks_per_batch:
            write_batch(buf, batch_idx)
            batch_idx += 1
            buf = []

    if max_workers <= 1:
        for i, job in enumerate(jobs, 1):
            handle_part(_extract_chunk_long_star(job))
            if i % 100 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} chunks read ({time.time() - t0:.0f}s elapsed)", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            for i, part in enumerate(ex.map(_extract_chunk_long_star, jobs, chunksize=4), 1):
                handle_part(part)
                if i % 100 == 0 or i == len(jobs):
                    print(f"  {i}/{len(jobs)} chunks read ({time.time() - t0:.0f}s elapsed)", flush=True)

    write_batch(buf, batch_idx)  # flush the remainder

    if n_total_rows == 0:
        raise ValueError(f"No rows matched gene_ids={gene_ids}, biosample_names={biosample_names} in {interval}")
    print(f"Extracted {n_total_rows} (variant, gene, track) rows total", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read AlphaGenome Atlas variant effect scores for a gene/region/track subset from a chunk cache, as a long-format parquet."
    )
    parser.add_argument(
        "--interval",
        required=True,
        help="Region to extract, as 'chrom:start-end' (0-based, half-open)",
    )
    parser.add_argument(
        "--chunk-cache-dir",
        required=True,
        help="Directory of cached chunks written by download_alphagenome_atlas_genexpr.py --interval",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Must match the chunk_size used when downloading",
    )
    parser.add_argument(
        "--gene-ids",
        nargs="+",
        default=None,
        help="Bare Ensembl gene IDs to restrict to (requires --gene-annotation)",
    )
    parser.add_argument(
        "--gene-annotation",
        default=None,
        help="GENCODE gtf parquet matching Atlas's gene annotation version (required with --gene-ids)",
    )
    parser.add_argument(
        "--biosample-name",
        nargs="+",
        default=None,
        help="Restrict tracks to these var['biosample_name'] values (e.g. K562)",
    )
    parser.add_argument("--output", required=True, help="Output directory of batch_*.parquet (zstd-compressed)")
    parser.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="Parallel processes for reading chunk files (useful for wide intervals with many chunks)",
    )
    parser.add_argument(
        "--chunks-per-batch",
        type=int,
        default=50,
        help="Source h5ad chunk files accumulated per output batch_*.parquet (bounds peak memory)",
    )

    args = parser.parse_args()

    extract_to_batches(
        args.chunk_cache_dir,
        parse_interval(args.interval),
        args.output,
        chunk_size=args.chunk_size,
        gene_ids=args.gene_ids,
        gene_annotation=args.gene_annotation,
        biosample_names=args.biosample_name,
        max_workers=args.max_workers,
        chunks_per_batch=args.chunks_per_batch,
    )
