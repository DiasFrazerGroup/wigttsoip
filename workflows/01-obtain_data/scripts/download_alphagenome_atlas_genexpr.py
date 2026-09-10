"""Download precomputed AlphaGenome Atlas variant effect scores on gene expression."""

import argparse
import os
import time
from pathlib import Path

import anndata
import grpc
import pandas as pd
from alphagenome.atlas import atlas
from alphagenome.data import genome

MAX_RETRIES = 5
RETRY_BASE_DELAY = 30  # seconds; empirically the quota needs ~30-60s to reset, shorter delays exhaust all retries and crash the batch

DEFAULT_CHUNK_SIZE = 2_048  # a single 16_384bp chunk near HBB measured ~12.8GB RSS (every overlapping gene, unfiltered); keep this small to bound per-chunk memory


def parse_interval(interval_str):
    """Parse a 'chrom:start-end' string into a genome.Interval (0-based, half-open)."""
    chrom, span = interval_str.split(":")
    start, end = span.split("-")
    return genome.Interval(chrom, int(start), int(end))


def load_gencode_gene_ids(gene_annotation_parquet, gene_ids):
    """Map bare Ensembl gene IDs to GENCODE-versioned ones (e.g. ENSG00000100365 -> ENSG00000100365.16).

    Atlas's gene_id values are versioned to the GENCODE release it was trained
    on (v46, per google-deepmind/alphagenome's own visualization colab), while
    our GTFs (Ensembl) carry bare IDs. Genes absent from this annotation can
    never be matched in Atlas responses, so we report and drop them upfront
    rather than silently returning zero rows for them.
    """
    genes = pd.read_parquet(
        gene_annotation_parquet, columns=["Feature", "gene_id"], filters=[("Feature", "==", "gene")]
    )
    versioned_by_bare = dict(
        zip(genes["gene_id"].str.split(".").str[0], genes["gene_id"])
    )
    missing = [g for g in gene_ids if g not in versioned_by_bare]
    if missing:
        print(
            f"Warning: {len(missing)}/{len(gene_ids)} requested genes not found in "
            f"{gene_annotation_parquet}, dropping: {missing}",
            flush=True,
        )
    return [versioned_by_bare[g] for g in gene_ids if g in versioned_by_bare], versioned_by_bare


def load_variants(variants_parquets):
    """Load one or more variants parquets (id_var, chrom, pos, ref, alt, gene_id) into genome.Variant objects.

    Returns the variants alongside each one's *assigned* gene_id (the gene it
    was actually sampled for, e.g. in popping_gnomad) so callers can restrict
    the (usually much broader) Atlas response down to just those intended
    variant-gene pairs.
    """
    df = pd.concat(
        [
            pd.read_parquet(f, columns=["id_var", "chrom", "pos", "ref", "alt", "gene_id"])
            for f in variants_parquets
        ],
        ignore_index=True,
    )
    variants = [
        genome.Variant(row.chrom, int(row.pos), row.ref, row.alt)
        for row in df.itertuples(index=False)
    ]
    return variants, df["id_var"].tolist(), df["gene_id"].tolist()


def query_with_retries(query_fn, *args, **kwargs):
    """Call an Atlas query function, retrying on RESOURCE_EXHAUSTED (Atlas rate limits)."""
    for attempt in range(MAX_RETRIES):
        try:
            return query_fn(*args, **kwargs)
        except grpc.RpcError as error:
            if error.code() != grpc.StatusCode.RESOURCE_EXHAUSTED or attempt == MAX_RETRIES - 1:
                raise
            delay = RETRY_BASE_DELAY * (2**attempt)
            print(
                f"Atlas quota exceeded (attempt {attempt + 1}/{MAX_RETRIES}), "
                f"retrying in {delay}s...",
                flush=True,
            )
            time.sleep(delay)


def postprocess(
    adata,
    versioned_gene_ids=None,
    id_var_by_variant=None,
    assigned_gene_by_id_var=None,
):
    """Restrict to the wanted genes (all RNA_SEQ tracks kept); attach id_var if given.

    Atlas's gene_id filter param can't be used server-side because its values
    are GENCODE-versioned (e.g. "ENSG00000100365.16") and don't match our bare
    Ensembl IDs, and it has no per-variant gene filter anyway. So we query
    without a gene filter (server returns every gene overlapping each variant)
    and filter here instead, on the exact versioned ID resolved from the
    GENCODE annotation (see load_gencode_gene_ids). Two modes:
      - assigned_gene_by_id_var: keep only each variant's own sampled gene
        (e.g. from popping_gnomad) - the narrow, intended case.
      - versioned_gene_ids: keep any row whose gene is in this set - the
        broad case, used when there's no per-variant assignment (interval mode).
    """
    adata = adata.copy()
    adata.obs["variant"] = adata.obs["variant"].astype(str)
    if id_var_by_variant is not None:
        adata.obs["id_var"] = adata.obs["variant"].map(id_var_by_variant)
    if assigned_gene_by_id_var is not None:
        bare_gene_id = adata.obs["gene_id"].str.split(".").str[0]
        wanted_gene_id = adata.obs["id_var"].map(assigned_gene_by_id_var)
        adata = adata[bare_gene_id == wanted_gene_id].copy()
    elif versioned_gene_ids is not None:
        adata = adata[adata.obs["gene_id"].isin(versioned_gene_ids)].copy()
    return adata


def already_done_id_vars(temp_dir):
    """Read id_var values already queried, whether they matched their assigned gene or not.

    A variant whose assigned gene is in GENCODE but simply has no overlap in
    Atlas's own annotation also produces zero rows - same as one whose gene
    is missing from GENCODE entirely. Relying on the h5ad's obs["id_var"]
    (matched rows only) would mean such a variant is never marked done and
    gets re-queried on every resume forever. Each batch also writes a plain
    "attempted" id list alongside its h5ad, covering both cases.
    """
    done = set()
    for attempted_file in sorted(temp_dir.glob("batch_*.attempted.txt")):
        done.update(attempted_file.read_text().split())

    # Backward compatibility with batches written before the .attempted.txt
    # companion existed: fall back to their matched rows.
    for batch_file in sorted(temp_dir.glob("batch_*.h5ad")):
        if (batch_file.parent / f"{batch_file.stem}.attempted.txt").exists():
            continue
        try:
            batch_obs = anndata.read_h5ad(batch_file, backed="r").obs
            done.update(batch_obs["id_var"].unique())
        except Exception as e:
            print(f"Warning: deleting corrupted batch file {batch_file}: {e}", flush=True)
            batch_file.unlink(missing_ok=True)
    return done


def grid_chunks(interval, chunk_size):
    """Sub-intervals of `chunk_size` covering `interval`, snapped to a fixed global grid.

    Snapping chunk boundaries to absolute multiples of `chunk_size` (rather than
    relative to this particular interval's start) means two differently
    positioned/sized interval queries against the same chromosome share the
    exact same chunk files wherever they overlap - e.g. a small TSS+/-5kb dev
    window and a later full-window run both hit the same handful of cached
    chunks around the TSS, so expanding the window never re-downloads them.
    """
    grid_start = (interval.start // chunk_size) * chunk_size
    for chunk_start in range(grid_start, interval.end, chunk_size):
        yield genome.Interval(interval.chromosome, chunk_start, chunk_start + chunk_size)


def query_interval_in_chunks(
    client,
    interval,
    requested_scorer,
    max_workers,
    chunk_size,
    chunk_cache_dir,
):
    """Query a genomic interval in resumable, grid-aligned chunks.

    Each chunk is cached *without* gene filtering (only restricted to GTEx
    tracks) so that a later run requesting a different/larger set of genes
    over the same region can reuse it as-is - gene filtering happens once,
    at merge time, in `main`.
    """
    chunk_cache_dir.mkdir(parents=True, exist_ok=True)
    chunks = list(grid_chunks(interval, chunk_size))

    cached_files = []
    n_done = 0
    for i, chunk in enumerate(chunks):
        stem = f"chunk_{chunk.chromosome}_{chunk.start}_{chunk.end}"
        chunk_path = chunk_cache_dir / f"{stem}.h5ad"
        if chunk_path.exists():
            n_done += 1
            cached_files.append(chunk_path)
            continue

        print(
            f"[{i + 1}/{len(chunks)}] Querying chunk {chunk}...",
            flush=True,
        )
        result = query_with_retries(
            client.query_interval,
            chunk,
            requested_scorers=[requested_scorer],
            max_workers=max_workers,
            progress_bar=False,
        )
        adata = postprocess(result[requested_scorer])

        chunk_path_tmp = chunk_cache_dir / f"{stem}.h5ad.tmp"
        adata.write_h5ad(chunk_path_tmp, compression="gzip")
        chunk_path_tmp.rename(chunk_path)
        cached_files.append(chunk_path)

    print(f"{n_done}/{len(chunks)} chunks were already cached in {chunk_cache_dir}", flush=True)
    return cached_files


def query_variants_in_batches(
    client,
    variants,
    id_vars,
    assigned_gene_by_id_var,
    requested_scorer,
    max_workers,
    batch_size,
    temp_dir,
):
    """Query variants in small resumable batches, skipping already-completed ones."""
    temp_dir.mkdir(parents=True, exist_ok=True)

    done = already_done_id_vars(temp_dir)
    if done:
        print(f"Found {len(done)} already-downloaded variants in {temp_dir}, skipping them", flush=True)

    pending = [
        (variant, id_var)
        for variant, id_var in zip(variants, id_vars)
        if id_var not in done
    ]
    print(f"{len(pending)} variants left to query ({len(done)} already done)", flush=True)

    existing_batch_ids = [
        int(f.stem.split("_")[1]) for f in temp_dir.glob("batch_*.h5ad")
    ]
    next_batch_id = max(existing_batch_ids, default=-1) + 1

    for batch_start in range(0, len(pending), batch_size):
        batch = pending[batch_start : batch_start + batch_size]
        batch_variants = [v for v, _ in batch]
        batch_id_vars = [i for _, i in batch]

        print(
            f"Querying batch {next_batch_id} ({len(batch_variants)} variants)...",
            flush=True,
        )
        result = query_with_retries(
            client.query_variants,
            batch_variants,
            requested_scorers=[requested_scorer],
            max_workers=max_workers,
            progress_bar=False,
        )
        adata = postprocess(
            result[requested_scorer],
            id_var_by_variant=dict(zip((str(v) for v in batch_variants), batch_id_vars)),
            assigned_gene_by_id_var=assigned_gene_by_id_var,
        )

        batch_stem = f"batch_{next_batch_id:06d}"
        batch_path = temp_dir / f"{batch_stem}.h5ad"
        batch_path_tmp = temp_dir / f"{batch_stem}.h5ad.tmp"
        adata.write_h5ad(batch_path_tmp, compression="gzip")
        batch_path_tmp.rename(batch_path)

        # Record every attempted id_var, matched or not, so an unmatched
        # variant-gene pair (no overlap in Atlas's own annotation) is never
        # re-queried on a future resume - it's already a known non-match.
        attempted_path = temp_dir / f"{batch_stem}.attempted.txt"
        attempted_path_tmp = temp_dir / f"{batch_stem}.attempted.txt.tmp"
        attempted_path_tmp.write_text("\n".join(batch_id_vars))
        attempted_path_tmp.rename(attempted_path)

        next_batch_id += 1

    return sorted(temp_dir.glob("batch_*.h5ad"))


def main(
    output,
    api_key,
    gene_annotation,
    interval=None,
    variants_parquets=None,
    gene_ids=None,
    requested_scorer="RNA_SEQ",
    max_workers=1,
    batch_size=100,
    chunk_size=DEFAULT_CHUNK_SIZE,
    chunk_cache_dir=None,
):
    """Query AlphaGenome Atlas and save GTEx RNA-seq tracks for the requested genes.

    Exactly one of `interval` (dense per-position scan) or `variants_parquet`
    (only the given variant-gene pairs) must be provided. Both modes are
    resumable via small cached pieces kept in a temp/cache directory: a quota
    failure or requeue only forces re-querying pieces not yet cached, never
    ones already saved.

    Args:
        output: Path to output file. In --variants mode, an .h5ad. In
            --interval mode, just an empty marker file (touched once every
            chunk is cached) - the interval's chunks are kept as separate,
            individually queryable h5ad files in chunk_cache_dir rather than
            merged into one; see extract_alphagenome_atlas_genexpr.py to
            read a gene/region's data on demand.
        api_key: AlphaGenome Atlas API key.
        gene_annotation: Path to the GENCODE gtf parquet matching Atlas's own
            gene annotation version (used to resolve bare gene_ids to the
            versioned ones Atlas returns, and to drop genes it doesn't know).
        interval: genome.Interval to query densely (mutually exclusive with variants_parquets).
        variants_parquets: List of parquet paths with id_var, chrom, pos, ref, alt columns.
        gene_ids: List of (bare) Ensembl gene IDs to restrict scoring to.
        requested_scorer: Atlas scorer name to request (default "RNA_SEQ").
        max_workers: Max concurrent Atlas requests (respects Snakemake threads).
        batch_size: Number of variants per Atlas query (only used in variants mode).
        chunk_size: Width in bp of each cached, grid-aligned interval chunk
            (only used in interval mode).
        chunk_cache_dir: Directory for cached interval chunks (only used in
            interval mode). Pass the *same* directory across runs targeting
            overlapping regions (e.g. a small dev window and the eventual
            full window) so they share already-downloaded chunks. Defaults
            to a sibling of `output` if not given.
    """
    client = atlas.create(api_key)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    versioned_gene_ids, versioned_by_bare = load_gencode_gene_ids(gene_annotation, gene_ids)
    if not versioned_gene_ids:
        raise ValueError(f"None of the requested genes {gene_ids} were found in {gene_annotation}")

    if variants_parquets is not None:
        variants, id_vars, assigned_gene_ids = load_variants(variants_parquets)

        # A variant assigned to a gene absent from GENCODE can never produce a
        # matching row (see postprocess), so it's also never marked "done" by
        # already_done_id_vars - querying it anyway wastes a request now AND
        # forces it to be re-queried on every future resume. Drop it upfront.
        n_total = len(variants)
        keep = [g in versioned_by_bare for g in assigned_gene_ids]
        variants = [v for v, k in zip(variants, keep) if k]
        id_vars = [i for i, k in zip(id_vars, keep) if k]
        assigned_gene_ids = [g for g, k in zip(assigned_gene_ids, keep) if k]
        print(
            f"Skipping {n_total - len(variants)}/{n_total} variants whose assigned gene "
            "isn't in GENCODE (would never match)",
            flush=True,
        )

        assigned_gene_by_id_var = dict(zip(id_vars, assigned_gene_ids))
        print(f"Querying {len(variants)} variants for {len(versioned_gene_ids)} genes (scorer={requested_scorer})...", flush=True)

        # Sibling of output, not inside it, so Snakemake rerunning the rule
        # doesn't delete already-downloaded batches.
        temp_dir = output.parent / f"{output.stem}.tmp"
        batch_files = query_variants_in_batches(
            client,
            variants,
            id_vars,
            assigned_gene_by_id_var,
            requested_scorer,
            max_workers,
            batch_size,
            temp_dir,
        )

        print(f"Merging {len(batch_files)} batches...", flush=True)
        adata = anndata.concat(
            [anndata.read_h5ad(f) for f in batch_files], join="outer"
        )
        n_unmatched = len(set(id_vars) - set(adata.obs["id_var"]))
        print(
            f"{n_unmatched}/{len(id_vars)} variants had no row for their assigned gene "
            "(gene absent from GENCODE, or no overlap in Atlas's own annotation)",
            flush=True,
        )
        print(f"Temp batches preserved at {temp_dir} for resumption", flush=True)
    else:
        print(
            f"Querying {interval} for {len(versioned_gene_ids)} genes "
            f"(scorer={requested_scorer}, chunk_size={chunk_size})...",
            flush=True,
        )
        # Sibling of output, not inside it, so Snakemake rerunning the rule
        # doesn't delete already-cached chunks.
        if chunk_cache_dir is None:
            chunk_cache_dir = output.parent / f"{output.stem}.chunks"
        else:
            chunk_cache_dir = Path(chunk_cache_dir)

        chunk_files = query_interval_in_chunks(
            client,
            interval,
            requested_scorer,
            max_workers,
            chunk_size,
            chunk_cache_dir,
        )

        # No merge: for a wide interval, concatenating every chunk into one
        # anndata (let alone writing it out) is itself a multi-hundred-GB,
        # very slow operation - and unnecessary, since each chunk is already
        # a queryable, self-contained h5ad. Downstream code should read only
        # the chunks it needs (see extract_alphagenome_atlas_genexpr.py) via
        # grid_chunks() over its own region/gene of interest, rather than
        # loading the whole window at once. `output` is just an empty marker
        # (inside chunk_cache_dir) recording that this interval is fully cached.
        output.touch()
        print(
            f"Touched {output} ({len(chunk_files)} chunks cached in {chunk_cache_dir})",
            flush=True,
        )
        return

    print(f"Final: {adata.n_obs} variant-gene pairs x {adata.n_vars} tracks", flush=True)
    adata.write_h5ad(output, compression="gzip")
    print(f"Wrote {output}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download AlphaGenome Atlas variant effect scores on gene expression (GTEx tracks)."
    )
    query_mode = parser.add_mutually_exclusive_group(required=True)
    query_mode.add_argument(
        "--interval",
        help="Genomic interval to densely query, as 'chrom:start-end' (0-based, half-open)",
    )
    query_mode.add_argument(
        "--variants",
        nargs="+",
        help="One or more parquets with id_var, chrom, pos, ref, alt columns (e.g. seen + unseen) to query only those variants",
    )
    gene_mode = parser.add_mutually_exclusive_group(required=True)
    gene_mode.add_argument(
        "--gene-ids",
        nargs="+",
        help="Bare Ensembl gene IDs to restrict scoring to",
    )
    gene_mode.add_argument(
        "--gene-ids-file",
        help="Text file with one bare Ensembl gene ID per line (for gene lists too long for argv)",
    )
    parser.add_argument(
        "--gene-annotation",
        required=True,
        help="GENCODE gtf parquet matching Atlas's gene annotation version (e.g. gencode.v46)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path: an .h5ad in --variants mode, or an empty marker file (touched once cached) in --interval mode",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="AlphaGenome Atlas API key (defaults to ALPHAGENOME_API_KEY env var)",
    )
    parser.add_argument(
        "--requested-scorer",
        default="RNA_SEQ",
        help="Atlas scorer to request",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="Max concurrent Atlas requests (respects Snakemake threads)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Variants per Atlas query in --variants mode (keeps quota failures cheap to retry)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Width in bp of each cached, grid-aligned interval chunk (only used in --interval mode)",
    )
    parser.add_argument(
        "--chunk-cache-dir",
        default=None,
        help=(
            "Directory for cached interval chunks (only used in --interval mode). "
            "Pass the same directory across runs over overlapping regions to reuse "
            "already-downloaded chunks instead of re-querying them. "
            "Defaults to a sibling of --output."
        ),
    )

    args = parser.parse_args()

    api_key = args.api_key or os.environ["ALPHAGENOME_API_KEY"]

    if args.gene_ids_file:
        gene_ids = Path(args.gene_ids_file).read_text().split()
    else:
        gene_ids = args.gene_ids

    main(
        output=args.output,
        api_key=api_key,
        gene_annotation=args.gene_annotation,
        interval=parse_interval(args.interval) if args.interval else None,
        variants_parquets=args.variants,
        gene_ids=gene_ids,
        requested_scorer=args.requested_scorer,
        max_workers=args.max_workers,
        batch_size=args.batch_size,
        chunk_size=args.chunk_size,
        chunk_cache_dir=args.chunk_cache_dir,
    )
