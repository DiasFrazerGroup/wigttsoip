"""Join, per sample, gene, and track, the model's predicted effect of an individual's whole
combination of HBB-window SNVs (genexpr_personalized, "whole_score") against the sum of
those same SNVs' own single-variant Atlas effects (alphagenome_atlas, "sum_parts_score") -
the additive/no-epistasis null. Keyed by gene_id too (not just sample/track): a variant can
score differently against different nearby genes, so when singles/combined cover more than
one gene (see the allgenes_* config paths), the join must match each combined row's own
gene_id, not just its track - for the HBB-only paths this is a no-op (single gene_id value).

A track is identified by track_name; per alphagenome.models.variant_scorers.tidy_anndata's
own docstring, (track_name, track_strand) is the minimal key that uniquely identifies a
track in general - e.g. some biosamples have two distinct "polyA plus RNA-seq" tracks (one
strand='.', one strand='-'). The join below matches on track_name ALONE (not also
track_strand) as a defensive default: for the GTEx whole blood polyA track used here, both
Atlas's singles table and alphagenome_genexpr.py's personalized-sequence model agree on
strand='.' (confirmed directly), so this doesn't currently collapse anything that should
have stayed distinct - but track_name alone is kept as the join key rather than reintroducing
a track_strand equality check, since Atlas's singles table never has more than one
track_strand value per track_name (so this can't create spurious duplicate matches via
fan-out either way). Output keeps combined's own track_strand value (not singles'), since
that's the one meaningful for match_gene_strand filtering. No averaging across duplicate
(sample, track) rows otherwise: each output row is one real track's whole_score/sum_parts_score pair.
n_variants vs n_variants_matched flags combinations with a variant outside the singles
window/track set, where sum_parts_score would otherwise silently understate the true sum.
"""

import argparse
from pathlib import Path

import duckdb

SQL = """
COPY (
    -- combined_rows carries the ORIGINAL variant_id string (a semicolon-joined list of
    -- every SNV a sample carries in the window - up to ~3,800 for the full 1,048,576bp
    -- window, tens of KB of text) exactly once per (sample, track) row. row_id is a
    -- lightweight surrogate key so the explode/join stages below never need to touch
    -- variant_id at all - carrying that giant string through the explode alongside the
    -- UNNEST is what actually OOM'd the full-scale run (DuckDB duplicates every other
    -- selected column across each row UNNEST produces, so that string got copied into
    -- every one of a sample's thousands of exploded rows: tens of GB of pure duplication,
    -- dwarfing the join itself). variant_id is reattached once at the very end instead.
    WITH combined_rows AS (
        SELECT
            ROW_NUMBER() OVER () AS row_id,
            sample,
            gene_id,
            variant_id,
            track_name,
            track_strand,
            raw_score AS whole_score,
            CASE WHEN variant_id = '' THEN 0 ELSE len(string_split(variant_id, ';')) END AS n_variants
        FROM read_parquet(?)
    ),
    exploded AS (
        SELECT
            row_id,
            gene_id,
            track_name,
            UNNEST(
                CASE WHEN variant_id = '' THEN []::VARCHAR[] ELSE string_split(variant_id, ';') END
            ) AS variant
        FROM combined_rows
    ),
    -- alphagenome_genexpr.py's variant_id labels are 0-based (pysam's rec.start,
    -- matching the 0-based sequence-splicing math it does internally), but
    -- AlphaGenome Atlas's own `variant` field is 1-based (standard VCF POS) -
    -- verified directly against the reference fasta. Converting here (rather than in
    -- alphagenome_genexpr.py) keeps every already-written batch_*.parquet from the
    -- expensive full GPU run consistent with every batch still to come.
    exploded_1based AS (
        SELECT
            row_id,
            gene_id,
            track_name,
            concat(
                split_part(variant, ':', 1), ':',
                (split_part(variant, ':', 2)::BIGINT + 1)::VARCHAR, ':',
                split_part(variant, ':', 3)
            ) AS variant
        FROM exploded
    ),
    matched AS (
        SELECT
            e.row_id,
            s.score AS single_score
        FROM exploded_1based e
        LEFT JOIN read_parquet(?) s
            ON e.variant = s.variant
            AND e.track_name = s.track_name
            -- combined's gene_id is unversioned ("ENSG...") but Atlas's singles gene_id is
            -- versioned ("ENSG....N", straight from the h5ad obs/gene_id column) - strip
            -- the version before comparing, confirmed directly (an unstripped join here
            -- silently matched zero rows: LEFT JOIN never errors, it just leaves every row
            -- with a NULL score).
            AND e.gene_id = split_part(s.gene_id, '.', 1)
    ),
    aggregated AS (
        SELECT
            row_id,
            SUM(single_score) AS sum_parts_score,
            COUNT(single_score) AS n_variants_matched
        FROM matched
        GROUP BY row_id
    )
    SELECT
        c.sample,
        c.gene_id,
        c.variant_id,
        c.track_name,
        c.track_strand,
        c.whole_score,
        a.sum_parts_score,
        c.n_variants,
        a.n_variants_matched
    FROM combined_rows c
    JOIN aggregated a USING (row_id)
) TO '{output}' (FORMAT parquet, COMPRESSION zstd)
"""


def main(combined_path, singles_path, output, memory_limit_mb=None):
    # singles_path is a directory of zstd-compressed batch_*.parquet files (see
    # extract_to_batches() in extract_alphagenome_atlas_genexpr.py) - duckdb's
    # read_parquet accepts a glob string directly, so expand the directory to one here.
    if Path(singles_path).is_dir():
        singles_path = str(Path(singles_path) / "batch_*.parquet")

    con = duckdb.connect()
    if memory_limit_mb is not None:
        # Match the SLURM job's own --mem request (see the rule's resources.mem_mb) so
        # duckdb's buffer manager caps itself and spills/errors gracefully instead of the
        # whole exploded-variant join getting OS OOM-killed (as happened on the full-scale
        # run before this was wired in - the UNNEST explosion over the full 1Mb window's
        # per-sample variant lists is what actually needs the memory, not the final
        # row counts, which stay modest).
        con.execute(f"SET memory_limit = '{int(memory_limit_mb)}MB'")

    # duckdb's COPY ... TO doesn't support parameterizing the destination path via `?`
    # (it silently mis-binds the other placeholders when tried) - the output path is a
    # trusted Snakemake-provided path, not user input, so a plain format() is fine here;
    # the two source paths are still passed as real bind parameters.
    con.execute(SQL.format(output=output), [combined_path, singles_path])
    n_rows = con.execute("SELECT count(*) FROM read_parquet(?)", [output]).fetchone()[0]
    print(f"Wrote {output} ({n_rows} rows)", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--combined", required=True, help="HBB/whole blood genexpr_personalized parquet (02-preprocess_data)")
    parser.add_argument("--singles", required=True, help="HBB/whole blood Atlas variant-effects long parquet (01-obtain_data)")
    parser.add_argument("--output", required=True)
    parser.add_argument("--memory-limit-mb", type=int, default=None, help="duckdb memory_limit, in MB - match the job's own --mem/resources.mem_mb")
    args = parser.parse_args()

    main(
        combined_path=args.combined,
        singles_path=args.singles,
        output=args.output,
        memory_limit_mb=args.memory_limit_mb,
    )
