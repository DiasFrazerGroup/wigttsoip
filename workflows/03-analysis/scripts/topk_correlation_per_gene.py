"""For every gene in the window, correlate each individual's combinatorial (whole_score)
effect on that gene against the sum of their top-k strongest (by |score|) carried variants'
own single-variant effects on that same gene - the same "effective dimensionality" question
notebooks/hbb.ipynb answers for HBB alone, generalized across genes.

Everything up through the per-(sample, gene)-ranked cumulative sum runs in duckdb (mirrors
summarize_whole_vs_sum_parts.py's explode/1-based-conversion/join, keyed by gene_id too, but
keeping the whole ranked list per row instead of collapsing straight to one sum) - only the
small (row_id, k, cumsum) result table is pulled into pandas, so this never holds the full
per-variant explosion (potentially tens of millions of rows across 94 genes) in memory at
once.
"""

import argparse
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy import stats

K_VALUES = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 1200, 1500, 2000, 2200, 2500, 3000, 3200]

SQL = """
WITH combined_rows AS (
    -- reads summarize_whole_vs_sum_parts.py's own output directly (already one row per
    -- sample/gene/track, already filtered to a single track), rather than re-deriving
    -- whole_score/variant_id from the raw genexpr_personalized parquet.
    SELECT
        ROW_NUMBER() OVER () AS row_id,
        sample,
        gene_id,
        variant_id,
        whole_score
    FROM read_parquet(?)
),
exploded AS (
    SELECT
        row_id,
        gene_id,
        UNNEST(
            CASE WHEN variant_id = '' THEN []::VARCHAR[] ELSE string_split(variant_id, ';') END
        ) AS variant
    FROM combined_rows
),
-- see summarize_whole_vs_sum_parts.py: alphagenome_genexpr.py's variant_id positions are
-- 0-based, Atlas's own `variant` field is 1-based.
exploded_1based AS (
    SELECT
        row_id,
        gene_id,
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
    JOIN read_parquet(?) s
        ON e.variant = s.variant
        -- combined's gene_id is unversioned but Atlas's singles gene_id is versioned
        -- ("ENSG....N", straight from the h5ad obs/gene_id column) - strip the version
        -- before comparing (see summarize_whole_vs_sum_parts.py for how this was found:
        -- an unstripped join here matches zero rows, silently, since this is an INNER
        -- JOIN it just returns an empty result rather than erroring).
        AND e.gene_id = split_part(s.gene_id, '.', 1)
        AND s.track_name = ?
        AND s.track_strand = ?
),
ranked AS (
    SELECT
        row_id,
        single_score,
        ROW_NUMBER() OVER (PARTITION BY row_id ORDER BY ABS(single_score) DESC) AS rn,
        SUM(single_score) OVER (
            PARTITION BY row_id ORDER BY ABS(single_score) DESC
            ROWS UNBOUNDED PRECEDING
        ) AS cumsum
    FROM matched
),
totals AS (
    SELECT row_id, SUM(single_score) AS total_sum, COUNT(*) AS n_matched
    FROM matched
    GROUP BY row_id
)
SELECT
    c.row_id,
    c.sample,
    c.gene_id,
    c.whole_score,
    t.total_sum,
    t.n_matched,
    r.rn,
    r.cumsum
FROM combined_rows c
JOIN totals t USING (row_id)
LEFT JOIN ranked r USING (row_id)
WHERE r.rn IS NULL OR r.rn IN ({k_list})
"""


def _correlations_for_group(group, k_values):
    """One row per k for a single gene: Pearson/Spearman of (top-k sum) vs. whole_score."""
    out = []
    y = group["whole_score"].to_numpy()
    for k in list(k_values) + ["all"]:
        if k == "all":
            x = group["total_sum"].to_numpy()
        else:
            # last cumsum at or before rank k, per sample; falls back to total_sum for
            # samples with fewer than k matched variants (their full sum already is the
            # "top-k" sum once k exceeds how many variants they have).
            x = group[f"cumsum_{k}"].to_numpy()
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 3:
            r = p = rho = p_rho = np.nan
        else:
            r, p = stats.pearsonr(x[mask], y[mask])
            rho, p_rho = stats.spearmanr(x[mask], y[mask])
        out.append({
            # stringified so the column stays a single type once "all" is mixed in
            # alongside the integer k values (parquet/arrow can't write a mixed column).
            "k": str(k), "pearson_r": r, "pearson_p": p,
            "spearman_rho": rho, "spearman_p": p_rho, "n": int(mask.sum()),
        })
    return pd.DataFrame(out)


def main(whole_vs_sum_path, singles_path, output, track_name, track_strand, memory_limit_mb=None):
    if Path(singles_path).is_dir():
        singles_path = str(Path(singles_path) / "batch_*.parquet")

    con = duckdb.connect()
    if memory_limit_mb is not None:
        con.execute(f"SET memory_limit = '{int(memory_limit_mb)}MB'")

    k_list = ", ".join(str(k) for k in K_VALUES)
    long_df = con.execute(
        SQL.format(k_list=k_list),
        [whole_vs_sum_path, singles_path, track_name, track_strand],
    ).df()
    print(f"Pulled {len(long_df)} (row_id, k) rows into pandas", flush=True)

    # one row per (row_id=sample,gene), columns cumsum_<k> - the running sum at the last
    # rank <= k (ffill handles samples whose n_matched < k, where no exact rn==k row
    # exists: their cumsum already reached their own total by their own last rank).
    base = long_df[["row_id", "sample", "gene_id", "whole_score", "total_sum", "n_matched"]].drop_duplicates("row_id").set_index("row_id")
    pivot = long_df.dropna(subset=["rn"]).pivot_table(index="row_id", columns="rn", values="cumsum")
    pivot = pivot.reindex(columns=K_VALUES).ffill(axis=1)
    pivot.columns = [f"cumsum_{k}" for k in K_VALUES]
    wide = base.join(pivot)
    for k in K_VALUES:
        col = f"cumsum_{k}"
        if col not in wide:
            wide[col] = np.nan
        wide[col] = wide[col].fillna(wide["total_sum"])

    results = (
        wide.groupby("gene_id", group_keys=True)
        .apply(lambda g: _correlations_for_group(g, K_VALUES), include_groups=False)
        .reset_index(level=0)
        .reset_index(drop=True)
    )
    results.to_parquet(output, compression="zstd", index=False)
    print(f"Wrote {output} ({len(results)} rows, {results['gene_id'].nunique()} genes)", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--whole-vs-sum", required=True, help="whole_vs_sum_parts_allgenes parquet (summarize_whole_vs_sum_parts.py's own output)")
    parser.add_argument("--singles", required=True, help="allgenes Atlas variant-effects long parquet/directory (01-obtain_data)")
    parser.add_argument("--output", required=True)
    parser.add_argument("--track-name", required=True)
    parser.add_argument("--track-strand", required=True)
    parser.add_argument("--memory-limit-mb", type=int, default=None)
    args = parser.parse_args()

    main(
        whole_vs_sum_path=args.whole_vs_sum,
        singles_path=args.singles,
        output=args.output,
        track_name=args.track_name,
        track_strand=args.track_strand,
        memory_limit_mb=args.memory_limit_mb,
    )
