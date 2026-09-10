"""Restrict AlphaGenome Atlas single-variant HBB/K562 scores to just the unique SNVs listed
in --variants (see get_unique_variants.py) - Atlas's saturation set scores every possible
alt allele at every position in the window, but only a fraction of those are ever
polymorphic in the sample set actually being joined against, so restricting the singles
table to this subset before the per-sample join is both correct (an unmatched Atlas row
never survives that join anyway) and much smaller/faster to join against.
"""

import argparse
from pathlib import Path

import duckdb


def main(variants_path, singles_glob, output):
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(
        """
        COPY (
            SELECT s.*
            FROM read_parquet(?) s
            WHERE s.variant IN (SELECT variant FROM read_parquet(?))
        ) TO '{output}' (FORMAT parquet, COMPRESSION zstd)
        """.format(output=output),
        [singles_glob, variants_path],
    )
    n_rows = con.execute("SELECT count(*) FROM read_parquet(?)", [output]).fetchone()[0]
    print(f"Wrote {output} ({n_rows} rows)", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", required=True, help="get_unique_variants.py output parquet")
    parser.add_argument("--singles", required=True, help="Atlas long directory of batch_*.parquet (or a glob/single file)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    singles_glob = args.singles
    if Path(singles_glob).is_dir():
        singles_glob = str(Path(singles_glob) / "batch_*.parquet")

    main(args.variants, singles_glob, args.output)
