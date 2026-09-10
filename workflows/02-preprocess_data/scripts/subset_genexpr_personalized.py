"""Subset a directory of batch_*.parquet files (from alphagenome_genexpr.py) into a
single parquet, filtered by gene name and/or biosample name, using duckdb so none of
the batch files (nor their union) ever needs to be loaded into pandas at once."""

import argparse

import duckdb


def main(input_dir, output, gene_name=None, biosample_name=None):
    con = duckdb.connect()
    glob_pattern = f"{input_dir}/batch_*.parquet"

    where = []
    params = []
    if gene_name is not None:
        where.append("gene_name = ?")
        params.append(gene_name)
    if biosample_name is not None:
        where.append("biosample_name = ?")
        params.append(biosample_name)
    where_clause = f"WHERE {' AND '.join(where)}" if where else ""

    query = f"""
        COPY (
            SELECT * FROM read_parquet(?)
            {where_clause}
        ) TO '{output}' (FORMAT parquet, COMPRESSION zstd)
    """
    con.execute(query, [glob_pattern] + params)

    n_rows = con.execute(f"SELECT count(*) FROM read_parquet('{output}')").fetchone()[0]
    print(f"Wrote {output} ({n_rows} rows)", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, help="Directory of batch_*.parquet files")
    parser.add_argument("--output", required=True)
    parser.add_argument("--gene-name", default=None)
    parser.add_argument("--biosample-name", default=None)
    args = parser.parse_args()

    main(
        input_dir=args.input_dir,
        output=args.output,
        gene_name=args.gene_name,
        biosample_name=args.biosample_name,
    )
