"""Annotate a set of unique HBB-window SNVs with gnomAD v3.1.1 joint allele frequency and
minor allele frequency (MAF = min(AF, 1-AF)).

gnomAD itself is only available as a (huge, but tabix-indexed) VCF, so pysam does the one
region fetch that's unavoidable - but it's a single indexed fetch covering the whole window
(~34s for the full 1,048,576bp HBB window, confirmed directly), followed by in-memory dict
lookups, not one query per variant. Everything that touches parquet (reading the distinct
variant list in, writing the annotated table out) goes through duckdb rather than pandas.
"""

import argparse
from pathlib import Path

import duckdb
import pysam


def gnomad_af_lookup(gnomad_vcf, chromosome, start, end):
    """{"chrom:pos:ref>alt": (af, ac, an)} for every SNV in the window - one indexed fetch,
    not one per variant."""
    lookup = {}
    with pysam.VariantFile(gnomad_vcf) as vf:
        for rec in vf.fetch(chromosome, start, end):
            afs = rec.info.get("AF")
            if afs is None:
                continue
            acs = rec.info.get("AC")
            an = rec.info.get("AN")
            for ai, alt in enumerate(rec.alts or ()):
                if not alt or len(alt) != 1 or len(rec.ref) != 1:
                    continue
                key = f"{chromosome}:{rec.pos}:{rec.ref.upper()}>{alt.upper()}"
                lookup[key] = (afs[ai], acs[ai] if acs is not None else None, an)
    return lookup


def main(variants_path, gnomad_vcf, chromosome, start, end, output):
    con = duckdb.connect()
    variants = [
        row[0]
        for row in con.execute(
            "SELECT DISTINCT variant FROM read_parquet(?)", [variants_path]
        ).fetchall()
    ]
    print(f"Annotating {len(variants)} unique variants with gnomAD AF...", flush=True)

    lookup = gnomad_af_lookup(gnomad_vcf, chromosome, start, end)

    rows = []
    for v in variants:
        af, ac, an = lookup.get(v, (None, None, None))
        maf = min(af, 1 - af) if af is not None else None
        rows.append((v, af, ac, an, maf))
    n_found = sum(1 for r in rows if r[1] is not None)
    print(f"{n_found}/{len(rows)} variants found in gnomAD ({len(rows) - n_found} missing)", flush=True)

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    con.execute(
        "CREATE TABLE annotated (variant VARCHAR, gnomad_af DOUBLE, gnomad_ac BIGINT, gnomad_an BIGINT, gnomad_maf DOUBLE)"
    )
    con.executemany("INSERT INTO annotated VALUES (?, ?, ?, ?, ?)", rows)
    con.execute("COPY annotated TO '{}' (FORMAT parquet, COMPRESSION zstd)".format(output))
    print(f"Wrote {output}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", required=True, help="subset_unique_variant_scores.py output parquet")
    parser.add_argument("--gnomad-vcf", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    main(args.variants, args.gnomad_vcf, args.chromosome, args.start, args.end, args.output)
