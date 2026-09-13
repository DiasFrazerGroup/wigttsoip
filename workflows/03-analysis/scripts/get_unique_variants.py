"""Distinct SNVs (chrom:pos:ref>alt, 1-based - matching AlphaGenome Atlas's own `variant`
convention directly) actually carried by a set of 1000 Genomes samples in a genomic window.

Deliberately independent of any downstream track/score data - just the variant identities -
so it can feed the Atlas singles-subsetting step without depending on anything else.
"""

import argparse
from pathlib import Path

import duckdb
import pysam


def unique_snvs(vcf_path, chromosome, start, end, samples=None):
    """Distinct (pos, ref, alt) for every SNV (ref/alt both length 1) actually carried by
    >=1 sample in `samples` (unphased: any ALT copy counts).

    If samples is None, every SNV row in the VCF counts as-is: a jointly-called multi-sample
    VCF only ever lists a site because >=1 of ALL its samples carries the ALT there, so "the
    full cohort" needs no per-sample GT check at all - only a genuine subset of samples does.
    """
    out = set()
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf.fetch(chromosome, start, end):
            ref = rec.ref
            if not ref or len(ref) != 1 or ref.upper() not in "ACGT":
                continue
            for ai, alt in enumerate(rec.alts or (), start=1):
                if not alt or len(alt) != 1 or alt.upper() not in "ACGT":
                    continue
                if samples is not None:
                    carried = any(
                        (gt := rec.samples[s].get("GT")) and ai in gt for s in samples
                    )
                    if not carried:
                        continue
                out.add((rec.pos, ref.upper(), alt.upper()))
    return out


def main(vcf_path, chromosome, start, end, output, samples=None):
    snvs = unique_snvs(vcf_path, chromosome, start, end, samples=samples)
    scope = f"{len(samples)} samples" if samples else "full cohort"
    print(f"{len(snvs)} unique SNVs carried in {chromosome}:{start}-{end} ({scope})", flush=True)
    variants = sorted(f"{chromosome}:{pos}:{ref}>{alt}" for pos, ref, alt in snvs)

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute("CREATE TABLE unique_variants (variant VARCHAR)")
    con.executemany("INSERT INTO unique_variants VALUES (?)", [(v,) for v in variants])
    con.execute("COPY unique_variants TO '{}' (FORMAT parquet, COMPRESSION zstd)".format(output))
    print(f"Wrote {output} ({len(variants)} unique variants)", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="1000 Genomes VCF")
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--samples", nargs="*", default=None, help="Restrict to SNVs carried by these samples; omit for the full VCF cohort")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    main(args.vcf, args.chromosome, args.start, args.end, args.output, samples=args.samples)
