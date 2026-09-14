"""Personalized-sequence AlphaGenome forward passes for real 1000 Genomes individuals.

For each sample, apply every SNV they carry (unphased - no maternal/paternal split,
indels/SVs skipped) onto the reference sequence, run both the reference and the
personalized sequence through the local AlphaGenome model (GPU), and score the
predicted RNA_SEQ change per gene per track using AlphaGenome's own official
gene-mask scorer (alphagenome_research.model.variant_scoring.gene_mask.GeneVariantScorer
+ variant_scorers.GeneMaskLFCScorer, the same machinery behind
RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']) - not a hand-rolled exon-sum - so scores are
directly comparable to the official Atlas/paper conventions: gene exon mask from
TSS-selected transcripts, gene/track strand-consistency filtering, and a masked-mean
log fold change (natural log, 1e-3 pseudocount).

No target/observed expression data involved - pure forward-pass inference.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from alphagenome.data.genome import Interval, Variant
from alphagenome.models import variant_scorers
from alphagenome.models.dna_model import Organism
from alphagenome.models.dna_output import OutputType
from alphagenome_research.model.dna_model import OrganismSettings
from alphagenome_research.model.dna_model import create as create_model
from alphagenome_research.model.variant_scoring import gene_mask as gene_mask_lib
from alphagenome_research.model.variant_scoring import (
    gene_mask_extractor as gene_mask_extractor_lib,
)

RNA_SEQ_SCORER = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]  # GeneMaskLFCScorer


def load_model(weights_dir, fasta_path):
    import jax

    weights_dir = str(Path(weights_dir).resolve())  # orbax requires an absolute path
    print(f"Checkpoint path: {weights_dir}", flush=True)
    print(f"Available JAX devices: {jax.devices()}", flush=True)
    organism_settings = {
        Organism.HOMO_SAPIENS: OrganismSettings(fasta_path=str(fasta_path))
    }
    model = create_model(weights_dir, organism_settings=organism_settings)
    print("Model loaded successfully", flush=True)
    return model


def build_gene_scorer(gtf_path):
    """The same GeneVariantScorer construction AlphaGenomeModel.create() uses
    internally for GENE_MASK_LFC/GENE_MASK_ACTIVE - exon mask, TSS-selected
    transcripts, INTERVAL_CONTAINED query."""
    gtf = pd.read_parquet(gtf_path)
    extractor = gene_mask_extractor_lib.GeneMaskExtractor(
        gtf=gtf,
        gene_mask_type=gene_mask_extractor_lib.GeneMaskType.EXONS,
        gene_query_type=gene_mask_extractor_lib.GeneQueryType.INTERVAL_CONTAINED,
    )
    return gene_mask_lib.GeneVariantScorer(gene_mask_extractor=extractor)


def load_sample_snvs(vcf_path, chrom, start, end, sample):
    """(pos, ref, alt) for every SNV where `sample` carries >=1 ALT copy.

    Genotypes here are unphased, so we don't distinguish maternal/paternal -
    just record every alternate allele actually carried. Indels/SVs are skipped
    (REF/ALT both length 1 only), so positions never shift.
    """
    out = []
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf.fetch(chrom, start, end):
            gt = rec.samples[sample].get("GT")
            if not gt or all(a in (0, None) for a in gt):
                continue
            for ai in set(a for a in gt if a):
                ref, alt = rec.ref, rec.alleles[ai]
                if ref and alt and len(ref) == 1 and len(alt) == 1 and alt.upper() in "ACGT":
                    out.append((rec.start, ref.upper(), alt.upper()))
    return out


def build_personal_sequence(ref_seq, start, snvs):
    """Apply every SNV in place. Position never shifts: SNV-only, ref/alt both len 1."""
    seq = list(ref_seq)
    for pos, ref, alt in snvs:
        local = pos - start
        if 0 <= local < len(seq) and seq[local] == ref:
            seq[local] = alt
    return "".join(seq)


def predict_rna(model, sequence, interval):
    return model.predict_sequence(
        sequence=sequence,
        requested_outputs=[OutputType.RNA_SEQ],
        ontology_terms=None,
        interval=interval,
    ).rna_seq


def score_ref_vs_alt(model, scorer, masks, mask_metadata, interval, ref_values, alt_values, variant_str):
    """Official GeneMaskLFCScorer scoring of ref vs. alt RNA_SEQ tracks, tidied to a
    long-format dataframe via variant_scorers.tidy_anndata (same utility the official
    pipeline uses to format AnnData scores as tidy rows)."""
    ref = {OutputType.RNA_SEQ: np.asarray(ref_values, dtype=np.float32)}
    alt = {OutputType.RNA_SEQ: np.asarray(alt_values, dtype=np.float32)}
    scores = scorer.score_variant(ref, alt, masks=masks, settings=RNA_SEQ_SCORER)
    result = scorer.finalize_variant(
        scores,
        track_metadata=model.output_metadata(Organism.HOMO_SAPIENS),
        mask_metadata=mask_metadata,
        settings=RNA_SEQ_SCORER,
    )
    result.uns["interval"] = interval
    result.uns["variant"] = variant_str
    result.uns["variant_scorer"] = RNA_SEQ_SCORER
    return variant_scorers.tidy_anndata(result, match_gene_strand=True, include_extended_metadata=True)


def already_processed_samples(output_dir):
    done = set()
    for f in Path(output_dir).glob("batch_*.parquet"):
        try:
            done.update(pd.read_parquet(f, columns=["sample"])["sample"].unique())
        except Exception as e:
            print(f"Warning: deleting corrupted file {f}: {e}", flush=True)
            f.unlink(missing_ok=True)
    return done


def next_batch_index(output_dir):
    """First unused batch_*.parquet index in output_dir - a resumed run must continue
    numbering from here, not restart at 0, or it silently overwrites (destroys) earlier
    completed batches with a different run's samples instead of adding to them."""
    existing = sorted(Path(output_dir).glob("batch_*.parquet"))
    if not existing:
        return 0
    return int(existing[-1].stem.split("_")[1]) + 1


def main(
    vcf_path,
    fasta_path,
    gtf_path,
    weights_dir,
    gene_ids,
    samples,
    chromosome,
    tss,
    window_size,
    output_dir,
    batch_size,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model = load_model(weights_dir, fasta_path)
    scorer = build_gene_scorer(gtf_path)

    half = window_size // 2
    start = tss - half
    end = start + window_size
    interval = Interval(chromosome, start, end)

    # Only used to satisfy GeneMaskExtractor's API (IndelMask.from_variant); the
    # EXONS gene mask itself never depends on the variant, and for a non-indel
    # placeholder (ref/alt both length 1) align_alternate() is a no-op regardless
    # of position - so this dummy variant never affects any score.
    placeholder_variant = Variant(chromosome, start + 1, "A", "T")
    masks, mask_metadata = scorer.get_masks_and_metadata(
        interval, placeholder_variant, settings=RNA_SEQ_SCORER, track_metadata=None
    )
    print(f"Official gene mask covers {len(mask_metadata)} genes (TSS-selected transcripts)", flush=True)

    with pysam.FastaFile(fasta_path) as fasta:
        ref_seq = fasta.fetch(chromosome, start, end).upper()

    print(f"Predicting reference sequence ({chromosome}:{start}-{end})...", flush=True)
    ref_track = predict_rna(model, ref_seq, interval)

    done = already_processed_samples(output_dir)
    samples_todo = [s for s in samples if s not in done]
    print(f"{len(done)} samples already processed, {len(samples_todo)} to go", flush=True)

    start_batch_idx = next_batch_index(output_dir)
    n_batches = (len(samples_todo) + batch_size - 1) // batch_size
    for i in range(n_batches):
        batch_idx = start_batch_idx + i
        batch_samples = samples_todo[i * batch_size : (i + 1) * batch_size]
        parts = []
        for sample in batch_samples:
            snvs = load_sample_snvs(vcf_path, chromosome, start, end, sample)
            personal_seq = build_personal_sequence(ref_seq, start, snvs)
            # NOTE: pos is pysam's 0-based rec.start, intentionally NOT Atlas's 1-based
            # convention - summarize_whole_vs_sum_parts.py converts at join time instead,
            # so already-written batch_*.parquet stay consistent with future ones.
            variant_str = ";".join(f"{chromosome}:{pos}:{ref}>{alt}" for pos, ref, alt in snvs)

            print(f"[{sample}] {len(snvs)} SNVs in window, scoring...", flush=True)
            alt_track = predict_rna(model, personal_seq, interval)

            tidy = score_ref_vs_alt(
                model, scorer, masks, mask_metadata, interval,
                ref_track.values, alt_track.values, variant_str,
            )
            if gene_ids is not None:
                tidy = tidy[tidy["gene_id"].isin(gene_ids)]
            tidy.insert(0, "sample", sample)
            parts.append(tidy)

        if parts:
            batch_df = pd.concat(parts, ignore_index=True)
            # tidy_anndata's `scored_interval` is a raw genome.Interval object -
            # not parquet-serializable, so stringify it.
            batch_df["scored_interval"] = batch_df["scored_interval"].astype(str)
            tmp = output_dir / f"batch_{batch_idx:06d}.parquet.tmp"
            final = output_dir / f"batch_{batch_idx:06d}.parquet"
            batch_df.to_parquet(tmp, compression="zstd", index=False)
            tmp.rename(final)
            print(f"Wrote {final} ({len(batch_df)} rows)", flush=True)

    print("Done!", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--weights-dir", required=True)
    parser.add_argument("--gene-ids-file", default=None, help="Optional: restrict output to these bare gene_ids")
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--tss", type=int, required=True)
    parser.add_argument("--window-size", type=int, required=True)
    parser.add_argument("--output", required=True, help="Output directory of batch_*.parquet")
    parser.add_argument("--batch-size", type=int, default=2)
    args = parser.parse_args()

    gene_ids = None
    if args.gene_ids_file:
        with open(args.gene_ids_file) as f:
            gene_ids = [line.strip() for line in f if line.strip()]

    main(
        vcf_path=args.vcf,
        fasta_path=args.fasta,
        gtf_path=args.gtf,
        weights_dir=args.weights_dir,
        gene_ids=gene_ids,
        samples=args.samples,
        chromosome=args.chromosome,
        tss=args.tss,
        window_size=args.window_size,
        output_dir=args.output,
        batch_size=args.batch_size,
    )
