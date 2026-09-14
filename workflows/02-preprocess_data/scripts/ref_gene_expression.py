"""One-time reference-sequence AlphaGenome forward pass to extract each gene's baseline
(reference-allele) RNA_seq expression level, via the official RNA_SEQ_ACTIVE gene-mask
scorer (variant_scorers.GeneMaskActiveScorer: max(ref_mean, alt_mean) per gene, over exon
mask) scored ref-against-itself - which collapses exactly to the reference-allele masked-mean
level, since max(x, x) = x.

alphagenome_genexpr.py's per-sample forward passes already compute this same reference
prediction internally (ref_track) once per run, but only ever save ref-vs-alt log-fold-change
scores (raw_score) - the absolute reference-expression level itself was never persisted.
This script recomputes that single reference-sequence forward pass (cheap: one forward
pass, not per-sample) and saves it, so downstream analysis can bin genes by baseline
expression independent of any individual's variants.
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

RNA_SEQ_ACTIVE_SCORER = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ_ACTIVE"]


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
    """Same GeneVariantScorer construction alphagenome_genexpr.py uses - exon mask,
    TSS-selected transcripts, INTERVAL_CONTAINED query."""
    gtf = pd.read_parquet(gtf_path)
    extractor = gene_mask_extractor_lib.GeneMaskExtractor(
        gtf=gtf,
        gene_mask_type=gene_mask_extractor_lib.GeneMaskType.EXONS,
        gene_query_type=gene_mask_extractor_lib.GeneQueryType.INTERVAL_CONTAINED,
    )
    return gene_mask_lib.GeneVariantScorer(gene_mask_extractor=extractor)


def main(fasta_path, gtf_path, weights_dir, gene_ids, chromosome, tss, window_size, output):
    model = load_model(weights_dir, fasta_path)
    scorer = build_gene_scorer(gtf_path)

    half = window_size // 2
    start = tss - half
    end = start + window_size
    interval = Interval(chromosome, start, end)

    # Only used to satisfy GeneMaskExtractor's API (IndelMask.from_variant); the EXONS
    # gene mask itself never depends on the variant, and for a non-indel placeholder
    # (ref/alt both length 1) align_alternate() is a no-op regardless of position - same
    # reasoning as alphagenome_genexpr.py's own placeholder_variant.
    placeholder_variant = Variant(chromosome, start + 1, "A", "T")
    masks, mask_metadata = scorer.get_masks_and_metadata(
        interval, placeholder_variant, settings=RNA_SEQ_ACTIVE_SCORER, track_metadata=None
    )
    print(f"Official gene mask covers {len(mask_metadata)} genes (TSS-selected transcripts)", flush=True)

    with pysam.FastaFile(fasta_path) as fasta:
        ref_seq = fasta.fetch(chromosome, start, end).upper()

    print(f"Predicting reference sequence ({chromosome}:{start}-{end})...", flush=True)
    ref_track = model.predict_sequence(
        sequence=ref_seq,
        requested_outputs=[OutputType.RNA_SEQ],
        ontology_terms=None,
        interval=interval,
    ).rna_seq

    # Score ref against itself: GeneMaskActiveScorer returns max(ref_mean, alt_mean) per
    # gene/track, which with ref==alt is exactly the reference-allele masked-mean level.
    ref_values = np.asarray(ref_track.values, dtype=np.float32)
    ref_input = {OutputType.RNA_SEQ: ref_values}
    scores = scorer.score_variant(ref_input, ref_input, masks=masks, settings=RNA_SEQ_ACTIVE_SCORER)
    result = scorer.finalize_variant(
        scores,
        track_metadata=model.output_metadata(Organism.HOMO_SAPIENS),
        mask_metadata=mask_metadata,
        settings=RNA_SEQ_ACTIVE_SCORER,
    )
    result.uns["interval"] = interval
    result.uns["variant"] = None
    result.uns["variant_scorer"] = RNA_SEQ_ACTIVE_SCORER

    tidy = variant_scorers.tidy_anndata(result, match_gene_strand=True, include_extended_metadata=True)
    if gene_ids is not None:
        tidy = tidy[tidy["gene_id"].isin(gene_ids)]
    tidy = tidy.rename(columns={"raw_score": "ref_expression"})
    # scored_interval is a raw genome.Interval object - not parquet-serializable.
    tidy["scored_interval"] = tidy["scored_interval"].astype(str)

    tidy.to_parquet(output, compression="zstd", index=False)
    print(f"Wrote {output} ({len(tidy)} rows, {tidy['gene_id'].nunique()} genes)", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--weights-dir", required=True)
    parser.add_argument("--gene-ids-file", default=None, help="Optional: restrict output to these bare gene_ids")
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--tss", type=int, required=True)
    parser.add_argument("--window-size", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    gene_ids = None
    if args.gene_ids_file:
        with open(args.gene_ids_file) as f:
            gene_ids = [line.strip() for line in f if line.strip()]

    main(
        fasta_path=args.fasta,
        gtf_path=args.gtf,
        weights_dir=args.weights_dir,
        gene_ids=gene_ids,
        chromosome=args.chromosome,
        tss=args.tss,
        window_size=args.window_size,
        output=args.output,
    )
