rule subset_hbb_wholeblood_genexpr_personalized_dev:
    """Subset the dev batch_*.parquet directory down to HBB gene rows in the GTEx whole blood track,
    across all samples, into a single parquet - via duckdb, so none of the batch files
    are ever loaded into pandas at once."""
    input:
        config["alphagenome_genexpr"]["paths"]["dev"],
    output:
        config["alphagenome_genexpr"]["paths"]["hbb_wholeblood_dev"],
    resources:
        runtime = 15,
        mem_mb = 4000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/02-preprocess_data/scripts/subset_genexpr_personalized.py \
            --input-dir {input} \
            --gene-name HBB \
            --biosample-name "venous blood" \
            --output {output}

        echo "Done!"
        """


rule subset_hbb_wholeblood_genexpr_personalized_full:
    """Same as subset_hbb_wholeblood_genexpr_personalized_dev, but over the full-scale (all 3,202
    samples, 1,048,576bp window) genexpr_personalized directory."""
    input:
        config["alphagenome_genexpr"]["paths"]["full"],
    output:
        config["alphagenome_genexpr"]["paths"]["hbb_wholeblood_full"],
    resources:
        runtime = 60,
        mem_mb = 16000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/02-preprocess_data/scripts/subset_genexpr_personalized.py \
            --input-dir {input} \
            --gene-name HBB \
            --biosample-name "venous blood" \
            --output {output}

        echo "Done!"
        """


rule predict_alphagenome_genexpr_dev:
    """Dev-scale: personalized-sequence AlphaGenome forward passes (local model, GPU) for a
    small set of 1000 Genomes individuals over the HBB TSS +/- 8192bp window (16,384bp,
    AlphaGenome's smallest supported input length). Fast to iterate on before scaling to
    the full 1,048,576bp window / full cohort (predict_alphagenome_genexpr_full)."""
    input:
        vcf = config["alphagenome_genexpr"]["vcf"],
        fasta = config["alphagenome_genexpr"]["fasta"],
        gtf = config["alphagenome_genexpr"]["gene_annotation"],
        weights = config["alphagenome_jax"]["paths"]["weights"],
        gene_ids_file = config["alphagenome_genexpr"]["gene_ids"],
    output:
        directory(config["alphagenome_genexpr"]["paths"]["dev"]),
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        tss = config["alphagenome_atlas"]["hbb_window"]["tss"],
        window_size = config["alphagenome_genexpr"]["dev"]["window_size"],
        batch_size = config["alphagenome_genexpr"]["dev"]["batch_size"],
        samples = config["alphagenome_genexpr"]["samples_dev"],
    resources:
        runtime = 60,
        mem_mb = 32000,
        # sbatch (unlike salloc) requires an exact GPU type; a small MIG slice
        # comfortably fits the 16,384bp dev window (verified manually).
        gres = "gpu:1g.24gb:1",
        partition = "gpu_diasfrazer",
        qos = "short",
    conda:
        "wigttsoip"
    retries: 3  # only unfinished samples are re-run on retry/resume - see output dir
    shell:
        """
        python workflows/02-preprocess_data/scripts/alphagenome_genexpr.py \
            --vcf {input.vcf} \
            --fasta {input.fasta} \
            --gtf {input.gtf} \
            --weights-dir {input.weights} \
            --gene-ids-file {input.gene_ids_file} \
            --samples {params.samples} \
            --chromosome {params.chromosome} \
            --tss {params.tss} \
            --window-size {params.window_size} \
            --batch-size {params.batch_size} \
            --output {output}

        echo "Done!"
        """


rule predict_alphagenome_genexpr_full:
    """Full-scale: personalized-sequence AlphaGenome forward passes for all 3,202 1000
    Genomes individuals over the full 1,048,576bp HBB window (matches
    alphagenome_atlas.hbb_window.full exactly, for apples-to-apples comparison with the
    Atlas per-variant scores). NOT in rule all - quota/time-bound, trigger explicitly,
    same as download_alphagenome_atlas_hbb_window_genexpr."""
    input:
        vcf = config["alphagenome_genexpr"]["vcf"],
        fasta = config["alphagenome_genexpr"]["fasta"],
        gtf = config["alphagenome_genexpr"]["gene_annotation"],
        weights = config["alphagenome_jax"]["paths"]["weights"],
        gene_ids_file = config["alphagenome_genexpr"]["gene_ids"],
    output:
        directory(config["alphagenome_genexpr"]["paths"]["full"]),
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        tss = config["alphagenome_atlas"]["hbb_window"]["tss"],
        window_size = config["alphagenome_genexpr"]["full"]["window_size"],
        batch_size = config["alphagenome_genexpr"]["full"]["batch_size"],
    resources:
        # Measured in production (job 28372308, 4 consecutive samples): a
        # steady ~50s/sample on H100 for the full 1,048,576bp forward pass +
        # official gene-mask scoring. 3202 samples * 50s ~= 44.5h (~1.85 days) -
        # runtime below is a ~30% margin over that, not a guess at "how slow
        # AlphaGenome is". No `threads`/many-CPU request: unlike a PyTorch
        # DataLoader or parallel remote-API calls (see e.g.
        # publication_variant_interpretation's enhancer_prediction.smk, which
        # only bumps `threads` for those two cases), this script makes one
        # serial predict_sequence() GPU call at a time, so extra CPU cores
        # would sit idle - the single default thread matches that repo's own
        # convention for pure local-GPU-forward-pass rules.
        runtime = int(60*60),  # 60h ceiling, ~1.35x the measured 44.5h estimate
        mem_mb = 64000,
        # sbatch (unlike salloc) requires an exact GPU type, not just "gpu:1".
        # H100 (80GB) verified to comfortably fit a 1,048,576bp forward pass.
        # This cluster's GresTypes name is "h100" (confirmed via sinfo -p
        # gpu_diasfrazer -o "%P %G" and successful salloc/sbatch allocations) -
        # NOT "nvidia_h100_80gb_hbm3" (that longer name is a different
        # vocabulary used by this cluster's REST/accounting validator, not
        # slurm.conf - publication_variant_interpretation's own GPU rules
        # confirm this: their only "nvidia_..."-style GRES strings are in
        # commented-out, abandoned rule variants).
        gres = "gpu:h100:1",
        partition = "gpu_diasfrazer",
        # No explicit qos: this cluster auto-derives QOS from the requested
        # runtime (confirmed: requesting 2 days here previously auto-selected
        # qos=vlong regardless of an explicit qos="marathon"), and
        # publication_variant_interpretation's GPU rules never set qos either.
    conda:
        "wigttsoip"
    retries: 3
    shell:
        """
        SAMPLES=$(python -c "import pysam; print(' '.join(pysam.VariantFile('{input.vcf}').header.samples))")

        python workflows/02-preprocess_data/scripts/alphagenome_genexpr.py \
            --vcf {input.vcf} \
            --fasta {input.fasta} \
            --gtf {input.gtf} \
            --weights-dir {input.weights} \
            --gene-ids-file {input.gene_ids_file} \
            --samples $SAMPLES \
            --chromosome {params.chromosome} \
            --tss {params.tss} \
            --window-size {params.window_size} \
            --batch-size {params.batch_size} \
            --output {output}

        echo "Done!"
        """


rule subset_allgenes_genexpr_personalized_full:
    """Same as subset_hbb_wholeblood_genexpr_personalized_full, but keeping every gene (no
    --gene-name filter) - the GPU forward pass already scored every gene in gene_ids.txt via
    the official gene-mask scorer, so this is just a re-subset of the same already-computed
    full/ batch directory, no new GPU compute."""
    input:
        config["alphagenome_genexpr"]["paths"]["full"],
    output:
        config["alphagenome_genexpr"]["paths"]["allgenes_full"],
    resources:
        runtime = 60,
        mem_mb = 16000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/02-preprocess_data/scripts/subset_genexpr_personalized.py \
            --input-dir {input} \
            --biosample-name "venous blood" \
            --output {output}

        echo "Done!"
        """
