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
        # ~50s/sample measured on H100 x 3202 samples ~= 44.5h; runtime below is a ~30%
        # margin over that. Single-threaded: one serial predict_sequence() GPU call, so
        # extra CPU cores would just idle.
        runtime = int(60*60),
        mem_mb = 64000,
        # This cluster's GresTypes name is "h100" (per slurm.conf), not the longer
        # "nvidia_h100_80gb_hbm3" some accounting tools use.
        gres = "gpu:h100:1",
        partition = "gpu_diasfrazer",
        # No explicit qos: this cluster auto-derives it from the requested runtime.
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


rule predict_ref_gene_expression_full:
    """One-time reference-sequence forward pass (no per-sample looping) over the same full
    1,048,576bp window/gene set as predict_alphagenome_genexpr_full, to get each gene's
    baseline (reference-allele) RNA_seq expression level - never saved by that per-sample
    rule, which only persists ref-vs-alt log-fold-change scores. Needed to bin genes by
    baseline expression in notebooks/hbb.ipynb, independent of any individual's variants."""
    input:
        fasta = config["alphagenome_genexpr"]["fasta"],
        gtf = config["alphagenome_genexpr"]["gene_annotation"],
        weights = config["alphagenome_jax"]["paths"]["weights"],
        gene_ids_file = config["alphagenome_genexpr"]["gene_ids"],
    output:
        config["alphagenome_genexpr"]["paths"]["ref_gene_expression_full"],
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        tss = config["alphagenome_atlas"]["hbb_window"]["tss"],
        window_size = config["alphagenome_genexpr"]["full"]["window_size"],
    resources:
        # A single forward pass over the full 1,048,576bp window, not 3,202 of them like
        # predict_alphagenome_genexpr_full - minutes, not days - but still needs the same
        # H100 (80GB) to fit the window at all (verified there).
        runtime = 30,
        mem_mb = 32000,
        gres = "gpu:h100:1",
        partition = "gpu_diasfrazer",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/02-preprocess_data/scripts/ref_gene_expression.py \
            --fasta {input.fasta} \
            --gtf {input.gtf} \
            --weights-dir {input.weights} \
            --gene-ids-file {input.gene_ids_file} \
            --chromosome {params.chromosome} \
            --tss {params.tss} \
            --window-size {params.window_size} \
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
