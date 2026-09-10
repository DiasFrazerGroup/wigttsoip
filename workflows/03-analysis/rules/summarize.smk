rule get_unique_hbb_k562_variants_dev:
    """Distinct SNVs actually carried by the dev sample list in the HBB dev window - see
    scripts/get_unique_variants.py. Independent of any track/score data, so it feeds both
    subset_hbb_k562_unique_variant_scores_dev and
    annotate_hbb_k562_unique_variant_gnomad_maf_dev without either depending on the other."""
    input:
        vcf = config["alphagenome_genexpr"]["vcf"],
    output:
        config["analysis"]["paths"]["unique_variants_dev"],
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["dev"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["dev"]["end"],
        samples = config["alphagenome_genexpr"]["samples_dev"],
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
        python workflows/03-analysis/scripts/get_unique_variants.py \
            --vcf {input.vcf} \
            --chromosome {params.chromosome} \
            --start {params.start} \
            --end {params.end} \
            --samples {params.samples} \
            --output {output}

        echo "Done!"
        """


rule get_unique_hbb_k562_variants_full:
    """Same as get_unique_hbb_k562_variants_dev, but over the full 1,048,576bp window and
    the full 3,202-sample cohort (no --samples filter needed: a jointly-called multi-sample
    VCF only lists a site because >=1 of ALL its samples carries the ALT there, so every SNV
    row already qualifies)."""
    input:
        vcf = config["alphagenome_genexpr"]["vcf"],
    output:
        config["analysis"]["paths"]["unique_variants_full"],
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["full"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["full"]["end"],
    resources:
        runtime = 30,
        mem_mb = 8000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/03-analysis/scripts/get_unique_variants.py \
            --vcf {input.vcf} \
            --chromosome {params.chromosome} \
            --start {params.start} \
            --end {params.end} \
            --output {output}

        echo "Done!"
        """


rule subset_hbb_k562_unique_variant_scores_dev:
    """Restrict Atlas's HBB/K562 singles scores (every possible alt allele at every
    position) to just the unique variant set from get_unique_hbb_k562_variants_dev - see
    scripts/subset_unique_variant_scores.py. Makes the whole_vs_sum_parts join below run
    against a much smaller table."""
    input:
        variants = config["analysis"]["paths"]["unique_variants_dev"],
        singles = config["alphagenome_atlas"]["paths"]["hbb_k562_long_dev"],
    output:
        config["analysis"]["paths"]["unique_variant_scores_dev"],
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
        python workflows/03-analysis/scripts/subset_unique_variant_scores.py \
            --variants {input.variants} \
            --singles {input.singles} \
            --output {output}

        echo "Done!"
        """


rule subset_hbb_k562_unique_variant_scores_full:
    """Same as subset_hbb_k562_unique_variant_scores_dev, but over the full window's unique
    variant set."""
    input:
        variants = config["analysis"]["paths"]["unique_variants_full"],
        singles = config["alphagenome_atlas"]["paths"]["hbb_k562_long_full"],
    output:
        config["analysis"]["paths"]["unique_variant_scores_full"],
    resources:
        runtime = 30,
        mem_mb = 8000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/03-analysis/scripts/subset_unique_variant_scores.py \
            --variants {input.variants} \
            --singles {input.singles} \
            --output {output}

        echo "Done!"
        """


rule annotate_hbb_k562_unique_variant_gnomad_maf_dev:
    """Annotate the dev-scale unique-variant set with gnomAD v3.1.1 allele frequency and
    minor allele frequency - see scripts/annotate_variants_gnomad_maf.py (one indexed
    region fetch of gnomAD's tabix'd VCF, then in-memory lookups, not one query/variant).
    Depends only on get_unique_hbb_k562_variants_dev, not on any Atlas/track data."""
    input:
        variants = config["analysis"]["paths"]["unique_variants_dev"],
        gnomad_vcf = config["gnomad"]["paths"]["chr11_vcf"],
    output:
        config["analysis"]["paths"]["unique_variant_gnomad_maf_dev"],
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["dev"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["dev"]["end"],
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
        python workflows/03-analysis/scripts/annotate_variants_gnomad_maf.py \
            --variants {input.variants} \
            --gnomad-vcf {input.gnomad_vcf} \
            --chromosome {params.chromosome} \
            --start {params.start} \
            --end {params.end} \
            --output {output}

        echo "Done!"
        """


rule annotate_hbb_k562_unique_variant_gnomad_maf_full:
    """Same as annotate_hbb_k562_unique_variant_gnomad_maf_dev, but over the full window's
    unique-variant set."""
    input:
        variants = config["analysis"]["paths"]["unique_variants_full"],
        gnomad_vcf = config["gnomad"]["paths"]["chr11_vcf"],
    output:
        config["analysis"]["paths"]["unique_variant_gnomad_maf_full"],
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["full"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["full"]["end"],
    resources:
        # A single fetch across the full 1,048,576bp window took ~34s/~296k gnomAD
        # records in direct testing - generous headroom, not a tight fit.
        runtime = 30,
        mem_mb = 8000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/03-analysis/scripts/annotate_variants_gnomad_maf.py \
            --variants {input.variants} \
            --gnomad-vcf {input.gnomad_vcf} \
            --chromosome {params.chromosome} \
            --start {params.start} \
            --end {params.end} \
            --output {output}

        echo "Done!"
        """


rule summarize_whole_vs_sum_parts_dev:
    """Join, per sample and track, the model's predicted effect of an individual's whole
    combination of HBB-window SNVs (genexpr_personalized, "whole_score") against the sum
    of those same SNVs' own single-variant Atlas effects (alphagenome_atlas,
    "sum_parts_score") - the additive/no-epistasis null. See
    scripts/summarize_whole_vs_sum_parts.py for the join logic (tracks keyed by
    (track_name, track_strand), no averaging across distinct tracks). Joins against the
    unique-variant-subsetted singles table (see subset_hbb_k562_unique_variant_scores_dev),
    not the raw Atlas extraction, for a smaller/faster join."""
    input:
        combined = config["alphagenome_genexpr"]["paths"]["hbb_k562_dev"],
        singles = config["analysis"]["paths"]["unique_variant_scores_dev"],
    output:
        config["analysis"]["paths"]["whole_vs_sum_parts_dev"],
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
        python workflows/03-analysis/scripts/summarize_whole_vs_sum_parts.py \
            --combined {input.combined} \
            --singles {input.singles} \
            --output {output} \
            --memory-limit-mb {resources.mem_mb}

        echo "Done!"
        """


rule summarize_whole_vs_sum_parts_full:
    """Same join as summarize_whole_vs_sum_parts_dev, but over all 3,202 samples and the
    full 1,048,576bp window on both sides.

    Without an explicit memory_limit, duckdb sizes its buffer pool off the node's total
    physical RAM, not this job's actual SLURM/cgroup --mem allocation - under-constrained,
    it over-allocates and gets OS OOM-killed rather than spilling to disk gracefully (hit
    this directly: the exploded-variant join over the full 1,048,576bp window's per-sample
    variant lists is memory-heavy well before the final row counts are). --memory-limit-mb
    is always passed matching resources.mem_mb for exactly this reason."""
    input:
        combined = config["alphagenome_genexpr"]["paths"]["hbb_k562_full"],
        singles = config["analysis"]["paths"]["unique_variant_scores_full"],
    output:
        config["analysis"]["paths"]["whole_vs_sum_parts_full"],
    resources:
        runtime = 60,
        mem_mb = 32000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/03-analysis/scripts/summarize_whole_vs_sum_parts.py \
            --combined {input.combined} \
            --singles {input.singles} \
            --output {output} \
            --memory-limit-mb {resources.mem_mb}

        echo "Done!"
        """
