rule list_hbb_window_genes:
    """Bare Ensembl gene IDs of every gene overlapping the full HBB window (a superset of the
    small dev window too - used as-is for both, since the merge-time gene/position filter in
    download_alphagenome_atlas_genexpr.py drops anything that isn't actually present)."""
    input:
        gtf = config["gencode"]["paths"]["gtf_parquet"],
    output:
        config["alphagenome_atlas"]["paths"]["gene_ids"],
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["full"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["full"]["end"],
    resources:
        runtime = 30,
        mem_mb = 4000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    run:
        import pandas as pd

        gtf = pd.read_parquet(
            input.gtf, columns=["Chromosome", "Feature", "Start", "End", "gene_id"]
        )
        chrom_bare = params.chromosome.replace("chr", "", 1)
        overlaps = (gtf["Start"] < params.end) & (gtf["End"] > params.start)
        genes = gtf.loc[
            (gtf["Feature"] == "gene")
            & gtf["Chromosome"].isin([params.chromosome, chrom_bare])
            & overlaps,
            "gene_id",
        ].str.split(".").str[0].unique()

        with open(output[0], "w") as f:
            f.write("\n".join(genes))

        print(f"Found {len(genes)} genes in {params.chromosome}:{params.start}-{params.end}")
        print("Done!")


rule download_alphagenome_atlas_hbb_tss_pm5kb_genexpr:
    """Dev-scale download: AlphaGenome Atlas variant effect scores on gene expression (all
    RNA_SEQ tracks) for every single-nucleotide variant in a 10kb window (HBB's TSS +/- 5kb),
    for every gene overlapping the full HBB window. Meant for developing/testing the downstream
    pipeline quickly; shares its chunk cache with download_alphagenome_atlas_hbb_window_genexpr,
    so running the full window afterwards never re-downloads these chunks. Output is just an
    empty marker, not a merged h5ad - see scripts/extract_alphagenome_atlas_genexpr.py to read data."""
    input:
        gene_ids_file = config["alphagenome_atlas"]["paths"]["gene_ids"],
        gene_annotation = config["gencode"]["paths"]["gtf_parquet"],
    output:
        config["alphagenome_atlas"]["paths"]["done_small"],
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["small"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["small"]["end"],
        chunk_size = config["alphagenome_atlas"]["chunk_size_bp"],
        chunk_cache_dir = config["alphagenome_atlas"]["paths"]["chunk_cache"],
    resources:
        runtime = 60,
        mem_mb = 16000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    threads: 5
    conda:
        "wigttsoip"
    retries: 3  # only unfetched chunks are re-queried on retry/resume - see chunk_cache_dir
    shell:
        """
        python workflows/01-obtain_data/scripts/download_alphagenome_atlas_genexpr.py \
            --interval {params.chromosome}:{params.start}-{params.end} \
            --gene-ids-file {input.gene_ids_file} \
            --gene-annotation {input.gene_annotation} \
            --output {output} \
            --chunk-size {params.chunk_size} \
            --chunk-cache-dir {params.chunk_cache_dir} \
            --max-workers {threads}

        echo "Done!"
        """


rule extract_hbb_wholeblood_variant_effects_long:
    """Example downstream extraction: how each variant in the HBB TSS+/-5kb window alters
    HBB's own expression specifically in GTEx whole blood tracks, as a long-format
    (variant, gene, track, score, quantile) parquet. Reads only the matching row/column
    slice of each overlapping chunk via h5py (scripts/extract_alphagenome_atlas_genexpr.py),
    never a whole chunk's full matrix."""
    input:
        done = config["alphagenome_atlas"]["paths"]["done_small"],
        gene_annotation = config["gencode"]["paths"]["gtf_parquet"],
    output:
        directory(config["alphagenome_atlas"]["paths"]["hbb_wholeblood_long"]),
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["small"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["small"]["end"],
        chunk_size = config["alphagenome_atlas"]["chunk_size_bp"],
        chunk_cache_dir = config["alphagenome_atlas"]["paths"]["chunk_cache"],
        gene_ids = ["ENSG00000244734"],  # HBB
        biosample_name = "venous blood",
    resources:
        runtime = 30,
        mem_mb = 4000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/01-obtain_data/scripts/extract_alphagenome_atlas_genexpr.py \
            --interval {params.chromosome}:{params.start}-{params.end} \
            --chunk-cache-dir {params.chunk_cache_dir} \
            --chunk-size {params.chunk_size} \
            --gene-ids {params.gene_ids} \
            --gene-annotation {input.gene_annotation} \
            --biosample-name "{params.biosample_name}" \
            --output {output}

        echo "Done!"
        """


rule extract_hbb_wholeblood_variant_effects_long_dev:
    """Same extraction as extract_hbb_wholeblood_variant_effects_long, but over hbb_window.dev
    (TSS +/- 8192bp) - matches alphagenome_genexpr.paths.dev's window exactly, so the
    dev-scale 03-analysis join can match every variant in a dev individual's combination
    (hbb_window.small, TSS +/- 5kb, is narrower and leaves most of that window's variants
    unmatched). Depends on done_full (not done_small): the small download's chunk cache
    doesn't cover the extra +/-3192bp this window needs, but the full-window download's
    does (same shared chunk_cache)."""
    input:
        done = config["alphagenome_atlas"]["paths"]["done_full"],
        gene_annotation = config["gencode"]["paths"]["gtf_parquet"],
    output:
        directory(config["alphagenome_atlas"]["paths"]["hbb_wholeblood_long_dev"]),
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["dev"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["dev"]["end"],
        chunk_size = config["alphagenome_atlas"]["chunk_size_bp"],
        chunk_cache_dir = config["alphagenome_atlas"]["paths"]["chunk_cache"],
        gene_ids = ["ENSG00000244734"],  # HBB
        biosample_name = "venous blood",
    resources:
        runtime = 30,
        mem_mb = 4000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/01-obtain_data/scripts/extract_alphagenome_atlas_genexpr.py \
            --interval {params.chromosome}:{params.start}-{params.end} \
            --chunk-cache-dir {params.chunk_cache_dir} \
            --chunk-size {params.chunk_size} \
            --gene-ids {params.gene_ids} \
            --gene-annotation {input.gene_annotation} \
            --biosample-name "{params.biosample_name}" \
            --output {output}

        echo "Done!"
        """


rule extract_hbb_wholeblood_variant_effects_long_full:
    """Same extraction as extract_hbb_wholeblood_variant_effects_long, but over every variant in
    the full 1,048,576bp HBB window rather than just the TSS+/-5kb dev window - reads only
    the matching row/column slice of each overlapping chunk via h5py, never a whole chunk's
    full matrix."""
    input:
        done = config["alphagenome_atlas"]["paths"]["done_full"],
        gene_annotation = config["gencode"]["paths"]["gtf_parquet"],
    output:
        directory(config["alphagenome_atlas"]["paths"]["hbb_wholeblood_long_full"]),
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["full"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["full"]["end"],
        chunk_size = config["alphagenome_atlas"]["chunk_size_bp"],
        chunk_cache_dir = config["alphagenome_atlas"]["paths"]["chunk_cache"],
        gene_ids = ["ENSG00000244734"],  # HBB
        biosample_name = "venous blood",
    resources:
        # ~2048 chunk files (1,048,576bp / 512bp chunk_size_bp) - parallelized across
        # `threads` (see --max-workers) instead of one serial process, which previously
        # both timed out and was near-OOM. Peak memory is bounded by chunks_per_batch
        # (default 50), not the whole interval, since extract_to_batches() flushes each
        # full batch to its own parquet immediately.
        runtime = 30,
        mem_mb = 24000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    threads: 32
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/01-obtain_data/scripts/extract_alphagenome_atlas_genexpr.py \
            --interval {params.chromosome}:{params.start}-{params.end} \
            --chunk-cache-dir {params.chunk_cache_dir} \
            --chunk-size {params.chunk_size} \
            --gene-ids {params.gene_ids} \
            --gene-annotation {input.gene_annotation} \
            --biosample-name "{params.biosample_name}" \
            --max-workers {threads} \
            --output {output}

        echo "Done!"
        """


rule extract_allgenes_variant_effects_long_full:
    """Same as extract_hbb_wholeblood_variant_effects_long_full, but keeping every gene in
    gene_ids.txt instead of filtering to just HBB - the raw chunk cache already has every
    gene overlapping the window (download_alphagenome_atlas_hbb_window_genexpr fetches
    every RNA_SEQ track for every gene, not just HBB), so this needs no new download, just
    a broader --gene-ids filter read from the same gene list alphagenome_genexpr.py's GPU
    forward pass already scores against."""
    input:
        done = config["alphagenome_atlas"]["paths"]["done_full"],
        gene_annotation = config["gencode"]["paths"]["gtf_parquet"],
        gene_ids_file = config["alphagenome_atlas"]["paths"]["gene_ids"],
    output:
        directory(config["alphagenome_atlas"]["paths"]["allgenes_long_full"]),
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["full"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["full"]["end"],
        chunk_size = config["alphagenome_atlas"]["chunk_size_bp"],
        chunk_cache_dir = config["alphagenome_atlas"]["paths"]["chunk_cache"],
        biosample_name = "venous blood",
    resources:
        # All 94 genes (not just HBB) means ~94x bigger row-explosion per chunk in this
        # gene-dense region - HBB-only memory settings OOM'd here. 128000MB/32 workers
        # gives headroom over the run that completed memory-safely at 16 workers; runtime
        # kept generous as a margin.
        runtime = 150,
        mem_mb = 128000,
        gres = "none",
        partition = "genoa64",
        qos = "short",
    threads: 32
    conda:
        "wigttsoip"
    shell:
        """
        python workflows/01-obtain_data/scripts/extract_alphagenome_atlas_genexpr.py \
            --interval {params.chromosome}:{params.start}-{params.end} \
            --chunk-cache-dir {params.chunk_cache_dir} \
            --chunk-size {params.chunk_size} \
            --gene-ids $(cat {input.gene_ids_file}) \
            --gene-annotation {input.gene_annotation} \
            --biosample-name "{params.biosample_name}" \
            --max-workers {threads} \
            --chunks-per-batch 10 \
            --output {output}

        echo "Done!"
        """


rule download_alphagenome_atlas_hbb_window_genexpr:
    """Download AlphaGenome Atlas variant effect scores on gene expression (all RNA_SEQ tracks)
    for every single-nucleotide variant in a 1,048,576 bp window (AlphaGenome's largest
    supported model input context) centered on HBB's TSS, for every gene overlapping that
    window. Shares its chunk cache with download_alphagenome_atlas_hbb_tss_pm5kb_genexpr, so
    chunks already fetched by that dev-scale run are never re-downloaded here. Output is just
    an empty marker, not a merged h5ad (would be hundreds of GB) - see
    scripts/extract_alphagenome_atlas_genexpr.py to read a gene/region's data on demand."""
    input:
        gene_ids_file = config["alphagenome_atlas"]["paths"]["gene_ids"],
        gene_annotation = config["gencode"]["paths"]["gtf_parquet"],
    output:
        config["alphagenome_atlas"]["paths"]["done_full"],
    params:
        chromosome = config["alphagenome_atlas"]["hbb_window"]["chromosome"],
        start = config["alphagenome_atlas"]["hbb_window"]["full"]["start"],
        end = config["alphagenome_atlas"]["hbb_window"]["full"]["end"],
        chunk_size = config["alphagenome_atlas"]["chunk_size_bp"],
        chunk_cache_dir = config["alphagenome_atlas"]["paths"]["chunk_cache"],
    resources:
        runtime = int(2*24*60),  # 2 days in minutes: generous ceiling for Atlas's real request-rate quota
        mem_mb = 24000,
        gres = "none",
        partition = "genoa64",
        qos = "marathon",  # this cluster's 7-day-max QOS; no GPU needed, network/quota-bound job
    threads: 5
    conda:
        "wigttsoip"
    retries: 3  # only unfetched chunks are re-queried on retry/resume - see chunk_cache_dir
    shell:
        """
        python workflows/01-obtain_data/scripts/download_alphagenome_atlas_genexpr.py \
            --interval {params.chromosome}:{params.start}-{params.end} \
            --gene-ids-file {input.gene_ids_file} \
            --gene-annotation {input.gene_annotation} \
            --output {output} \
            --chunk-size {params.chunk_size} \
            --chunk-cache-dir {params.chunk_cache_dir} \
            --max-workers {threads}

        echo "Done!"
        """
