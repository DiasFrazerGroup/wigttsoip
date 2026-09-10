rule render_hbb_notebook:
    """Execute notebooks/hbb.ipynb in place via jupyter nbconvert - a simple scatter-plot
    scaffold for the whole-vs-sum-of-parts analysis (single-variant effects across the HBB
    window, single-variant effects vs. gnomAD MAF, and combinatorial vs. sum-of-parts).
    Reads full-scale outputs (SCALE = "full" inside the notebook, matching the input paths
    below) - flip both back to "dev" if iterating on the notebook itself again.

    The rule's declared `output:` is a completion marker, NOT the notebook itself: Snakemake
    deletes a rule's declared outputs before running its shell command (to detect failed
    jobs reliably), and this rule's shell command executes the notebook IN PLACE - if the
    notebook were the declared output, Snakemake would delete it right before nbconvert
    tried to open it (hit this directly: "pattern 'notebooks/hbb.ipynb' matched no files").
    notebooks/hbb.ipynb is listed as an `input:` instead, so its own prose/plotting-code
    edits (source changes) still trigger a rerun, but Snakemake never deletes it."""
    input:
        notebook = "notebooks/hbb.ipynb",
        unique_variant_scores = config["analysis"]["paths"]["unique_variant_scores_full"],
        unique_variant_maf = config["analysis"]["paths"]["unique_variant_gnomad_maf_full"],
        whole_vs_sum_parts = config["analysis"]["paths"]["whole_vs_sum_parts_full"],
        gtf = config["gencode"]["paths"]["gtf_parquet"],
    output:
        touch("notebooks/.hbb_rendered"),
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
        python -m ipykernel install --user --name wigttsoip --display-name wigttsoip

        jupyter nbconvert --to notebook --execute --inplace \
            --ExecutePreprocessor.kernel_name=wigttsoip \
            --ExecutePreprocessor.timeout=600 \
            {input.notebook}

        echo "Done!"
        """
