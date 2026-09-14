rule render_hbb_notebook:
    """Execute notebooks/hbb.ipynb in place via jupyter nbconvert - gene annotation around
    HBB (1Mb overview + zoom), variant burden per individual, combinatorial-vs-summed effect
    distributions, a combined-vs-sum-of-parts scatter, effective-dimensionality plots
    (correlation vs. number of top variants summed, plus the same distribution/scatter pair
    restricted to just the top 3 variants), a correlation-vs-k comparison across every gene
    in the window (scripts/topk_correlation_per_gene.py's output), and the same comparison
    faceted into low/medium/high tertiles of each gene's baseline (reference) expression
    (scripts/ref_gene_expression.py's output). Figures are also saved as PDF under
    notebooks/pdfs/hbb/ (see notebooks/figutils.py, vendored from
    ../alphagenome_finetuning_rna/figures/figutils.py for consistent styling). Reads
    full-scale outputs (SCALE = "full" inside the notebook) - flip to "dev" if iterating on
    the notebook itself again. Needs unique_variant_scores (per-variant singles scores, for
    the effective-dimensionality plots), topk_correlation_per_gene_full, and
    ref_gene_expression_full in addition to whole_vs_sum_parts and the GTF.

    The rule's declared `output:` is a completion marker, NOT the notebook itself: Snakemake
    deletes a rule's declared outputs before running its shell command (to detect failed
    jobs reliably), and this rule's shell command executes the notebook IN PLACE - if the
    notebook were the declared output, Snakemake would delete it right before nbconvert
    tried to open it (hit this directly: "pattern 'notebooks/hbb.ipynb' matched no files").
    notebooks/hbb.ipynb is listed as an `input:` instead, so its own prose/plotting-code
    edits (source changes) still trigger a rerun, but Snakemake never deletes it."""
    input:
        notebook = "notebooks/hbb.ipynb",
        whole_vs_sum_parts = config["analysis"]["paths"]["whole_vs_sum_parts_full"],
        unique_variant_scores = config["analysis"]["paths"]["unique_variant_scores_full"],
        topk_correlation_per_gene = config["analysis"]["paths"]["topk_correlation_per_gene_full"],
        ref_gene_expression = config["alphagenome_genexpr"]["paths"]["ref_gene_expression_full"],
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
