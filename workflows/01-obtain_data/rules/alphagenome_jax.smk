rule download_alphagenome_jax_weights:
    """Local JAX weights for running the actual AlphaGenome model (as opposed to the
    remote Atlas API), used by 02-preprocess_data's personalized-sequence forward pass."""
    params:
        weights = config["alphagenome_jax"]["urls"]["weights"],
    output:
        weights = directory(config["alphagenome_jax"]["paths"]["weights"]),
    conda:
        "wigttsoip"
    shell:
        """
        hf download {params.weights} --local-dir {output.weights}

        echo "Done!"
        """
