rule download_1000genomes_vcf:
    """1000 Genomes 30x high-coverage rerelease, raw unphased GT + annotations,
    one chromosome at a time (only chr11 is requested, scoped to the HBB
    window - see thousand_genomes.chromosome in config.yaml)."""
    params:
        url = lambda wildcards: config["thousand_genomes"]["urls"]["vcf"].format(chromosome=wildcards.chromosome),
    output:
        vcf = config["thousand_genomes"]["paths"]["vcf"],
        tbi = config["thousand_genomes"]["paths"]["vcf"] + ".tbi",
    conda:
        "wihttsoip"
    shell:
        """
        wget --user-agent="Chrome" --no-check-certificate {params.url} -O {output.vcf}
        tabix -p vcf {output.vcf}

        echo "Done!"
        """
