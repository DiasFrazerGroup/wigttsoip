rule download_gencode_annotation:
    params:
        url = config["gencode"]["urls"]["gtf"],
    output:
        gtf = config["gencode"]["paths"]["gtf"]
    conda:
        "wihttsoip"
    shell:
        """
        wget --user-agent="Chrome" --no-check-certificate {params.url} -O {output.gtf}

        echo "Done!"
        """

rule gencode_gtf_to_parquet:
    input:
        gtf = config["gencode"]["paths"]["gtf"]
    output:
        parquet = config["gencode"]["paths"]["gtf_parquet"]
    run:
        import pyranges as pr
        gtf = pr.read_gtf(input.gtf)
        gtf.df.to_parquet(output.parquet, compression="zstd", index=False)

        print("Done!")
