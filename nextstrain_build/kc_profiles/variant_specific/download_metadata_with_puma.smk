

rule download_metadata_s3:
    message: "Downloading metadata from {params.address} -> {output.metadata}"
    output:
        metadata = "data/downloaded_meta.tsv"
    conda: config["conda_environment"]
    params:
        address = "s3://nextstrain-ncov-private/metadata.tsv.gz"
    shell:
        """
        aws s3 cp {params.address} - | gunzip -cq >{output.metadata:q}
        """

rule add_puma:
    message:
        """
        adds King County PUMA to metadata. If not in King County, PUMA = Washington or Region if != Washington
        """
    input:
        metadata = rules.download_metadata_s3.output,
        puma_raw = "data/fullkc_puma.tsv"
    output:
        # Note: the command doesn't use these, but adding them here makes snakemake
        # aware that this rule produces them
        kc_meta = "data/metadata_kc_puma.tsv",

    shell:
        """
        python3 scripts/puma_kc_seq_merge.py \
        --metadata {input.metadata} \
        --puma_raw {input.puma_raw} \
        --output {output.kc_meta}

        """
