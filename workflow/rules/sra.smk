# Rule: Prefetch SRA data
rule sra_prefetch:
    output:
        temp("sra-data/{SRA}/{SRA}.sra")
    conda:
        "../envs/sra.yaml"
    shell:
        """
        prefetch -O sra-data {wildcards.SRA}
        """

# Rule: Convert SRA to FASTQ using parallel-fastq-dump
rule parallel_fastq_dump:
    input:
        "sra-data/{srr}/{srr}.sra"
    output:
        r1 = "sra-data/{srr}_1.fastq.gz",
        r2 = "sra-data/{srr}_2.fastq.gz"
    conda:
        "../envs/sra.yaml"
    threads:
        8
    run:
        # Determine library type (Single or Paired)
        lib = samples.loc[samples["Fastq1"] == wildcards.srr, "Library"].unique()[0]

        # Run parallel-fastq-dump based on library type
        if lib == "Single":
            shell("""
                parallel-fastq-dump -t {threads} --split-files --gzip -s {input} -O sra-data
                touch {output.r2}  # Ensure r2 exists even for Single-end libraries
            """)
        elif lib == "Paired":
            shell("""
                parallel-fastq-dump -t {threads} --split-files --gzip -s {input} -O sra-data
            """)