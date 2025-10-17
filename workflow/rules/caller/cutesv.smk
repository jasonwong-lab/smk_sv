rule cutesv:
    conda:
        "../../envs/cutesv.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
    output:
        vcf=protected("cutesv/{sample}/cutesv.vcf"),
    params:
        dir=directory("cutesv/{sample}"),
        min_length_reads=config["min_length_reads"],
        min_reads=config["min_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_size=config["min_size"],
        max_size=config["max_size"],
    threads: 1
    log:
        "logs/{sample}/cutesv.log",
    shell:
        """
        cuteSV \\
            {input.bam} {input.fasta} {output.vcf} {params.dir} \\
            --sample {wildcards.sample} \\
            --threads {threads} \\
            --max_cluster_bias_INS 100 \\
            --diff_ratio_merging_INS 0.3 \\
            --max_cluster_bias_DEL 100 \\
            --diff_ratio_merging_DEL 0.3 \\
            --min_support {params.min_reads} \\
            --min_read_len {params.min_length_reads} \\
            --min_mapq {params.min_quality_mapping} \\
            --min_size {params.min_size} \\
            --max_size {params.max_size} \\
            --genotype \\
            --report_readid \\
            1> {log} 2>&1
        """
