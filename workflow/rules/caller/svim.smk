rule svim:
    conda:
        "../../envs/svim.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
    output:
        vcf=protected("svim/{sample}/variants.vcf"),
        vcf_renamed="svim/{sample}/svim.vcf",
    params:
        dir=directory("svim/{sample}"),
        min_reads=config["min_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_size=config["min_size"],
    log:
        "logs/{sample}/svim.log",
    shell:
        """
        {{ svim alignment \\
            {params.dir} {input.bam} {input.fasta} \\
            --sample {wildcards.sample} \\
            --read_names \\
            --min_mapq {params.min_quality_mapping} \\
            --min_sv_size {params.min_size} \\
            --minimum_depth {params.min_reads}

        ln -r -s {output.vcf} {output.vcf_renamed}; }} \\
        1> {log} 2>&1
        """
