rule delly:
    conda:
        "../../envs/delly.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
    output:
        bcf=temp("delly/{sample}/delly.bcf"),
    params:
        min_quality_mapping=config["min_quality_mapping"],
    log:
        "logs/{sample}/delly.log",
    shell:
        """
        delly lr \\
            -y ont \\
            -o {output.bcf} \\
            -g {input.fasta} \\
            -q {params.min_quality_mapping} \\
            {input.bam} \\
            1> {log} 2>&1
        """


rule format_delly:
    conda:
        "../../envs/bcftools.yaml"
    input:
        bcf="delly/{sample}/delly.bcf",
    output:
        tab=temp("debreak/{sample}/rename.tab"),
        vcf=protected("delly/{sample}/delly.vcf"),
    params:
        header="{sample}.sorted\\t{sample}",
    log:
        "logs/{sample}/format_delly.log",
    shell:
        """
        {{ echo -e "{params.header}" > {output.tab}
        bcftools view -O v {input.bcf} | \\
            bcftools reheader -s {output.tab} \\
            > {output.vcf}; }} \\
        1> {log} 2>&1
        """
