rule debreak:
    conda:
        "../../envs/debreak.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
    output:
        vcf=temp("debreak/{sample}/tmp.vcf"),
    params:
        dir=directory("debreak/{sample}"),
        min_reads=config["min_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_size=config["min_size"],
    threads: 1
    log:
        "logs/{sample}/debreak.log",
    shell:
        """
        {{ debreak \\
            --thread {threads} \\
            --min_quality {params.min_quality_mapping} \\
            --bam {input.bam} \\
            --outpath {params.dir} \\
            --min_size {params.min_size} \\
            --max_size 999999999999 \\
            --min_support {params.min_reads} \\
            --rescue_dup \\
            --rescue_large_ins \\
            --poa \\
            --tumor \\
            --ref {input.fasta}

        mv {params.dir}/debreak.vcf {output.vcf}; }} \\
        1> {log} 2>&1
        """


rule format_debreak:
    conda:
        "../../envs/bcftools.yaml"
    input:
        vcf="debreak/{sample}/tmp.vcf",
    output:
        vcf=protected("debreak/{sample}/debreak.vcf"),
        tab=temp("debreak/{sample}/rename.tab"),
    params:
        header=f"{MAPPER}{{sample}}/{{sample}}.sorted.bam\\t{{sample}}",
    log:
        "logs/{sample}/format_debreak.log",
    shell:
        """
        {{ echo -e "{params.header}" > {output.tab}
        bcftools reheader -s {output.tab} {input.vcf} > {output.vcf}; }} \\
        1> {log} 2>&1
        """
