rule nanovar:
    conda:
        "../../envs/nanovar.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
    output:
        vcf=temp("nanovar/{sample}/{sample}.sorted.nanovar.pass.vcf"),
    params:
        dir=directory("nanovar/{sample}"),
        genome=config["genome"],
        min_reads=config["min_reads"],
        min_length_reads=config["min_length_reads"],
        min_size=config["min_size"],
    threads: 1
    log:
        "logs/{sample}/nanovar.log",
    shell:
        """
        {{ [ {params.genome} == GRCh37 ] && genome=hg19
        [ {params.genome} == GRCh38 ] && genome=hg38

        nanovar \\
            --threads {threads} \\
            --filter_bed ${{genome}} \\
            --data_type ont \\
            --minalign {params.min_length_reads} \\
            --mincov {params.min_reads} \\
            --minlen {params.min_size} \\
            --homo 0.75 \\
            --hetero 0.1 \\
            {input.bam} {input.fasta} {params.dir}; }} \\
        1> {log} 2>&1
        """


rule format_nanosv:
    conda:
        "../../envs/bcftools.yaml"
    input:
        vcf="nanovar/{sample}/{sample}.sorted.nanovar.pass.vcf",
    output:
        tab=temp("nanovar/{sample}/rename.tab"),
        vcf=protected("nanovar/{sample}/nanovar.vcf"),
    params:
        header="{sample}.sorted\\t{sample}",
    log:
        "logs/{sample}/format_nanosv.log",
    shell:
        """
        {{ echo -e "{params.header}" > {output.tab}
        bcftools reheader -s {output.tab} {input.vcf} > {output.vcf}; }} \\
        1> {log} 2>&1
        """
