rule severus:
    conda:
        "../../envs/severus.yaml"
    input:
        bam=ancient("clair3/{sample}/{sample}.sorted.haplotagged.bam"),
        vcf=ancient("clair3/{sample}/phased_merge_output.vcf.gz"),
        bed=config["bed_nvtr"],
    output:
        vcf=protected("severus/{sample}/all_SVs/severus_all.vcf"),
    params:
        dir=directory("severus/{sample}"),
        min_reads=config["min_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_size=config["min_size"],
    threads: 1
    log:
        "logs/{sample}/severus.log",
    shell:
        """
        severus \\
            --target-bam {input.bam} \\
            --out-dir {params.dir} \\
            -t {threads} \\
            --min-mapq {params.min_quality_mapping} \\
            --min-support {params.min_reads} \\
            --min-sv-size {params.min_size} \\
            --between-junction-ins \\
            --output-read-ids \\
            --phasing-vcf {input.vcf} \\
            --vntr-bed {input.bed} \\
            1> {log} 2>&1
        """


rule format_severus:
    conda:
        "../../envs/bcftools.yaml"
    input:
        vcf="severus/{sample}/all_SVs/severus_all.vcf",
    output:
        tab=temp("debreak/{sample}/rename.tab"),
        vcf=protected("severus/{sample}/severus.vcf"),
    params:
        header="{sample}.sorted.haplotagged\\t{sample}",
    log:
        "logs/{sample}/format_severus.log",
    shell:
        """
        {{ echo -e "{params.header}" > {output.tab}

        bcftools reheader -s {output.tab} {input.vcf} > {output.vcf}; }} \\
        1> {log} 2>&1
        """
