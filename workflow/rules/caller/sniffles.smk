rule sniffles:
    conda:
        "../../envs/sniffles.yaml"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
        bed=config["bed_tandem_repeats"],
    output:
        vcf=protected("sniffles/{sample}/sniffles.vcf"),
    params:
        min_length_reads=config["min_length_reads"],
        min_reads=config["min_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_size=config["min_size"],
    threads: 1
    log:
        "logs/{sample}/sniffles.log",
    shell:
        """
        sniffles \\
            --input {input.bam} \\
            --vcf {output.vcf} \\
            --sample-id {wildcards.sample} \\
            --threads {threads} \\
            --reference {input.fasta} \\
            --tandem-repeats {input.bed} \\
            --minsupport {params.min_reads} \\
            --mapq {params.min_quality_mapping} \\
            --minsvlen {params.min_size} \\
            --min-alignment-length {params.min_length_reads} \\
            --output-rnames \\
            --allow-overwrite \\
            1> {log} 2>&1
        """
