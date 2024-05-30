rule call_sv_sniffles:
    input:
        bam=ancient("minimap2/{sample}/{sample}.sorted.bam"),
        fasta=config["fasta"],
        bed_tandem_repeats=config["bed_tandem_repeats"],
    output:
        vcf=protected("sniffles/{sample}/sniffles.vcf"),
    params:
        min_length_reads=config["min_length_reads"],
        min_num_reads=config["min_num_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_length_sv=config["min_length_sv"],
    threads: config["threads"]
    log:
        "logs/{sample}/call_sv_sniffles.log",
    shell:
        """
        {{ sniffles \\
        --input {input.bam} \\
        --vcf {output.vcf} \\
        --sample-id {wildcards.sample} \\
        --threads {threads} \\
        --reference {input.fasta} \\
        --tandem-repeats {input.bed_tandem_repeats} \\
        --minsupport {params.min_num_reads} \\
        --mapq {params.min_quality_mapping} \\
        --minsvlen {params.min_length_sv} \\
        --min-alignment-length {params.min_length_reads} \\
        --output-rnames \\
        --allow-overwrite

        echo -e "[INFO] Sniffles is done!"; }} \\
        1> {log} 2>&1
        """
