rule call_sv_cutesv:
    input:
        bam=ancient("minimap2/{sample}/{sample}.sorted.bam"),
        fasta=config["fasta"],
    output:
        vcf=protected("cutesv/{sample}/cutesv.vcf"),
    params:
        dir_out=directory("cutesv/{sample}"),
        min_length_reads=config["min_length_reads"],
        min_num_reads=config["min_num_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_length_sv=config["min_length_sv"],
    threads: math.floor(workflow.cores / NUMBER_SAMPLES / NUMBER_CALLERS)
    log:
        "logs/{sample}/call_sv_cutesv.log",
    shell:
        """
        {{ cuteSV \\
        {input.bam} {input.fasta} {output.vcf} {params.dir_out} \\
        --sample {wildcards.sample} \\
        --threads {threads} \\
        --max_cluster_bias_INS 100 \\
        --diff_ratio_merging_INS 0.3 \\
        --max_cluster_bias_DEL 100 \\
        --diff_ratio_merging_DEL 0.3 \\
        --min_support {params.min_num_reads} \\
        --min_read_len {params.min_length_reads} \\
        --min_mapq {params.min_quality_mapping} \\
        --min_size {params.min_length_sv} \\
        --max_size -1 \\
        --genotype \\
        --report_readid

        echo -e "[INFO] cuteSV is done!"; }} \\
        1> {log} 2>&1
        """
