rule call_sv_svim:
    input:
        bam=ancient("minimap2/{sample}/{sample}.sorted.bam"),
        fasta=config["fasta"],
    output:
        vcf_origin=protected("svim/{sample}/variants.vcf"),
        vcf="svim/{sample}/svim.vcf",
    params:
        dir_out=directory("svim/{sample}"),
        min_num_reads=config["min_num_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_length_sv=config["min_length_sv"],
    log:
        "logs/{sample}/call_sv_svim.log",
    shell:
        """
        {{ svim alignment \\
        {params.dir_out} {input.bam} {input.fasta} \\
        --sample {wildcards.sample} \\
        --read_names \\
        --min_mapq {params.min_quality_mapping} \\
        --min_sv_size {params.min_length_sv} \\
        --minimum_depth {params.min_num_reads}

        cd {params.dir_out}
        ln -s variants.vcf $(basename {output.vcf})

        echo -e "[INFO] SVIM is done!"; }} \\
        1> {log} 2>&1
        """
