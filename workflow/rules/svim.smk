rule call_sv_svim:
    input:
        bam="minimap2/{sample}/{sample}.sorted.bam",
        fasta=config["fasta"],
    output:
        dir_out=protected(directory("svim/{sample}/")),
        vcf=protected("svim/{sample}/svim.vcf"),
    params:
        min_num_reads=config["min_num_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_length_sv=config["min_length_sv"],
    log:
        "logs/{sample}/call_sv_svim.log",
    shell:
        """
        {{ svim alignment \\
        {output.dir_out} {input.bam} {input.fasta} \\
        --sample {wildcards.sample} \\
        --read_names \\
        --min_mapq {params.min_quality_mapping} \\
        --min_sv_size {params.min_length_sv} \\
        --minimum_depth {params.min_num_reads} && \\

        cd {output.dir_out} \\
            && ln -s variants.vcf $(basename {output.vcf})

        echo -e "[INFO] SVIM is done!"; }} \\
        1> {log} 2>&1
        """
