rule call_sv_nanovar:
    input:
        bam="minimap2/{sample}/{sample}.sorted.bam",
        fasta=config["fasta"],
    output:
        dir_out=directory("nanovar/{sample}/"),
        vcf=protected("nanovar/{sample}/nanovar.vcf"),
        tab_rename=temp("nanovar/{sample}/rename.tab"),
    params:
        genome=config["genome"],
        min_num_reads=config["min_num_reads"],
        min_length_reads=config["min_length_reads"],
        min_length_sv=config["min_length_sv"],
    threads: config["threads"]
    log:
        "logs/{sample}/call_sv_nanovar.log",
    shell:
        """
        {{ [ {params.genome} == GRCh37 ] && genome=hg19
        [ {params.genome} == GRCh38 ] && genome=hg38
        BAM={input.bam}
        nanovar \\
        --threads {threads} \\
        --filter_bed ${{genome}} \\
        --data_type ont \\
        --minalign {params.min_length_reads} \\
        --mincov {params.min_num_reads} \\
        --minlen {params.min_length_sv} \\
        --homo 0.75 \\
        --hetero 0.1 \\
        {input.bam} {input.fasta} {output.dir_out}
        echo -e "$(basename ${{BAM%.*}})\\t{wildcards.sample}" > {output.tab_rename}
        bcftools reheader -s {output.tab_rename} {output.dir_out}/"$(basename ${{BAM%.*}})".nanovar.pass.vcf > {output.vcf}; }} \\
        1> {log} 2>&1
        """
