rule call_sv_debreak:
    input:
        bam="minimap2/{sample}/{sample}.sorted.bam",
        fasta=config["fasta"],
    output:
        dir_out=directory("debreak/{sample}/"),
        vcf=protected("debreak/{sample}/debreak.vcf"),
        vcf_tmp=temp("debreak/{sample}/debreak.vcf.tmp"),
        tab_rename=temp("debreak/{sample}/rename.tab"),
    params:
        min_num_reads=config["min_num_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_length_sv=config["min_length_sv"],
    threads: config["threads"]
    log:
        "logs/{sample}/call_sv_debreak.log",
    shell:
        """
        {{ debreak \\
        --thread {threads} \\
        --min_quality {params.min_quality_mapping} \\
        --bam {input.bam} \\
        --outpath {output.dir_out} \\
        --min_size {params.min_length_sv} \\
        --max_size 999999999999 \\
        --min_support {params.min_num_reads} \\
        --rescue_dup \\
        --rescue_large_ins \\
        --poa \\
        --tumor \\
        --ref {input.fasta}

        mv {output.vcf} {output.vcf_tmp}
        echo -e "{input.bam}\\t{wildcards.sample}" > {output.tab_rename}
        bcftools reheader -s {output.tab_rename} {output.vcf_tmp} > {output.vcf}

        echo -e "[INFO] debreak is done!"; }} \\
        1> {log} 2>&1
        """
