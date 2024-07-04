rule call_sv_delly:
    container:
        "docker://mhjiang97/smk_sv:latest"
    input:
        bam=ancient("minimap2/{sample}/{sample}.sorted.bam"),
        fasta=config["fasta"],
    output:
        vcf=protected("delly/{sample}/delly.vcf"),
        bcf=protected("delly/{sample}/delly.bcf"),
        tab_rename=temp("debreak/{sample}/rename.tab"),
    params:
        min_quality_mapping=config["min_quality_mapping"],
    log:
        "logs/{sample}/call_sv_delly.log",
    shell:
        """
        {{ delly lr \\
        -y ont \\
        -o {output.bcf} \\
        -g {input.fasta} \\
        -q {params.min_quality_mapping} \\
        {input.bam}

        BAM={input.bam}
        echo -e "$(basename ${{BAM%.*}})\\t{wildcards.sample}" > {output.tab_rename}
        bcftools view -O v {output.bcf} | bcftools reheader -s {output.tab_rename} > {output.vcf}

        echo -e "[INFO] delly is done!"; }} \\
        1> {log} 2>&1
        """
