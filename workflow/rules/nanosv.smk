rule call_sv_nanosv:
    input:
        bam="minimap2/{sample}/{sample}.sorted.bam",
        bed_nanosv=config["bed_nanosv"],
        config_nanosv=config["config_nanosv"],
    output:
        vcf=protected("nanosv/{sample}/nanosv.vcf"),
        vcf_tmp=temp("nanosv/{sample}/nanosv.vcf.tmp"),
        tab_rename=temp("nanosv/{sample}/rename.tab"),
    threads: config["threads"]
    log:
        "logs/{sample}/call_sv_nanosv.log",
    shell:
        """
        {{ NanoSV -t {threads} -s samtools -c {input.config_nanosv} -b {input.bed_nanosv} -o {output.vcf} {input.bam}
        mv {output.vcf} {output.vcf_tmp}
        BAM={input.bam}
        echo -e "${{BAM%%.*}}\\t{wildcards.sample}" > {output.tab_rename}
        bcftools reheader -s {output.tab_rename} {output.vcf_tmp} > {output.vcf}; }} \\
        1> {log} 2>&1
        """
