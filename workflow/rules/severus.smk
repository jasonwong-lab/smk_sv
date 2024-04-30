rule phasenhaplotag_bam_clair3:
    input:
        bam="minimap2/{sample}/{sample}.sorted.bam",
        fasta=config["fasta"],
        model_clair3=config["model_clair3"],
    output:
        dir_out=directory("clair3/{sample}"),
        vcf_phased=protected("clair3/{sample}/phased_merge_output.vcf.gz"),
        bam_haplotagged=protected("minimap2/{sample}/{sample}.sorted.haplotagged.bam"),
    threads: config["threads"]
    log:
        "logs/{sample}/phasenhaplotag_bam_clair3.log",
    shell:
        """
        {{ run_clair3.sh \\
        --use_whatshap_for_final_output_haplotagging \\
        --use_whatshap_for_final_output_phasing \\
        --enable_phasing \\
        --bam_fn={input.bam} \\
        --ref_fn={input.fasta} \\
        --threads={threads} \\
        --platform=ont \\
        --model_path={input.model_clair3} \\
        --output={output.dir_out}

        mv {output.dir_out}/phased_output.bam {output.bam_haplotagged}
        samtools index -@ {threads} {output.bam_haplotagged}

        echo -e "[INFO] Clair3 is done!"; }} \\
        1> {log} 2>&1
        """


rule call_sv_severus:
    input:
        bam_haplotagged="minimap2/{sample}/{sample}.sorted.haplotagged.bam",
        vcf_phased="clair3/{sample}/phased_merge_output.vcf.gz",
        fasta=config["fasta"],
        bed_nvtr=config["bed_nvtr"],
    output:
        dir_out=directory("severus/{sample}"),
        vcf=protected("severus/{sample}/severus.vcf"),
        tab_rename=temp("severus/{sample}/rename.tab"),
    params:
        min_num_reads=config["min_num_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_length_sv=config["min_length_sv"],
    threads: config["threads"]
    log:
        "logs/{sample}/call_sv_severus.log",
    shell:
        """
        {{ echo -e "{wildcards.sample}.sorted.haplotagged\\t{wildcards.sample}" > {output.tab_rename}
        severus \\
        --target-bam {input.bam_haplotagged} \\
        --out-dir {output.dir_out} \\
        -t {threads} \\
        --min-mapq {params.min_quality_mapping} \\
        --min-support {params.min_num_reads} \\
        --min-sv-size {params.min_length_sv} \\
        --between-junction-ins \\
        --output-read-ids \\
        --phasing-vcf {input.vcf_phased} \\
        --vntr-bed {input.bed_nvtr}

        awk '/^##FILTER/ && !f {{print "##FILTER=<ID=FAIL_LONG,Description=\\"FAIL_LONG\\">\\n##FILTER=<ID=FAIL,Description=\\"FAIL\\">"; f=1}} 1' {output.dir_out}/all_SVs/severus_all.vcf \\
            | awk -F'\\t' -v OFS="\\t" '$3 ~ /.*BND.*/ {{gsub(/END=[0-9]+;/, ""); print $0; next}} {{print $0}}' \\
            | bcftools reheader -s {output.tab_rename} \\
            > {output.vcf}

        echo -e "[INFO] Severus is done!"; }} \\
        1> {log} 2>&1
        """
