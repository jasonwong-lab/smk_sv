rule phasenhaplotag_bam_clair3:
    input:
        bam=ancient("minimap2/{sample}/{sample}.sorted.bam"),
        fasta=config["fasta"],
        model_clair3=config["model_clair3"],
    output:
        vcf_phased=protected("clair3/{sample}/phased_merge_output.vcf.gz"),
        bam_haplotagged=protected("minimap2/{sample}/{sample}.sorted.haplotagged.bam"),
        bai_haplotagged=protected(
            "minimap2/{sample}/{sample}.sorted.haplotagged.bam.bai"
        ),
    params:
        dir_out=directory("clair3/{sample}"),
    threads: math.floor(workflow.cores / NUMBER_SAMPLES) # config["threads"]
    log:
        "logs/{sample}/phasenhaplotag_bam_clair3.log",
    shell:
        """
        {{ [ ! -f {input.fasta}.fai ] && samtools faidx {input.fasta}
        source activate clair3
        run_clair3.sh \\
        --use_whatshap_for_final_output_haplotagging \\
        --use_whatshap_for_final_output_phasing \\
        --enable_phasing \\
        --bam_fn={input.bam} \\
        --ref_fn={input.fasta} \\
        --threads={threads} \\
        --platform=ont \\
        --model_path={input.model_clair3} \\
        --output={params.dir_out}

        mv {params.dir_out}/phased_output.bam {output.bam_haplotagged}
        mv {params.dir_out}/phased_output.bam.bai {output.bai_haplotagged}
        # samtools index -@ {threads} {output.bam_haplotagged}

        echo -e "[INFO] Clair3 is done!"; }} \\
        1> {log} 2>&1
        """


rule call_sv_severus:
    input:
        bam_haplotagged=ancient("minimap2/{sample}/{sample}.sorted.haplotagged.bam"),
        vcf_phased=ancient("clair3/{sample}/phased_merge_output.vcf.gz"),
        fasta=config["fasta"],
        bed_nvtr=config["bed_nvtr"],
    output:
        vcf=protected("severus/{sample}/severus.vcf"),
        tab_rename=temp("severus/{sample}/rename.tab"),
    params:
        dir_out=directory("severus/{sample}"),
        min_num_reads=config["min_num_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_length_sv=config["min_length_sv"],
    threads: math.floor(workflow.cores / NUMBER_SAMPLES / NUMBER_CALLERS)
    log:
        "logs/{sample}/call_sv_severus.log",
    shell:
        """
        {{ echo -e "{wildcards.sample}.sorted.haplotagged\\t{wildcards.sample}" > {output.tab_rename}
        severus \\
        --target-bam {input.bam_haplotagged} \\
        --out-dir {params.dir_out} \\
        -t {threads} \\
        --min-mapq {params.min_quality_mapping} \\
        --min-support {params.min_num_reads} \\
        --min-sv-size {params.min_length_sv} \\
        --between-junction-ins \\
        --output-read-ids \\
        --phasing-vcf {input.vcf_phased} \\
        --vntr-bed {input.bed_nvtr}

        awk '/^##FILTER/ && !f {{print "##FILTER=<ID=FAIL_LONG,Description=\\"FAIL_LONG\\">\\n##FILTER=<ID=FAIL,Description=\\"FAIL\\">"; f=1}} 1' {params.dir_out}/all_SVs/severus_all.vcf \\
            | awk -F'\\t' -v OFS="\\t" '$3 ~ /.*BND.*/ {{gsub(/END=[0-9]+;/, ""); print $0; next}} {{print $0}}' \\
            | bcftools reheader -s {output.tab_rename} \\
            > {output.vcf}

        echo -e "[INFO] Severus is done!"; }} \\
        1> {log} 2>&1
        """
