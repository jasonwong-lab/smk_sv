rule clair3:
    container:
        "docker://hkubal/clair3:latest"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
        fai=f"{config['fasta']}.fai",
        model_clair3=config["model_clair3"],
    output:
        vcf=protected("clair3/{sample}/phased_merge_output.vcf.gz"),
        bam=protected("clair3/{sample}/phased_output.bam"),
        bai=protected("clair3/{sample}/phased_output.bam.bai"),
        bam_renamed="clair3/{sample}/{sample}.sorted.haplotagged.bam",
        bai_renamed="clair3/{sample}/{sample}.sorted.haplotagged.bam.bai",
    params:
        dir=directory("clair3/{sample}"),
    threads: 1
    log:
        "logs/{sample}/clair3.log",
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
            --output={params.dir}

        ln -r -s {output.bam} {output.bam_renamed}
        ln -r -s {output.bai} {output.bai_renamed}; }} \\
        1> {log} 2>&1
        """
