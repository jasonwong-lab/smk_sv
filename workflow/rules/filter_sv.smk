rule filter_sv:
    container:
        "docker://mhjiang97/smk_sv:latest"
    input:
        vcf="{caller}/{sample}/{caller}.vcf",
        bam="minimap2/{sample}/{sample}.sorted.bam",
        fasta=config["fasta"],
    output:
        vcf_filtered=temp("{caller}/{sample}/{caller}.filtered.vcf"),
        vcf_duphold=temp("{caller}/{sample}/{caller}.filtered.duphold.vcf"),
        vcf_duphold_filtered="{caller}/{sample}/{caller}.duphold.filtered.vcf",
    params:
        min_num_reads=config["min_num_reads"],
        min_length_sv=config["min_length_sv"],
        min_coverage=config["min_coverage"],
        min_dhffc=config["min_dhffc"],
        max_dhbfc=config["max_dhbfc"],
    threads: config["threads"]
    log:
        "logs/{sample}/filter_sv.{caller}.log",
    script:
        "../scripts/filter_sv.sh"
