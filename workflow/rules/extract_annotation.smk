rule extract_annotation:
    conda:
        "../envs/bcftools.yaml"
    input:
        vcf_vep="{caller}/{sample}/{caller}.{type_sv}.snpeff.vep.vcf",
        vcf_merged="survivor/{sample}/{sample}.{type_sv}.merged.vcf",
    output:
        tab="survivor/{sample}/{sample}.{caller}.{type_sv}.tab",
        ids=touch(temp("survivor/{sample}/{caller}.{type_sv}.id")),
        vcf_extracted="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vep.vcf",
    params:
        callers=list(config["callers"].keys()),
    log:
        "logs/{sample}/extract_annotation.{caller}.{type_sv}.log",
    script:
        "../scripts/extract_annotation.sh"
