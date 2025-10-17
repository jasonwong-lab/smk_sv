rule separate_types:
    conda:
        "../../envs/bcftools.yaml"
    input:
        vcf="{caller}/{sample}/{caller}.duphold.filtered.vcf",
    output:
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
    params:
        max_size=config["max_size"],
    log:
        "logs/{sample}/separate_types.{caller}.{type_sv}.log",
    script:
        "../../scripts/separate_types.sh"
