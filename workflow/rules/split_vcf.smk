rule split_vcf:
    input:
        vcf_filtered="{caller}/{sample}/{caller}.duphold.filtered.vcf",
    output:
        vcf_splitted="{caller}/{sample}/{caller}.{type_sv}.vcf",
    log:
        "logs/{sample}/split_vcf.{caller}.{type_sv}.log",
    script:
        "../scripts/split_vcf.sh"
