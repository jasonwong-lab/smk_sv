rule hard_filter:
    conda:
        "../../envs/bcftools.yaml"
    input:
        vcf="{caller}/{sample}/{caller}.duphold.vcf",
    output:
        vcf="{caller}/{sample}/{caller}.duphold.filtered.vcf",
    params:
        min_size=config["min_size"],
        min_reads=config["min_reads"],
        min_coverage=config["min_coverage"],
        min_dhffc=config["min_dhffc"],
        max_dhbfc=config["max_dhbfc"],
    log:
        "logs/{sample}/hard_filter.{caller}.log",
    script:
        "../../scripts/hard_filter.sh"
