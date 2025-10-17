rule filter_annotations:
    conda:
        "../../envs/r.yaml"
    input:
        annotsv="{caller}/{sample}/{caller}.{type_sv}.annotsv.tsv",
        maf="{caller}/{sample}/merged/{caller}.{type_sv}.vep.maf",
        tsv="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.tsv",
        tab="survivor/{sample}/{sample}.{caller}.{type_sv}.tab",
        vep="{caller}/{sample}/merged/{caller}.{type_sv}.vep.vcf",
        snpeff="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vcf",
        vcf="{caller}/{sample}/merged/{caller}.{type_sv}.vcf",
    output:
        vep="{caller}/{sample}/merged/filtered/{caller}.{type_sv}.vep.vcf",
        snpeff="{caller}/{sample}/merged/filtered/{caller}.{type_sv}.snpeff.vcf",
        vcf="{caller}/{sample}/merged/filtered/{caller}.{type_sv}.vcf",
        ids=temp("{caller}/{sample}/merged/filtered/tmp.{type_sv}"),
        rdata="{caller}/{sample}/merged/filtered/{caller}.{type_sv}.RData",
        table="{caller}/{sample}/merged/filtered/{caller}.{type_sv}.tsv",
    params:
        callers=CALLERS,
        terms_relative=config["terms_relative"],
    threads: 1
    script:
        "../../scripts/filter_annotations.R"
