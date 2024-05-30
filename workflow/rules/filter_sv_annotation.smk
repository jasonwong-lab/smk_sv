# TODO: need modification
rule filter_sv_annotation:
    input:
        annotsv="{caller}/{sample}/{caller}.{type_sv}.annotsv.tsv",
        maf="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vep.maf",
        tsv="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vep.tsv",
        tab="survivor/{sample}/{sample}.{caller}.{type_sv}.tab",
        vcf="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vep.vcf",
    output:
        vcf_final="{caller}/{sample}/merged/filtered/{caller}.{type_sv}.snpeff.vep.vcf",
        ids=temp("{caller}/{sample}/merged/filtered/tmp.{type_sv}"),
        rdata="{caller}/{sample}/merged/filtered/{caller}.{type_sv}.RData",
        table="{caller}/{sample}/merged/filtered/{caller}.{type_sv}.tsv",
    params:
        libs_r=config["libs_r"],
        caller=config["caller"],
    script:
        "scripts/filter_sv_annotation.R"
