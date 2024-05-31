rule convert_annotation:
    container:
        None
    input:
        vcf_extracted="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vep.vcf",
        fasta=config["fasta"],
        dir_db_vep=config["dir_db_vep"],
    output:
        maf="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vep.maf",
        tsv="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vep.tsv",
    params:
        genome=config["genome"],
        version_cache_vep=config["version_cache_vep"],
    log:
        "logs/{sample}/convert_annotation.{caller}.{type_sv}.log",
    script:
        "../scripts/convert_annotation.sh"
