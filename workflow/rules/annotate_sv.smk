rule annotate_sv_snpeffnvep:
    container:
        None
    input:
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
        fasta=config["fasta"],
        dir_db_snpeff=config["dir_db_snpeff"],
        dir_db_vep=config["dir_db_vep"],
    output:
        vcf_snpeff=protected("{caller}/{sample}/{caller}.{type_sv}.snpeff.vcf"),
        vcf_vep=protected("{caller}/{sample}/{caller}.{type_sv}.snpeff.vep.vcf"),
    params:
        genome=config["genome"],
        version_cache_snpeff=config["version_cache_snpeff"],
        version_cache_vep=config["version_cache_vep"],
    threads: config["threads"]
    log:
        "logs/{sample}/annotate_sv_snpeffnvep.{caller}.{type_sv}.log",
    script:
        "../scripts/annotate_sv.sh"


rule annotate_sv_annotsv:
    container:
        None
    input:
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
    output:
        tsv=touch(protected("{caller}/{sample}/{caller}.{type_sv}.annotsv.tsv")),
    params:
        genome=config["genome"],
    threads: workflow.cores
    log:
        "logs/{sample}/annotate_sv_annotsv.{caller}.{type_sv}.log",
    script:
        "../scripts/annotate_sv_annotsv.sh"
