rule annotate_sv_snpeffnvep:
    conda:
        "../../envs/annotators.yaml"
    input:
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
        fasta=config["fasta"],
        cache_snpeff=config["cache_snpeff"],
        cache_vep=config["cache_vep"],
    output:
        vcf_snpeff=protected("{caller}/{sample}/{caller}.{type_sv}.snpeff.vcf"),
        vcf_vep=protected("{caller}/{sample}/{caller}.{type_sv}.snpeff.vep.vcf"),
    params:
        genome=config["genome"],
        cache_snpeff=config["cache_snpeff"],
        version_vep=config["version_vep"],
    threads: 1
    log:
        "logs/{sample}/annotate_sv_snpeffnvep.{caller}.{type_sv}.log",
    script:
        "../../scripts/annotate_sv.sh"


rule annotate_sv_annotsv:
    input:
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
    output:
        tsv=touch(protected("{caller}/{sample}/{caller}.{type_sv}.annotsv.tsv")),
    params:
        genome=config["genome"],
    resources:
        n_instance=1,
    log:
        "logs/{sample}/annotate_sv_annotsv.{caller}.{type_sv}.log",
    script:
        "../../scripts/annotate_sv_annotsv.sh"
