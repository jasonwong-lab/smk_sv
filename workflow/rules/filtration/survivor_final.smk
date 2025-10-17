rule survivor_final:
    priority: 10
    conda:
        "../../envs/survivor.yaml"
    input:
        vcfs=expand(
            "{caller}/{{sample}}/merged/filtered/{caller}.{{type_sv}}.vcf",
            caller=CALLERS,
        ),
    output:
        list="survivor/{sample}/final/vcfs.{type_sv}.list",
        vcf_tmp=temp("survivor/{sample}/final/{sample}.{type_sv}.merged.vcf.tmp"),
        tab=temp("survivor/{sample}/final/{sample}.{type_sv}.rename.tab"),
        vcf="survivor/{sample}/final/{sample}.{type_sv}.merged.vcf",
    params:
        min_size=config["min_size"],
        callers=CALLERS,
        distance_sv=lambda wildcards: config["distance_sv"][wildcards.type_sv],
        n_callers=lambda wildcards: config["n_callers"][wildcards.type_sv] - 1,
        consider_type=lambda wildcards: (
            1 if config["consider_type"][wildcards.type_sv] else 0
        ),
        consider_strand=lambda wildcards: (
            1 if config["consider_strand"][wildcards.type_sv] else 0
        ),
        estimate_distance=lambda wildcards: (
            1 if config["estimate_distance"][wildcards.type_sv] else 0
        ),
    log:
        "logs/{sample}/survivor_final.{type_sv}.log",
    script:
        "../../scripts/survivor.sh"
