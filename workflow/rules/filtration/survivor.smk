rule survivor:
    conda:
        "../../envs/survivor.yaml"
    input:
        vcfs=expand(
            "{caller}/{{sample}}/{caller}.{{type_sv}}.vcf",
            caller=CALLERS,
        ),
    output:
        list="survivor/{sample}/vcfs.{type_sv}.list",
        vcf_tmp=temp("survivor/{sample}/{sample}.{type_sv}.merged.vcf.tmp"),
        tab=temp("survivor/{sample}/{sample}.{type_sv}.rename.tab"),
        vcf="survivor/{sample}/{sample}.{type_sv}.merged.vcf",
    params:
        min_size=config["min_size"],
        callers=CALLERS,
        distance_sv=lambda wildcards: config["distance_sv"][wildcards.type_sv],
        n_callers=lambda wildcards: config["n_callers"][wildcards.type_sv],
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
        "logs/{sample}/survivor.{type_sv}.log",
    script:
        "../../scripts/survivor.sh"
