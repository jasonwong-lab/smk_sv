rule merge_sv_survivor:
    input:
        vcfs=expand(
            "{caller}/{{sample}}/{caller}.{{type_sv}}.vcf",
            caller=list(config["callers"].keys()),
        ),
    output:
        list_vcfs="survivor/{sample}/vcfs.{type_sv}.list",
        vcf_merged_tmp=temp("survivor/{sample}/{sample}.{type_sv}.merged.vcf.tmp"),
        tab_rename=temp("survivor/{sample}/{sample}.{type_sv}.rename.tab"),
        vcf_merged="survivor/{sample}/{sample}.{type_sv}.merged.vcf",
    params:
        min_length_sv=config["min_length_sv"],
        callers=list(config["callers"].keys()),
        merge_distance_sv=lambda wildcards: config["merge_distance_sv"][wildcards.type_sv],
        merge_nbr_callers=lambda wildcards: config["merge_nbr_callers"][wildcards.type_sv],
        merge_type_sv=lambda wildcards: 1 if config["merge_type_sv"][wildcards.type_sv] else 0,
        merge_strand=lambda wildcards: 1 if config["merge_strand"][wildcards.type_sv] else 0,
        merge_estimate_distance=lambda wildcards: 1 if config["merge_estimate_distance"][wildcards.type_sv] else 0,
    log:
        "logs/{sample}/merge_sv_survivor.{type_sv}.log",
    script:
        "../scripts/merge_sv_survivor.sh"
