rule annotsv:
    shadow:
        "minimal"
    conda:
        "../../envs/annotsv.yaml"
    input:
        **get_annotsv_cache_outputs(),
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
    output:
        tsv=touch(protected("{caller}/{sample}/{caller}.{type_sv}.annotsv.tsv")),
    params:
        genome=config["genome"],
        dir_cache=config["cache_annotsv"],
    log:
        "logs/{sample}/annotsv.{caller}.{type_sv}.log",
    shell:
        """
        {{ g={params.genome}
        if [[ "${{g}}" == "hg19" ]]; then
            g="GRCh37"
        elif [[ "${{g}}" == "hg38" ]]; then
            g="GRCh38"
        fi

        max_size_mb=300
        size_input=$(du -sm {input.vcf} | cut -f1)

        input={input.vcf}
        if {{ [ {wildcards.caller} == "cutesv" ] && [ {wildcards.type_sv} == "DEL" ]; }} \\
            || {{ [ {wildcards.caller} == "debreak" ] && [ {wildcards.type_sv} == "DEL" ]; }} \\
            && [ ${{size_input}} -ge ${{size_max}} ]; then
            awk 'BEGIN {{OFS=FS="\\t"}} !/^#/ {{$4 = "N"; $5 = "<DEL>"}} 1' ${{input}} > ${{input%.*}}.4_5.vcf
            input=${{input%.*}}.4_5.vcf
        fi

        AnnotSV \\
            -genomeBuild ${{g}} \\
            -annotationsDir {params.dir_cache} \\
            -SVinputFile ${{input}} \\
            -outputFile {output.tsv} \\
            -SVminSize 1 \\
            -overwrite 1; }} \\
        1> {log} 2>&1
        """
