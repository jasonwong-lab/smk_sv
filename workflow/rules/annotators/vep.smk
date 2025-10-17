rule vep:
    conda:
        "../../envs/vep.yaml"
    input:
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
        fasta=config["fasta"],
        dir_cache=path_cache_vep,
    output:
        vcf=protected("{caller}/{sample}/{caller}.{type_sv}.vep.vcf"),
        html="{caller}/{sample}/{caller}.{type_sv}.vep.html",
    params:
        cache=config["cache_vep"],
        version=config["version_vep"],
        genome=config["genome"],
        species=config["species"],
        max_size=config["max_size_vep"],
        arg_buffer_size=lambda wildcards: (
            "--buffer_size 50" if wildcards.type_sv == "BND" else ""
        ),
    threads: 1
    log:
        "logs/{sample}/vep.{caller}.{type_sv}.log",
    shell:
        """
        max_size_mb=300
        size_input=$(du -sm {input.vcf} | cut -f1)

        input={input.vcf}
        if [ {wildcards.caller} == "svision" ] && [ {wildcards.type_sv} == "INS" ] ; then
            awk 'BEGIN {{FS=OFS="\\t"}} $5 == "N" {{$5 = "<INS>"}} 1' ${{input}} > ${{input%.*}}.5.vcf
            input_vep=${{input%.*}}.5.vcf
        elif [ {wildcards.caller} == "cutesv" ] && [ {wildcards.type_sv} == "DEL" ] && [ ${{size_input}} -ge ${{max_size_mb}} ]; then
            awk 'BEGIN {{OFS=FS="\\t"}} !/^#/ {{$4 = "N"; $5 = "<DEL>"}} 1' ${{input}} > ${{input_vep%.*}}.4_5.vcf
            input_vep=${{input%.*}}.4_5.vcf
        fi
        if [ {wildcards.type_sv} == "INS" ] || [ {wildcards.type_sv} == "DEL" ]; then
            sed -e 's/SVTYPE=[^;]*;//g' ${{input}} > ${{input%.*}}.SVTYPE.vcf
            input=${{input%.*}}.SVTYPE.vcf
        fi

        vep \\
            -i {input.vcf} -o {output.vcf} --stats_file {output.html} \\
            --species {params.species} --assembly {params.genome} \\
            --cache_version {params.version} --fasta {input.fasta} \\
            --dir_cache {params.cache} --fork {threads} \\
            --force_overwrite --cache --vcf --everything --filter_common \\
            --per_gene --total_length --offline --format vcf --dont_skip \\
            --max_sv_size {params.max_size} \\
            {params.arg_buffer_size} \\
            1> {log} 2>&1
        """
