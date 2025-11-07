rule snpeff:
    conda:
        "../../envs/snpeff.yaml"
    input:
        vcf="{caller}/{sample}/{caller}.{type_sv}.vcf",
        fasta=config["fasta"],
        dir_cache=path_cache_snpeff,
    output:
        vcf=protected("{caller}/{sample}/{caller}.{type_sv}.snpeff.vcf"),
        html="{caller}/{sample}/{caller}.{type_sv}.snpeff.html",
    params:
        cache=config["cache_snpeff"],
        version=config["version_snpeff"],
        genome=config["genome"],
    resources:
        mem_mb=1,
    log:
        "logs/{sample}/snpeff.{caller}.{type_sv}.log",
    shell:
        """
        {{ max_size_mb=300
        size_input=$(du -sm {input.vcf} | cut -f1)

        input={input.vcf}
        if [ {wildcards.caller} == "svision" ] && [ {wildcards.type_sv} == "INS" ]; then
            awk 'BEGIN {{OFS=FS="\\t"}} !/^#/ {{$5 = "INS"}} 1' ${{input}} > ${{input%.*}}.5.vcf
            input=${{input%.*}}.5.vcf
        elif [ {wildcards.caller} == "svision" ] && [ {wildcards.type_sv} != "INS" ]; then
            awk -v var={wildcards.type_sv} 'BEGIN {{OFS=FS="\\t"}} !/^#/ {{$5 = "<"var">"}} 1' ${{input}} > ${{input%.*}}.5.vcf
            input=${{input%.*}}.5.vcf
        elif [ {wildcards.caller} == "debreak" ] && [ {wildcards.type_sv} == "DEL" ] && [ ${{size_input}} -ge ${{max_size_mb}} ]; then
            awk 'BEGIN {{OFS=FS="\\t"}} !/^#/ {{$4 = "N"; $5 = "<DEL>"}} 1' ${{input}} > ${{input%.*}}.4_5.vcf
            input=${{input%.*}}.4_5.vcf
        fi

        snpEff -Xmx{resources.mem_mb}M \\
            -nodownload -v -lof -canon \\
            -dataDir {params.cache} -s {output.html} \\
            {params.genome}.{params.version} ${{input}} \\
            1> {output.vcf}; }} \\
        1> {log} 2>&1
        """
