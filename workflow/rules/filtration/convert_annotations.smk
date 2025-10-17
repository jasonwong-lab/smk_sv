rule convert_vep:
    conda:
        "../../envs/vcf2maf.yaml"
    input:
        vcf="{caller}/{sample}/merged/{caller}.{type_sv}.vep.vcf",
        fasta=config["fasta"],
    output:
        maf="{caller}/{sample}/merged/{caller}.{type_sv}.vep.maf",
    params:
        version=config["version_vep"],
        genome=config["genome"],
        cache=config["cache_vep"],
        species=config["species"],
    log:
        "logs/{sample}/convert_vep.{caller}.{type_sv}.log",
    shell:
        """
        {{ input={input.vcf}
        if [ {wildcards.caller} == "svim" ] && [ {wildcards.type_sv} == "DUP" ]; then
            sed -e 's/DUP:TANDEM/DUP/g' -e 's/DUP:INT/DUP/g' ${{input}} > ${{input%.*}}.DUP.vcf
            input=${{input%.*}}.DUP.vcf
        elif [ {wildcards.caller} == "svision" ]; then
            perl -pe "s/SVTYPE=.*?;/SVTYPE={wildcards.type_sv};/g" ${{input}} > ${{input%.*}}.{wildcards.type_sv}.vcf
            input=${{input%.*}}.{wildcards.type_sv}.vcf
        fi

        vcf2maf.pl \\
            --input-vcf ${{input}} --output-maf {output.maf} \\
            --ncbi-build {params.genome} --cache-version {params.version} \\
            --ref-fasta {input.fasta} --vcf-tumor-id {wildcards.sample} --tumor-id {wildcards.sample} --vep-data {params.cache} \\
            --species {params.species} --inhibit-vep; }} \\
        1> {log} 2>&1
        """


rule convert_snpeff:
    conda:
        "../../envs/snpeff.yaml"
    input:
        vcf="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vcf",
    output:
        tsv="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.tsv",
    params:
        fields_common=FIELDS_COMMON,
        fields_fmt=get_convert_snpeff_arguments,
    log:
        "logs/{sample}/convert_snpeff.{caller}.{type_sv}.log",
    shell:
        """
        {{ SnpSift extractFields \\
            -s ";" -e "." \\
            {input.vcf} {params.fields_common} {params.fields_fmt} \\
                | sed '1s/GEN\\[\\*\\]\\.//g ; 1s/ANN\\[\\*\\]\\.//g ; 1s/\\[\\*\\]//g' \\
                    > {output.tsv}; }} \\
        > {log} 2>&1
        """
