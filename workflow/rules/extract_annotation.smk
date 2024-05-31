# TODO: Need modification
rule extract_annotation:
    input:
        vcf_vep="{caller}/{sample}/{caller}.{type_sv}.snpeff.vep.vcf",
        vcf_merged="survivor/{sample}/{sample}.{type_sv}.merged.vcf",
    output:
        tab="survivor/{sample}/{sample}.{caller}.{type_sv}.tab",
        ids=touch(temp("survivor/{sample}/{caller}.{type_sv}.id")),
        vcf_extracted="{caller}/{sample}/merged/{caller}.{type_sv}.snpeff.vep.vcf",
    params:
        caller=config["caller"],
    log:
        "logs/{sample}/extract_annotation.{caller}.{type_sv}.log",
    shell:
        """
        {{
        # caller_unsorted=({params.caller})
        # list_string=$(printf "%s\\n" "${{caller_unsorted[@]}}")
        # caller=($(echo "${{list_string}}" | sort))

        declare -A map_caller_column
        map_caller_column=(["cutesv"]=5 ["severus"]=6 ["sniffles"]=7 ["svim"]=8 ["svision"]=9)
        caller_column=${{map_caller_column["{wildcards.caller}"]}}

        bcftools query -f'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%ID]' {input.vcf_merged} \\
            | awk -F'\\t' -v OFS='\\t' '{{
                split($5, cutesv, "=");
                split($6, severus, "=");
                split($7, sniffles, "=");
                split($8, svim, "=");
                split($9, svision, "=");
                print($1, $2, $3, $4, cutesv[2], severus[2], sniffles[2], svim[2], svision[2])
            }}' \\
        > {output.tab}

        cut -f ${{caller_column}} {output.tab} | {{ grep -v 'NaN' || true; }} > {output.ids}

        bcftools view -i "ID=@{output.ids}" {input.vcf_vep} > {output.vcf_extracted}; }} \\
        1> {log} 2>&1
        """