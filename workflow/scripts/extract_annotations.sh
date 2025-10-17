#!/usr/bin/env bash
# shellcheck disable=SC2154

set -x

{ vep="${snakemake_input[vep]}"
snpeff="${snakemake_input[snpeff]}"
vcf="${snakemake_input[vcf]}"
survivor="${snakemake_input[survivor]}"
tab="${snakemake_output[tab]}"
ids="${snakemake_output[ids]}"
vep_extracted="${snakemake_output[vep]}"
snpeff_extracted="${snakemake_output[snpeff]}"
vcf_extracted="${snakemake_output[vcf]}"
callers="${snakemake_params[callers]}"

caller="${snakemake_wildcards[caller]}"


IFS=" " read -r -a callers_unsorted <<< "${callers}"
list_string=$(printf "%s\n" "${callers_unsorted[@]}")
mapfile -t callers < <(echo "${list_string}" | sort)

caller_column=$(python -c "callers = \"${callers[*]}\".split(); map_caller_column = {callers[i]: i+5 for i in range(len(callers))}; print(map_caller_column['${caller}'])")

bcftools query -f'%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%ID]' "${survivor}" \
    | awk -F'\t' -v OFS='\t' '{
        printf "%s\t%s\t%s\t%s", $1, $2, $3, $4
        for (i = 5; i <= NF; i++) {split($i, arr, "="); printf "\t%s", arr[2]}
        printf "\n"
    }' \
    > "${tab}"

cut -f "${caller_column}" "${tab}" | { grep -v 'NaN' || true; } > "${ids}"

bcftools view -i "ID=@${ids}" "${vep}" > "${vep_extracted}"
bcftools view -i "ID=@${ids}" "${snpeff}" > "${snpeff_extracted}"
bcftools view -i "ID=@${ids}" "${vcf}" > "${vcf_extracted}"; } \
1> "${snakemake_log[0]}" 2>&1
