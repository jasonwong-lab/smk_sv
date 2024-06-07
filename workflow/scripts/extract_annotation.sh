#!/usr/bin/env bash
# shellcheck disable=SC2154

set -x

vcf_vep="${snakemake_input[vcf_vep]}"
vcf_merged="${snakemake_input[vcf_merged]}"
tab="${snakemake_output[tab]}"
ids="${snakemake_output[ids]}"
vcf_extracted="${snakemake_output[vcf_extracted]}"

callers="${snakemake_params[callers]}"

log="${snakemake_log[0]}"

caller="${snakemake_wildcards[caller]}"


{ IFS=" " read -r -a callers_unsorted <<< "${callers}"
list_string=$(printf "%s\n" "${callers_unsorted[@]}")
mapfile -t callers < <(echo "${list_string}" | sort)

caller_column=$(python -c "callers = \"${callers[*]}\".split(); map_caller_column = {callers[i]: i+5 for i in range(len(callers))}; print(map_caller_column['${caller}'])")

bcftools query -f'%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%ID]' "${vcf_merged}" \
    | awk -F'\t' -v OFS='\t' '{
        printf "%s\t%s\t%s\t%s\t", $1, $2, $3, $4
        for (i = 5; i <= NF; i++) {split($i, arr, "="); printf "%s\t", arr[2]}
        printf "\n"
    }' \
> "${tab}"

cut -f "${caller_column}" "${tab}" | { grep -v 'NaN' || true; } > "${ids}"

bcftools view -i "ID=@${ids}" "${vcf_vep}" > "${vcf_extracted}"; } \
1> "${log}" 2>&1
