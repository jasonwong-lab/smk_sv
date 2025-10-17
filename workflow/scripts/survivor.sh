#!/usr/bin/env bash
# shellcheck disable=SC2154

set -x

{ vcfs="${snakemake_input[vcfs]}"
list="${snakemake_output[list]}"
vcf_tmp="${snakemake_output[vcf_tmp]}"
tab="${snakemake_output[tab]}"
vcf="${snakemake_output[vcf]}"
min_size="${snakemake_params[min_size]}"
callers="${snakemake_params[callers]}"
distance_sv="${snakemake_params[distance_sv]}"
n_callers="${snakemake_params[n_callers]}"
consider_type="${snakemake_params[consider_type]}"
consider_strand="${snakemake_params[consider_strand]}"
estimate_distance="${snakemake_params[estimate_distance]}"
sample="${snakemake_wildcards[sample]}"
type_sv="${snakemake_wildcards[type_sv]}"


IFS=" " read -r -a callers_unsorted <<< "${callers}"
list_string=$(printf "%s\n" "${callers_unsorted[@]}")
mapfile -t callers < <(echo "${list_string}" | sort)

IFS=" " read -r -a vcfs_unsorted <<< "${vcfs}"
list_string=$(printf "%s\n" "${vcfs_unsorted[@]}")
mapfile -t vcfs < <(echo "${list_string}" | sort)
printf "%s\n" "${vcfs[@]}" | xargs realpath > "${list}"


formula=("${distance_sv}" "${n_callers}" "${consider_type}" "${consider_strand}" "${estimate_distance}" "${min_size}")

SURVIVOR merge "${list}" "${formula[@]}" "${vcf_tmp}"

for i in "${!callers[@]}"; do
    line="${sample}_${i}\t${sample}_${callers[${i}]}\n"
    output+=$line
done
echo -e "${output}" | sed 's/^ *//' | sed 's/_0//' > "${tab}"
bcftools reheader -s "${tab}" "${vcf_tmp}" > "${vcf}"; } \
1> "${snakemake_log[0]}" 2>&1
