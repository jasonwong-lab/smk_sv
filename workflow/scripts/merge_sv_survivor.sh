#!/usr/bin/env bash
# shellcheck disable=SC2154

set -x

vcfs="${snakemake_input[vcfs]}"
list_vcfs="${snakemake_output[list_vcfs]}"
vcf_merged_tmp="${snakemake_output[vcf_merged_tmp]}"
tab_rename="${snakemake_output[tab_rename]}"
vcf_merged="${snakemake_output[vcf_merged]}"

min_length_sv="${snakemake_params[min_length_sv]}"
callers="${snakemake_params[callers]}"
merge_distance_sv="${snakemake_params[merge_distance_sv]}"
merge_nbr_callers="${snakemake_params[merge_nbr_callers]}"
merge_type_sv="${snakemake_params[merge_type_sv]}"
merge_strand="${snakemake_params[merge_strand]}"
merge_estimate_distance="${snakemake_params[merge_estimate_distance]}"

log="${snakemake_log[0]}"

sample="${snakemake_wildcards[sample]}"
type_sv="${snakemake_wildcards[type_sv]}"


{ IFS=" " read -r -a callers_unsorted <<< "${callers}"
list_string=$(printf "%s\n" "${callers_unsorted[@]}")
mapfile -t callers < <(echo "${list_string}" | sort)

IFS=" " read -r -a vcfs_unsorted <<< "${vcfs}"
list_string=$(printf "%s\n" "${vcfs_unsorted[@]}")
mapfile -t vcfs < <(echo "${list_string}" | sort)
printf "%s\n" "${vcfs[@]}" | xargs realpath > "${list_vcfs}"


formula=("${merge_distance_sv}" "${merge_nbr_callers}" "${merge_type_sv}" "${merge_strand}" "${merge_estimate_distance}" "${min_length_sv}")

SURVIVOR merge "${list_vcfs}" "${formula[@]}" "${vcf_merged_tmp}"

for i in "${!callers[@]}"; do
    line="${sample}_${i}\t${sample}_${callers[${i}]}\n"
    output+=$line
done
echo -e "${output}" | sed 's/^ *//' | sed 's/_0//' > "${tab_rename}"
bcftools reheader -s "${tab_rename}" "${vcf_merged_tmp}" > "${vcf_merged}"; } \
1> "${log}" 2>&1
