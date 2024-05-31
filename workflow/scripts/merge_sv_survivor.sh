#!/usr/bin/env bash
# shellcheck disable=SC2154

set -x

vcfs="${snakemake_input[vcfs]}"
list_vcfs="${snakemake_output[list_vcfs]}"
vcf_merged_tmp="${snakemake_output[vcf_merged_tmp]}"
tab_rename="${snakemake_output[tab_rename]}"
vcf_merged="${snakemake_output[vcf_merged]}"

min_length_sv="${snakemake_params[min_length_sv]}"
caller="${snakemake_params[caller]}"
merge_distance_sv="${snakemake_params[merge_distance_sv]}"
merge_nbr_callers="${snakemake_params[merge_nbr_callers]}"
merge_type_sv="${snakemake_params[merge_type_sv]}"
merge_strand="${snakemake_params[merge_strand]}"
merge_estimate_distance="${snakemake_params[merge_estimate_distance]}"

log="${snakemake_log[0]}"

sample="${snakemake_wildcards[sample]}"
type_sv="${snakemake_wildcards[type_sv]}"


{ IFS=" " read -r -a caller_unsorted <<< "${caller}"
list_string=$(printf "%s\n" "${caller_unsorted[@]}")
mapfile -t caller < <(echo "${list_string}" | sort)
# for clr in "${caller[@]}"; do
#     realpath "${clr}"/"${sample}"/"${clr}"."${type_sv}".vcf
# done > "${list_vcfs}"

IFS=" " read -r -a vcfs_unsorted <<< "${vcfs}"
list_string=$(printf "%s\n" "${vcfs_unsorted[@]}")
mapfile -t vcfs < <(echo "${list_string}" | sort)
printf "%s\n" "${vcfs[@]}" | xargs realpath > "${list_vcfs}"


formula=("${merge_distance_sv}" "${merge_nbr_callers}" "${merge_type_sv}" "${merge_strand}" "${merge_estimate_distance}" "${min_length_sv}")

SURVIVOR merge "${list_vcfs}" "${formula[@]}" "${vcf_merged_tmp}"

for i in "${!caller[@]}"; do
    line="${sample}_${i}\t${sample}_${caller[${i}]}\n"
    output+=$line
done
echo -e "${output}" | sed 's/^ *//' | sed 's/_0//' > "${tab_rename}"
bcftools reheader -s "${tab_rename}" "${vcf_merged_tmp}" > "${vcf_merged}"; } \
1> "${log}" 2>&1
