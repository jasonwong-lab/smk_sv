#! /bin/bash

# """
# {{ caller_unsorted=({params.caller})
# list_string=$(printf "%s\\n" "${{caller_unsorted[@]}}")
# caller=($(echo "${{list_string}}" | sort))

# for clr in ${{caller[@]}}; do
#     echo $(realpath ${{clr}}/{wildcards.sample}/${{clr}}.{wildcards.type_sv}.vcf)
# done > {output.list_vcfs}

# # declare -A map_formula
# # map_formula=(["BND"]="10 3 0 0 1 1" ["DEL"]="10 3 0 0 1 {params.min_length_sv}" ["INS"]="10 3 0 0 1 {params.min_length_sv}" ["INV"]="10 3 0 0 1 {params.min_length_sv}" ["DUP"]="10 3 0 0 1 {params.min_length_sv}")
# # formula=${{map_formula["{wildcards.type_sv}"]}}

# formula="{params.merge_distance_sv} {params.merge_nbr_callers} {params.merge_type_sv} {params.merge_strand} {params.merge_estimate_distance} {params.min_length_sv}"

# SURVIVOR merge {output.list_vcfs} ${{formula}} {output.vcf_merged_tmp}

# for i in "${{!caller[@]}}"; do
#     line="{wildcards.sample}_${{i}}\\t{wildcards.sample}_${{caller[$i]}}\\n"
#     output+=$line
# done
# echo -e ${{output}} | sed 's/^ *//' | sed 's/_0//' > {output.tab_rename}
# bcftools reheader -s {output.tab_rename} {output.vcf_merged_tmp} > {output.vcf_merged}; }} \\
# 1> {log} 2>&1
# """
