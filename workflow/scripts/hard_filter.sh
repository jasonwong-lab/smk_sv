#!/usr/bin/env bash
# shellcheck disable=SC2154

set -x

{ vcf=${snakemake_input[vcf]}
vcf_filtered=${snakemake_output[vcf]}
min_size=${snakemake_params[min_size]}
min_reads=${snakemake_params[min_reads]}
min_coverage=${snakemake_params[min_coverage]}
min_dhffc=${snakemake_params[min_dhffc]}
max_dhbfc=${snakemake_params[max_dhbfc]}
caller=${snakemake_wildcards[caller]}


formula_length="INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < ${min_size}"
formula_del="INFO/SVTYPE == \"DEL\" & FMT/DHFFC[0] > ${min_dhffc}"
formula_dup="INFO/SVTYPE == \"DUP\" & FMT/DHBFC[0] < ${max_dhbfc}"

if [ "${caller}" == "cutesv" ]; then

    formula_reads="INFO/RE < ${min_reads} | FMT/DV[0] < ${min_reads}"
    formula_coverage="FMT/DR[0] < ${min_reads} & (FMT/DV[0] + FMT/DR[0]) < ${min_coverage}"

elif [ "${caller}" == "svim" ]; then

    formula_reads="INFO/SUPPORT < ${min_reads} | (SVTYPE = \"DEL,INS\" & FMT/AD[0:1] < ${min_reads})"
    formula_coverage="SVTYPE = \"DEL,INS\" & FMT/AD[0:0] < ${min_reads} & FMT/DP[0] < ${min_coverage}"

elif [ "${caller}" == "severus" ]; then

    formula_reads="FMT/DV[0] < ${min_reads} "
    formula_coverage="FMT/DR[0] < ${min_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}"

elif [ "${caller}" == "debreak" ]; then

    formula_reads="INFO/SUPPREAD < ${min_reads}"
    formula_coverage=""

elif [ "${caller}" == "nanovar" ]; then

    formula_reads="INFO/SR < ${min_reads} | FMT/AD[0:1] < ${min_reads} "
    formula_coverage="FMT/AD[0:0] < ${min_reads} & FMT/DP[0] < ${min_coverage}"

elif [ "${caller}" == "delly" ]; then

    formula_reads="INFO/SR < ${min_reads} | FMT/RV[0] < ${min_reads}"
    formula_coverage="FMT/RR[0] < ${min_reads} & (FMT/RV[0] + FMT/RR[0]) < ${min_coverage}"

elif [ "${caller}" == "nanosv" ]; then

    formula_reads="FMT/DV[0:0] < ${min_reads}"
    formula_coverage="SUM(FMT/DR[0:]) < ${min_reads} & (FMT/DV[0:0] + SUM(FMT/DR[0:])) < ${min_coverage}"

elif [ "${caller}" == "svision" ] || [ "${caller}" == "sniffles" ]; then

    formula_reads="INFO/SUPPORT < ${min_reads} | FMT/DV[0] < ${min_reads}"
    formula_coverage="FMT/DR[0] < ${min_reads} & (FMT/DV[0] + FMT/DR[0]) < ${min_coverage}"

else
    echo "$(date +"%Y-%m-%d %H:%M:%S") [ERROR] Unknown caller: ${caller}"
    exit 1
fi

formula="(${formula_length}) | (${formula_coverage}) | (${formula_reads}) | (${formula_del}) | (${formula_dup})"

grep -E "^#|^chr" "${vcf}" \
    | bcftools filter -e "${formula}" -Ov - \
    > "${vcf_filtered}"; } \
1> "${snakemake_log[0]}" 2>&1
