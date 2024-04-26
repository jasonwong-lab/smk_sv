#!/usr/bin/env bash
# shellcheck disable=SC2154
# leave svs with genotype of 0/0

vcf="${snakemake_input[vcf]}"
vcf_filtered="${snakemake_output[vcf_filtered]}"
min_num_reads="${snakemake_params[min_num_reads]}"
min_coverage="${snakemake_params[min_coverage]}"
min_length_sv="${snakemake_params[min_length_sv]}"

log="${snakemake_log[0]}"
caller="${snakemake_wildcards[caller]}"
type_sv="${snakemake_wildcards[type_sv]}"


{ if [ "${caller}" == "cutesv" ]; then
    case "${type_sv}" in
        BND)
            formula=("INFO/IMPRECISE = 1 | FILTER = \"q5\" | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\"")
        ;;
        DEL)
            formula=("ABS(INFO/SVLEN) < ${min_length_sv} | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\" | FMT/DHFFC[0] > 0.7 | INFO/IMPRECISE = 1 | FILTER = \"q5\"")
        ;;
        INS)
            formula=("ABS(INFO/SVLEN) < ${min_length_sv} | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\" | INFO/IMPRECISE = 1 | FILTER = \"q5\"")
        ;;
        INV)
            formula=("ABS(INFO/SVLEN) < ${min_length_sv} | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\" | INFO/IMPRECISE = 1 | FILTER = \"q5\"")
        ;;
        DUP)
            formula=("ABS(INFO/SVLEN) < ${min_length_sv} | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\" | FMT/DHBFC[0] < 1.3 | INFO/IMPRECISE = 1 | FILTER = \"q5\"")
        ;;
    esac
elif [ "${caller}" == "sniffles" ]; then
    case "${type_sv}" in
        BND)
            formula=("INFO/SUPPORT < ${min_num_reads} | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\"")
        ;;
        DEL)
            formula=("INFO/SUPPORT < ${min_num_reads} | ABS(INFO/SVLEN) < ${min_length_sv} | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\" | FMT/DHFFC[0] > 0.7")
        ;;
        INS)
            formula=("INFO/SUPPORT < ${min_num_reads} | ABS(INFO/SVLEN) < ${min_length_sv} | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\"")
        ;;
        INV)
            formula=("INFO/SUPPORT < ${min_num_reads} | ABS(INFO/SVLEN) < ${min_length_sv} | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\"")
        ;;
        DUP)
            formula=("INFO/SUPPORT < ${min_num_reads} | ABS(INFO/SVLEN) < ${min_length_sv} | FMT/DV[0] < ${min_num_reads} | (FMT/DR[0] < ${min_num_reads} & FMT/DV[0] + FMT/DR[0] < ${min_coverage}) | FMT/GT[0] = \"./.\" | FMT/DHBFC[0] < 1.3")
        ;;
    esac
elif [ "${caller}" == "svim" ]; then
    case "${type_sv}" in
        BND)
            formula=("INFO/SUPPORT < ${min_num_reads} | QUAL < ${min_num_reads} | FILTER = \"hom_ref\"")
        ;;
        DEL)
            formula=("INFO/SUPPORT < ${min_num_reads} | ABS(INFO/SVLEN) < ${min_length_sv} | QUAL < ${min_num_reads} | FILTER = \"hom_ref\" | FMT/DHFFC[0] > 0.7 | FMT/GT[0] = \"./.\"")
        ;;
        INS)
            formula=("INFO/SUPPORT < ${min_num_reads} | ABS(INFO/SVLEN) < ${min_length_sv} | QUAL < ${min_num_reads} | FILTER = \"hom_ref\" | FMT/GT[0] = \"./.\"")
        ;;
        INV)
            formula=("INFO/SUPPORT < ${min_num_reads} | ABS(INFO/SVLEN) < ${min_length_sv} | QUAL < ${min_num_reads} | FILTER = \"hom_ref\"")
        ;;
        DUP)
            formula=("INFO/SUPPORT < ${min_num_reads} | ABS(INFO/SVLEN) < ${min_length_sv} | QUAL < ${min_num_reads} | FILTER = \"hom_ref\" | FILTER = \"not_fully_covered\" | FMT/DHBFC[0] < 1.3")
        ;;
    esac
fi
bcftools view -e "${formula[*]}" "${vcf}" > "${vcf_filtered}"; } \
1> "${log}" 2>&1
