#!/usr/bin/env bash
# shellcheck disable=SC2154

vcf="${snakemake_input[vcf]}"
bam="${snakemake_input[bam]}"
fasta="${snakemake_input[fasta]}"
vcf_filtered="${snakemake_output[vcf_filtered]}"
vcf_duphold="${snakemake_output[vcf_duphold]}"
vcf_duphold_filtered="${snakemake_output[vcf_duphold_filtered]}"

min_num_reads="${snakemake_params[min_num_reads]}"
min_length_sv="${snakemake_params[min_length_sv]}"
min_coverage="${snakemake_params[min_coverage]}"
threads="${snakemake[threads]}"
log="${snakemake_log[0]}"

sample="${snakemake_wildcards[sample]}"
caller="${snakemake_wildcards[caller]}"

function filter_sv() {
    if [ $# -lt 4 ]; then
        echo "Usage: filter_sv <vcf> <bam> <sample name> <caller name>" >&2
        return 1
    fi
    local vcf=$1
    local bam=$2
    local s=$3
    local t=$4
    if [ ! -f "${vcf}" ]; then
        echo "Error: VCF file ${vcf} not found." >&2
        return 1
    fi
    if [ ! -f "${bam}" ]; then
        echo "Error: BAM file ${bam} not found." >&2
        return 1
    fi
    if [ -z "${fasta}" ] || [ -z "${threads}" ] || [ -z "${vcf_filtered}" ] || [ -z "${vcf_duphold}" ] || [ -z "${vcf_duphold_filtered}" ] || [ -z "${min_num_reads}" ] || [ -z "${min_length_sv}" ] || [ -z "${min_coverage}" ]; then
        echo "Error: Required environment variables fasta, threads, vcf_filtered, vcf_duphold, vcf_duphold_filtered, min_num_reads, min_length_sv, or min_coverage are not set." >&2
        return 1
    fi

    if [ "${t}" == "cutesv" ]; then
        local formula="INFO/RE < ${min_num_reads} | (INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < ${min_length_sv}) | FMT/DV[0] < 3 | (FMT/DR[0] < 3 & FMT/DV[0] + FMT/DR[0] < 5)"
    elif [ "${t}" == "svim" ]; then
        local formula="INFO/SUPPORT < 3 | (INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < 50) | (SVTYPE = \"DEL,INS\" & FMT/AD[0:1] < 3 | (FMT/AD[0:0] < 3 & FMT/DP[0] < 5))"
    elif [ "${t}" == "severus" ]; then
        local formula="FMT/DV[0] < 3 | (INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < 50) | (FMT/DR[0] < 3 & FMT/DV[0] + FMT/DR[0] < 5)"
    elif [ "${t}" == "debreak" ]; then
        local formula="INFO/SUPPREAD < 3 | (INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < 50)"
    elif [ "${t}" == "nanovar" ]; then
        local formula="INFO/SR < 3 | (INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < 50) | FMT/AD[0:1] < 3 | (FMT/AD[0:0] < 3 & FMT/DP[0] < 5)"
    elif [ "${t}" == "delly" ]; then
        local formula="INFO/SR < 3 | (INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < 50) | FMT/RV[0] < 3 | (FMT/RR[0] < 3 & FMT/RV[0] + FMT/RR[0] < 5)"
    elif [ "${t}" == "nanosv" ]; then
        local formula="FMT/DV[0:0] < 3 | (INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < 50) | (SUM(FMT/DR[0:]) < 3 & FMT/DV[0:0] + SUM(FMT/DR[0:]) < 5)"
    else # Sniffles SVision
        local formula="INFO/SUPPORT < 3 | (INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < 50) | FMT/DV[0] < 3 | (FMT/DR[0] < 3 & FMT/DV[0] + FMT/DR[0] < 5)"
    fi

    bcftools sort -Ov "${vcf}" \
        | grep -Ev '^alt|^random|^fix|^chrM|^chrUn|^chrEBV|^HPV|^CMV|^HBV|^KSHV|^HTLV|^MCV|^SV40|^HIV|^HCV|^GL' - \
        | bcftools view -e "${formula}" - > "${vcf_filtered}"

    export DUPHOLD_SAMPLE_NAME=${s}
    duphold -v "${vcf_filtered}" -b "${bam}" -f "${fasta}" -t "${threads}" -o "${vcf_duphold}"
    unset DUPHOLD_SAMPLE_NAME

    bcftools view -e '(SVTYPE = "DEL" & FMT/DHFFC[0] > 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0] < 1.3)' "${vcf_duphold}" > "${vcf_duphold_filtered}"
}


filter_sv ${vcf} ${bam} ${sample} ${caller} 1> ${log} 2>&1

