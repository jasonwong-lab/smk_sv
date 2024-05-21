#!/usr/bin/env bash
# shellcheck disable=SC2154

vcf="${snakemake_input[vcf]}"
fasta="${snakemake_input[fasta]}"
dir_db_snpeff="${snakemake_input[dir_db_snpeff]}"
dir_db_vep="${snakemake_input[dir_db_vep]}"
vcf_snpeff="${snakemake_output[vcf_snpeff]}"
vcf_vep="${snakemake_output[vcf_vep]}"

genome="${snakemake_params[genome]}"
version_cache_snpeff="${snakemake_params[version_cache_snpeff]}"
version_cache_vep="${snakemake_params[version_cache_vep]}"
threads="${snakemake[threads]}"
log="${snakemake_log[0]}"

type_sv="${snakemake_wildcards[type_sv]}"
caller="${snakemake_wildcards[caller]}"

function annotate_sv() {
    if [ $# -lt 3 ]; then
        cat << EOF >&2
Usage: annotate_sv <vcf> <sv type> <caller name> [<genome version>] [<version snpeff>] [<version vep>]
    <vcf>: The input VCF file.
    <sv type>: The type of structural variant. Must be one of BND, DEL, DUP, INS, INV.
    <caller name>: The name of caller.
    [<genome version>]: The genome version. Defaults to GRCh38.
    [<version snpeff>]: The cache version of snpeff. Defaults to 105.
    [<version vep>]: The cache version of vep. Defaults to 110.
EOF
        return 1
    fi
    local vcf=$1
    local v=$2
    local t=$3
    local genome_version=${4:-GRCh38}
    local version_snpeff=${5:-105}
    local version_vep=${6:-110}
    if [[ ! "BND DEL DUP INS INV" =~ ${v} ]]; then
        echo "Error: Invalid SV type. Must be one of BND, DEL, DUP, INS, INV." >&2
        return 1
    fi
    if [ -z "${dir_db_snpeff}" ] || [ -z "${dir_db_vep}" ] || [ -z "${fasta}" ] || [ -z "${threads}" ] || [ -z "${vcf_snpeff}" ] || [ -z "${vcf_vep}" ]; then
        echo "Error: Required environment variables dir_db_snpeff, dir_db_vep, fasta, threads, vcf_snpeff, or vcf_vep are not set." >&2
        return 1
    fi

    local l_input
    l_input=$(wc -l "${vcf}" | awk '{print $1}')
    local l_output_s=0
    local l_output_v=0
    [ -f "${vcf_snpeff}" ] && l_output_s=$(wc -l "${vcf_snpeff}" | awk '{print $1}')
    [ -f "${vcf_vep}" ] && l_output_v=$(wc -l "${vcf_vep}" | awk '{print $1}')
    local size_max=300
    local size_input
    size_input=$(du -m "${vcf}" | cut -f1)

    if [ "${l_input}" -ne $((l_output_s - 5)) ]; then
        local input_snpeff
        input_snpeff="${vcf}"
        if [ "${t}" == "svision" ] && [ "${v}" == "INS" ]; then
            awk 'BEGIN {OFS=FS="\t"} !/^#/ && $5 == "<INS>" {$5 = "N"} 1' "${input_snpeff}" > "${input_snpeff%.*}".5.vcf
            input_snpeff="${input_snpeff%.*}".5.vcf
        elif [ "${t}" == "svision" ] && [ "${v}" != "INS" ]; then
            awk -v var="$v" 'BEGIN {OFS=FS="\t"} !/^#/ {$5 = "<"var">"} 1' "${input_snpeff}" > "${input_snpeff%.*}".5.vcf
            input_snpeff="${input_snpeff%.*}".5.vcf
        elif [ "${t}" == "debreak" ] && [ "${v}" == "DEL" ] && [ "${size_input}" -ge ${size_max} ]; then
            awk 'BEGIN {OFS=FS="\t"} !/^#/ {$4 = "N"; $5 = "<DEL>"} 1' "${input_snpeff}" > "${input_snpeff%.*}".4_5.vcf
            input_snpeff="${input_snpeff%.*}".4_5.vcf
        fi
        snpeff -Xmx81920M "${genome_version}"."${version_snpeff}" -nodownload -canon -v -lof -s "${vcf_snpeff%.*}".html -dataDir "${dir_db_snpeff}" "${input_snpeff}" > "${vcf_snpeff}"
    else
        echo "[INFO] SnpEff has annotated this VCF..."
    fi

    if [ "${l_input}" -ne $((l_output_v - 8)) ]; then
        local input_vep
        input_vep="${vcf_snpeff}"
        if [ "${t}" == "svision" ] && [ "${v}" == "INS" ] ; then
            awk 'BEGIN {FS=OFS="\t"} $5 == "N" {$5 = "<INS>"} 1' "${input_vep}" > "${input_vep%.*}".5.vcf
            input_vep="${input_vep%.*}".5.vcf
        elif [ "${t}" == "cutesv" ] && [ "${v}" == "DEL" ] && [ "${size_input}" -ge ${size_max} ]; then
            awk 'BEGIN {OFS=FS="\t"} !/^#/ {$4 = "N"; $5 = "<DEL>"} 1' "${input_vep}" > "${input_vep%.*}".4_5.vcf
            input_vep="${input_vep%.*}".4_5.vcf
        fi
        local vep_command=(vep -i "${input_vep}" -o "${vcf_vep}" --vcf --everything --filter_common --per_gene --total_length --offline --format vcf --assembly "${genome_version}" --cache_version "${version_vep}" --species homo_sapiens --fasta "${fasta}" --cache --dir_cache "${dir_db_vep}" --fork "${threads}" --max_sv_size 999999999999999999 --stats_file "${vcf_vep%.*}".html --force_overwrite)
        if [ "${v}" == BND ]; then
            "${vep_command[@]}" --buffer_size 50
        else
            "${vep_command[@]}"
        fi
    else
        echo "[INFO] VEP has annotated this VCF..."
    fi
}


annotate_sv "${vcf}" "${type_sv}" "${caller}" "${genome}" "${version_cache_snpeff}" "${version_cache_vep}" 1> "${log}" 2>&1
