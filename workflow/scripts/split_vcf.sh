#!/usr/bin/env bash
# shellcheck disable=SC2154

vcf="${snakemake_input[vcf]}"
vcf_splitted="${snakemake_output[vcf_splitted]}"

log="${snakemake_log[0]}"
caller="${snakemake_wildcards[caller]}"
type_sv="${snakemake_wildcards[type_sv]}"

function split_vcf() {
    if [ $# -lt 3 ]; then
        echo "Usage: split_vcf <vcf> <sv type> <caller name>" >&2
        return 1
    fi
    local vcf=$1
    local v=$2
    local t=$3
    if [ ! -f "${vcf}" ]; then
        echo "Error: VCF file ${vcf} not found." >&2
        return 1
    fi
    if [[ ! "BND DEL DUP INS INV" =~ ${v} ]]; then
        echo "Error: Invalid SV type. Must be one of BND, DEL, DUP, INS, INV." >&2
        return 1
    fi
    if [ -z "${vcf_splitted}" ]; then
        echo "Error: Required environment variable vcf_output is not set." >&2
        return 1
    fi

    if { [ "${t}" == "cutesv" ] || [ "${t}" == "nanovar" ]; } && [ "${v}" == "BND" ]; then
        bcftools view -i 'SVTYPE ~ "BND"' "${vcf}" \
            | sed -e 's/END=[0-9]\+;//g' > "${vcf_splitted}"
    elif [ "${t}" == "svim" ] && [ "${v}" == "DUP" ]; then
        bcftools view -i '(SVTYPE ~ "DUP:INT") | (SVTYPE ~ "DUP:TANDEM") | (SVTYPE ~ "DUP")' "${vcf}" > "${vcf_splitted}"
    elif [ "${t}" == "svim" ] && [ "${v}" == "INV" ]; then
        bcftools view -i 'SVTYPE ~ "INV"' "${vcf}" \
            | awk -F'\t' -v OFS='\t' '{
                    if ($0 ~ /^#/) {print $0;} else {
                        split($8, INFO, ";");
                        split(INFO[2], array, "=");
                        sub(/;/, ";SVLEN="array[2] - $2";", $8);
                        $5="<INV>";
                        print $0
                    }
                }' > "${vcf_splitted}"
    elif [ "${t}" == "debreak" ] && [ "${v}" == "BND" ]; then
        bcftools view -i 'SVTYPE ~ "TRA"' "${vcf}" \
            | sed -e 's/SVTYPE=TRA/SVTYPE=BND/g' \
            | awk '/^##ALT/ && !f {print "##ALT=<ID=BND,Description=\"BND\">"; f=1} 1' \
            | awk 'BEGIN {OFS=FS="\t"} !/^#/ {$5 = "<BND>"} 1' > "${vcf_splitted}"
    else
        bcftools view -i "SVTYPE ~ \"${v}\"" "${vcf}" > "${vcf_splitted}"
    fi
}


split_vcf ${vcf} ${type_sv} ${caller} 1> {log} 2>&1

