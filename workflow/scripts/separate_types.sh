#!/usr/bin/env bash
# shellcheck disable=SC2154

set -x

{ vcf="${snakemake_input[vcf]}"
vcf_separated="${snakemake_output[vcf]}"
max_size="${snakemake_params[max_size]}"

caller="${snakemake_wildcards[caller]}"
type_sv="${snakemake_wildcards[type_sv]}"


if { [ "${caller}" == "cutesv" ] || [ "${caller}" == "nanovar" ]; } && [ "${type_sv}" == "BND" ]; then

    bcftools filter -i 'SVTYPE ~ "BND"' "${vcf}" \
        | sed -e 's/END=[0-9]\+;//g' \
        | awk '!/[\[\]]chr[1-22XYM]:0[\[\]]/' > "${vcf_separated}"

elif [ "${caller}" == "cutesv" ] && [ "${type_sv}" == "DEL" ]; then

    bcftools filter -i "SVTYPE ~ \"DEL\" & ABS(INFO/SVLEN) <= ${max_size}" "${vcf}" > "${vcf_separated}"

elif [ "${caller}" == "svim" ] && [ "${type_sv}" == "DUP" ]; then

    bcftools filter -i '(SVTYPE ~ "DUP:INT") | (SVTYPE ~ "DUP:TANDEM") | (SVTYPE ~ "DUP")' "${vcf}" > "${vcf_separated}"

elif [ "${caller}" == "svim" ] && [ "${type_sv}" == "INV" ]; then

    bcftools filter -i 'SVTYPE ~ "INV"' "${vcf}" \
        | awk -F'\t' -v OFS='\t' '{
            if ($0 ~ /^#/) {print $0;} else {
                split($8, INFO, ";");
                split(INFO[2], array, "=");
                sub(/;/, ";SVLEN="array[2] - $2";", $8);
                $5="<INV>";
                print $0
            }
        }' > "${vcf_separated}"

elif [ "${caller}" == "debreak" ] && [ "${type_sv}" == "BND" ]; then

    bcftools filter -i 'SVTYPE ~ "TRA"' "${vcf}" \
        | sed -e 's/SVTYPE=TRA/SVTYPE=BND/g' \
        | awk '/^##ALT/ && !f {print "##ALT=<ID=BND,Description=\"BND\">"; f=1} 1' \
        | awk 'BEGIN {OFS=FS="\t"} !/^#/ {$5 = "<BND>"} 1' > "${vcf_separated}"

else

    bcftools filter -i "SVTYPE ~ \"${type_sv}\"" "${vcf}" > "${vcf_separated}"

fi; } \
1> "${snakemake_log[0]}" 2>&1
