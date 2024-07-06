#! /bin/bash
# shellcheck disable=SC2154

set -x

vcf="${snakemake_input[vcf]}"
tsv="${snakemake_output[tsv]}"

genome="${snakemake_params[genome]}"
threads="${snakemake[threads]}"
log="${snakemake_log[0]}"

sample="${snakemake_wildcards[sample]}"
type_sv="${snakemake_wildcards[type_sv]}"
caller="${snakemake_wildcards[caller]}"


# sleep $((RANDOM % 10 + 1))

# lockfile=annotsv.${sample}.${caller}.lock
# while [ -f "${lockfile}" ]; do
#     sleep 10
# done
# touch "${lockfile}"

input_annotsv=${vcf}
size_max=300
size_input=$(du -m "${vcf}" | cut -f1)
if { [ "${caller}" == "cutesv" ] && [ "${type_sv}" == "DEL" ]; } \
    || { [ "${caller}" == "debreak" ] && [ "${type_sv}" == "DEL" ]; } \
    && [ "${size_input}" -ge ${size_max} ]; then
    awk 'BEGIN {OFS=FS="\t"} !/^#/ {$4 = "N"; $5 = "<DEL>"} 1' "${input_annotsv}" > "${input_annotsv%.*}".4_5.vcf
    input_annotsv="${input_annotsv%.*}".4_5.vcf
fi
AnnotSV -genomeBuild "${genome}" -SVinputFile "${input_annotsv}" -outputFile "${tsv}" -SVminSize 1 -overwrite 1 1> "${log}" 2>&1

# sleep 5
# rm "${lockfile}"
