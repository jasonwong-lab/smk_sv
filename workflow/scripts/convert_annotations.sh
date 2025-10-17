#!/usr/bin/env bash
# shellcheck disable=SC2154

set -x

vcf_extracted="${snakemake_input[vcf_extracted]}"
fasta="${snakemake_input[fasta]}"
cache_vep="${snakemake_input[cache_vep]}"
maf="${snakemake_output[maf]}"
tsv="${snakemake_output[tsv]}"

genome="${snakemake_params[genome]}"
version_vep="${snakemake_params[version_vep]}"

log="${snakemake_log[0]}"

sample="${snakemake_wildcards[sample]}"
caller="${snakemake_wildcards[caller]}"
type_sv="${snakemake_wildcards[type_sv]}"


{ input_vcf2maf=${vcf_extracted}

if [ "${caller}" == "svim" ] && [ "${type_sv}" == "DUP" ]; then
    sed -e 's/DUP:TANDEM/DUP/g' -e 's/DUP:INT/DUP/g' "${input_vcf2maf}" > "${input_vcf2maf%.*}".DUP.vcf
    input_vcf2maf=${input_vcf2maf%.*}.DUP.vcf
elif [ "${caller}" == "svision" ]; then
    perl -pe "s/SVTYPE=.*?;/SVTYPE=${type_sv};/g" "${input_vcf2maf}" > "${input_vcf2maf%.*}"."${type_sv}".vcf
    input_vcf2maf="${input_vcf2maf%.*}"."${type_sv}".vcf
fi

vcf2maf.pl --input-vcf "${input_vcf2maf}" --output-maf "${maf}" --inhibit-vep --ncbi-build "${genome}" --cache-version "${version_vep}" --ref-fasta "${fasta}" --vcf-tumor-id "${sample}" --tumor-id "${sample}" --vep-data "${cache_vep}" --species homo_sapiens

common_fields="CHROM POS ID REF ALT QUAL \
ANN[*].ALLELE ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID \
ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE ANN[*].RANK ANN[*].HGVS_C \
ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].CDS_POS \
ANN[*].CDS_LEN ANN[*].AA_POS ANN[*].AA_LEN ANN[*].DISTANCE ANN[*].ERRORS \
LOF[*].GENE LOF[*].GENEID LOF[*].NUMTR LOF[*].PERC \
NMD[*].GENE NMD[*].GENEID NMD[*].NUMTR NMD[*].PERC"

if [ "${caller}" == "svim" ]; then
    fmt_fields="GEN[*].GT GEN[*].DP GEN[*].AD"
else
    fmt_fields="GEN[*].GT GEN[*].DR GEN[*].DV"
fi

snpsift extractFields -s "," -e "." "${vcf_extracted}" "${common_fields}" "${fmt_fields}" \
    | grep -v 'Duplicate cpuset controllers detected' \
    | sed '1s/GEN\[\*\]\.//g ; 1s/ANN\[\*\]\.//g ; 1s/\[\*\]//g' > "${tsv}"; } \
1> "${log}" 2>&1
