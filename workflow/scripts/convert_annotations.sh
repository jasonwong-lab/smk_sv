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


function modify_vcf2maf() {
    file_vcf2maf=$(which vcf2maf.pl)
    line_410=$(sed -n '410p' "${file_vcf2maf}")
    [[ ! ${line_410} =~ ^# ]] && sed -i '410s/^/# /' "${file_vcf2maf}"
}

function download_snpeff() {
    if command -v snpsift &> /dev/null && [ -x "$(command -v snpsift)" ]; then
        echo "[INFO] SnpEff is already installed."
    else
        if [ ! -d snpEff ]; then
            lock_file="snpeff_download.lock"
            while ! mkdir "${lock_file}" 2>/dev/null; do
                echo "[INFO] Another instance is downloading SnpEff. Waiting..."
                sleep 5
            done
            if [ -d snpEff ]; then
                rmdir "${lock_file}"
                echo "[INFO] SnpEff is already downloaded."
                return
            fi
            echo "[INFO] Downloading SnpEff..."
            max_attempts=3
            attempt=0
            success=0
            while [ $attempt -lt $max_attempts ] && [ $success -eq 0 ]; do
                if curl -o snpEff_latest_core.zip 'https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip'; then
                    success=1
                    unzip snpEff_latest_core.zip && rm -f snpEff_latest_core.zip
                else
                    echo "[WARNING] Failed to download SnpEff, attempt $((attempt + 1)) of $max_attempts."
                    attempt=$((attempt + 1))
                    sleep 10 # wait before retrying
                fi
            done
            if [ $success -eq 0 ]; then
                echo "[ERROR] Failed to download SnpEff after $max_attempts attempts. Cleaning up..."
                [[ -f snpEff_latest_core.zip ]] && rm -f snpEff_latest_core.zip
                [[ -d snpEff ]] && rm -rf snpEff
                rmdir "${lock_file}"
                exit 1
            fi
            rmdir "${lock_file}"
        else
            echo "[INFO] SnpEff is already downloaded."
        fi
    fi
}


{ modify_vcf2maf
download_snpeff
# shellcheck disable=SC2155
export PATH="${PATH}:$(pwd)/snpEff/exec"

input_vcf2maf=${vcf_extracted}

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
