#! /bin/bash

# {{ input_vcf2maf="{input.vcf_extracted}"

# if [ "{wildcards.caller}" == "svim" ] && [ "{wildcards.type_sv}" == "DUP" ]; then
#     sed -e 's/DUP:TANDEM/DUP/g' -e 's/DUP:INT/DUP/g' ${{input_vcf2maf}} > ${{input_vcf2maf%.*}}.DUP.vcf
#     input_vcf2maf="${{input_vcf2maf%.*}}.DUP.vcf"
# elif [ "{wildcards.caller}" == "svision" ]; then
#     perl -pe 's/SVTYPE=.*?;/SVTYPE={wildcards.type_sv};/g' ${{input_vcf2maf}} > ${{input_vcf2maf%.*}}.{wildcards.type_sv}.vcf
#     input_vcf2maf="${{input_vcf2maf%.*}}.{wildcards.type_sv}.vcf"
# fi

# vcf2maf.pl --input-vcf ${{input_vcf2maf}} --output-maf {output.maf} --inhibit-vep --ncbi-build {params.genome} --cache-version {params.version_cache_vep} --ref-fasta {input.fasta} --vcf-tumor-id {wildcards.sample} --tumor-id {wildcards.sample} --vep-data {input.dir_db_vep} --species homo_sapiens

# common_fields="CHROM POS ID REF ALT QUAL \\
# ANN[*].ALLELE ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID \\
# ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE ANN[*].RANK ANN[*].HGVS_C \\
# ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].CDS_POS \\
# ANN[*].CDS_LEN ANN[*].AA_POS ANN[*].AA_LEN ANN[*].DISTANCE ANN[*].ERRORS \\
# LOF[*].GENE LOF[*].GENEID LOF[*].NUMTR LOF[*].PERC \\
# NMD[*].GENE NMD[*].GENEID NMD[*].NUMTR NMD[*].PERC"

# if [ "{wildcards.caller}" == "svim" ]; then
#     fmt_fields="GEN[*].GT GEN[*].DP GEN[*].AD"
# else
#     fmt_fields="GEN[*].GT GEN[*].DR GEN[*].DV"
# fi

# snpsift extractFields -s "," -e "." {input.vcf_extracted} ${{common_fields}} ${{fmt_fields}} > {output.tsv}

# sed -i '1s/GEN\[\*\]\.//g ; 1s/ANN\[\*\]\.//g ; 1s/\[\*\]//g' {output.tsv}; }} \\
# 1> {log} 2>&1
