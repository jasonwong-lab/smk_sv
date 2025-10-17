TYPES_SV = ["DEL", "DUP", "INV", "BND", "INS"]

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

FIELDS_COMMON = (
    "CHROM POS ID REF ALT QUAL "
    "ANN[*].ALLELE ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID "
    "ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE ANN[*].RANK ANN[*].HGVS_C "
    "ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].CDS_POS "
    "ANN[*].CDS_LEN ANN[*].AA_POS ANN[*].AA_LEN ANN[*].DISTANCE ANN[*].ERRORS "
    "LOF[*].GENE LOF[*].GENEID LOF[*].NUMTR LOF[*].PERC "
    "NMD[*].GENE NMD[*].GENEID NMD[*].NUMTR NMD[*].PERC"
)

CALLER2FMTS = {
    "cutesv": ["GT", "DR", "DV", "PL", "GQ"],
    "debreak": ["GT"],
    "delly": ["GT", "GL", "GQ", "FT", "RCL", "RC", "RCR", "RDCN", "DR", "DV", "RR", "RV"],
    "nanosv": ["GT", "DR", "DV", "GQ", "HR", "PL"],
    "nanovar": ["GT", "DP", "AD"],
    "severus": ["GT", "GQ", "VAF", "hVAF", "DR", "DV"],
    "sniffles": ["GT", "GQ", "DR", "DV"],
    "svim": ["GT", "DP", "AD"],
    "svision": ["GT", "DR", "DV"],
}
