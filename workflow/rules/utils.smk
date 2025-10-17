from pathlib import Path


def get_targets():
    targets = []

    targets += [
        f"survivor/{sample}/final/{sample}.{type_sv}.merged.vcf"
        for sample in SAMPLES
        for type_sv in TYPES_SV
    ]

    return targets


def get_annotsv_cache_outputs():
    if SPECIES in ["homo_sapiens"]:
        return {
            "dir_1": f"{config['cache_annotsv']}/Annotations_Human",
            "dir_2": f"{config['cache_annotsv']}/Annotations_Exomiser",
        }
    else:
        raise ValueError("Unsupported species")


def get_annotsv_cache_parameters():
    if SPECIES in ["homo_sapiens"]:
        return {
            "arg_install": "install-human-annotation",
            "dirs": [
                "share/AnnotSV/Annotations_Human",
                "share/AnnotSV/Annotations_Exomiser",
            ],
        }
    else:
        raise ValueError("Unsupported species")


def get_convert_snpeff_arguments(wildcards):
    caller = wildcards.caller

    fields = CALLER2FMTS.get(caller)
    if fields is None:
        raise ValueError("Unsupported caller")

    arg = " ".join(f"GEN[*].{field}" for field in fields)

    return arg


def get_format_svision_parameters(wildcards):
    sample = wildcards.sample

    suffixes = [
        f".{chrom}.svision.s{config['min_reads']}.graph.vcf" for chrom in CHROMS
    ]
    vcfs = [f"svision/{sample}/chroms/{sample}{suffix}" for suffix in suffixes]

    for index, (chrom, vcf) in enumerate(zip(CHROMS, vcfs)):
        with checkpoints.svision.get(sample=sample).output[index].open("r") as f:
            if Path(vcf).exists() and Path(vcf).stat().st_size > 0:
                return {
                    "chrom_lead": chrom,
                    "vcf_lead": vcf,
                }

    raise ValueError(f"No non-empty VCF found for sample {sample}.")
