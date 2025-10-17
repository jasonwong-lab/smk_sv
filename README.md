<!-- markdownlint-configure-file {"no-inline-html": {"allowed_elements": ["code", "details", "h2", "summary"]}} -->

# A snakemake pipeline to call structural variants from tumor-only ONT data

![License GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)

> :warning: This pipeline is in its early stages. Please use with caution.

## Author

Minghao Jiang, <jiang01@icloud.com>

## Tools Used

- SV callers

   [cuteSV](https://github.com/tjiangHIT/cuteSV), [Sniffles](https://github.com/fritzsedlazeck/Sniffles), [SVIM](https://github.com/eldariont/svim), [SVision](https://github.com/xjtu-omics/SVision), [Severus](https://github.com/KolmogorovLab/Severus), [NanoSV](https://github.com/mroosmalen/nanosv), [NanoVar](https://github.com/cytham/nanovar), [Delly](https://github.com/dellytools/delly), [Debreak](https://github.com/Maggi-Chen/DeBreak)

- Annotation tools

   [AnnotSV](https://github.com/lgmgeo/AnnotSV), [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), [SnpEff](http://pcingola.github.io/SnpEff/snpeff/introduction/)

- R packages

   [vroom](https://www.tidyverse.org/tags/vroom/), [tibble](https://tibble.tidyverse.org/reference/tibble-package.html), [glue](https://glue.tidyverse.org), [dplyr](https://dplyr.tidyverse.org), [tidyr](https://tidyr.tidyverse.org), [purrr](https://purrr.tidyverse.org), [GenomicRanges](https://github.com/Bioconductor/GenomicRanges), [stringr](https://stringr.tidyverse.org), [BiocParallel](https://github.com/Bioconductor/BiocParallel), parallel

- Other tools

   [Minimap2](https://github.com/lh3/minimap2), [SAMtools](https://github.com/samtools/samtools), [BCFtools](http://samtools.github.io/bcftools/bcftools.html), [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR), [vcf2maf](https://github.com/mskcc/vcf2maf), [SnpSift](http://pcingola.github.io/SnpEff/snpsift/introduction/), [duphold](https://github.com/brentp/duphold)

## Pipeline Structure

```mermaid
%%{
  init: {
    'theme': 'base',
    'themeVariables': {
      'fontFamily': 'Comic Sans MS',
      'primaryColor': '#ACD98DFF',
      'primaryTextColor': 'black',
      'lineColor': 'black',
      'secondaryColor': 'grey',
      'tertiaryColor': '#EEEEEE'
    }
  }
}%%

flowchart TD

  classDef norm stroke-width: 1px, stroke: black;
  classDef tool fill: #FFB977FF, padding: 0px, margin: 0px, stroke-width: 1px, stroke: black;
  classDef other_file stroke-dasharray: 5 5, stroke-width: 1px, stroke: black;
  classDef other_tool fill: #FFB977FF, padding: 0px, margin: 0px, stroke-dasharray: 5 5, stroke-width: 1px, stroke: black;
  classDef output fill: #98D9E4FF, stroke-width: 1px, stroke: black;

  fastq@{ shape: lean-r, label: "FASTQ" }
  bam([BAM])
  cutesv_vcf([VCF])
  sniffles2_vcf([VCF])
  cutesv_filtered_vcf([filtered
  VCF])
  sniffles2_filtered_vcf([filtered
  VCF])
  merged_vcf([merged VCF])
  annotated_vcf([annotated VCF/TSV])
  somatic_vcf@{ shape: lean-l, label: "somatic SVs" }
  germline_vcf@{ shape: lean-l, label: "germline SVs" }
  other_vcfs([VCFs])
  other_filtered_vcf([filtered
  VCFs])
  minimap2@{ shape: tag-rect, label: "Minimap2" }
  cutesv@{ shape: tag-rect, label: "cuteSV"  }
  sniffles2@{ shape: tag-rect, label: " Sniffles2 " }
  other_callers@{ shape: tag-rect, label: "other\ncallers" }
  bcftools@{ shape: tag-rect, label: "BCFtools" }
  survivor@{ shape: tag-rect, label: "SURVIVOR" }
  anno_tools@{ shape: tag-rect, label: "VEP, AnnotSV, SnpEff" }

  fastq:::norm
  bam:::norm
  cutesv_vcf:::norm
  sniffles2_vcf:::norm
  other_vcfs:::other_file
  cutesv_filtered_vcf:::norm
  sniffles2_filtered_vcf:::norm
  other_filtered_vcf:::other_file
  merged_vcf:::norm
  annotated_vcf:::norm
  somatic_vcf:::output
  germline_vcf:::output

  minimap2:::tool
  cutesv:::tool
  sniffles2:::tool
  other_callers:::other_tool
  bcftools:::tool
  survivor:::tool
  anno_tools:::tool

  fastq --> minimap2 --> bam
  bam --> cutesv --> cutesv_vcf --> bcftools --> cutesv_filtered_vcf --> survivor
  bam --> sniffles2 --> sniffles2_vcf --> bcftools --> sniffles2_filtered_vcf --> survivor
  bam -.-> other_callers -.-> other_vcfs -.-> bcftools -.-> other_filtered_vcf -.-> survivor
  survivor --> merged_vcf --> anno_tools --> annotated_vcf
  annotated_vcf --> somatic_vcf
  annotated_vcf --> germline_vcf
```

<details>

<summary><h2>Recommended Project Structure</h2></summary>

```text
project/
├── data/
│   └── wgs/
│       └── sample_xx.fq.gz
├── code/
│   └── wgs/
│       └── smk_sv/
└── analysis/
    └── wgs/
        ├── cutesv/
        |   └── sample_xx/
        |       └── merged/
        |           └── filtered/
        ├── sniffles/
        |   └── sample_xx/
        |       └── merged/
        |           └── filtered/
        └── survivor/
            └── sample_xx/
                └── final/
```

</details>

## Prerequisites

- [**Python**](https://www.python.org)
- [**Snakemake**](https://snakemake.github.io)
- [**eido**](https://pep.databio.org/eido/)
- [**SAMtools**](https://www.htslib.org)
- [**Mamba**](https://mamba.readthedocs.io/en/latest/) (recommended) or [**conda**](https://docs.conda.io/projects/conda/en/stable/)

Additional dependencies are automatically installed by **Mamba** or **conda**. Environments are defined in yaml files under `workflow/envs/`.

```shell
# ---------------------------------------------------------------------------- #
# Install Mamba and SAMtools manually. Since conda-packaged SAMtools           #
# occasionally encounters issues, this workflow presumes that samtools is      #
# executable within your system's PATH.                                        #
# ---------------------------------------------------------------------------- #
if ! command -v mamba &> /dev/null; then
    "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
    source ~/.bashrc
fi

if ! command -v samtools &> /dev/null; then
    [ -d ${HOME}/.local/opt ] || mkdir -p ${HOME}/.local/opt
    wget 'https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2'
    tar -xvf samtools-1.22.1.tar.bz2
    cd samtools-1.22.1
    ./configure --prefix=${HOME}/.local/opt/samtools
    make
    make install
    echo "export PATH=\${HOME}/.local/opt/samtools/bin:\${PATH}" >> ~/.bashrc
    source ~/.bashrc
    rm -rf samtools-1.22.1 samtools-1.22.1.tar.bz2
fi

# Install Snakemake and eido using pipx (https://pipx.pypa.io/stable/)
pipx install snakemake
pipx inject snakemake eido

# Clone the repository
git clone https://github.com/jasonwong-lab/smk_sv.git
cd smk_sv/

# Initialize configuration
cp config/.config.yaml config/config.yaml
cp config/pep/.config.yaml config/pep/config.yaml
cp workflow/profiles/default/.config.yaml workflow/profiles/default/config.yaml
```

## Configuration

### Main Configuration

<details>

<summary>Edit <code>config/config.yaml</code></summary>

```yaml
dir_run: /projects/SV/analysis/wgs
mapper: minimap2
callers:
  - cutesv
  - severus
  - sniffles
  - svim
  - svision
annotators:
  - vep
  - snpeff
  - annotsv
fasta: /doc/reference/fasta/GRCh37.primary_assembly.genome.fa
index_minimap2: /doc/reference/minimap2/GRCh37.primary_assembly.genome.mmi
dir_data: /projects/SV/data/wgs
suffix_fastq: .fq.gz
bed_tandem_repeats: /doc/sniffles/human_hs37d5.trf.chr.bed
min_reads: 3
min_length_reads: 1000
min_quality_mapping: 20
min_coverage: 6
min_size: 100
max_size: 10000000
min_dhffc: 0.7
max_dhbfc: 1.3
distance_sv:
  BND: 10
  DEL: 10
  INS: 10
  INV: 10
  DUP: 10
n_callers:
  BND: 3
  DEL: 3
  INS: 3
  INV: 3
  DUP: 3
consider_type:
  BND: false
  DEL: false
  INS: false
  INV: false
  DUP: false
consider_strand:
  BND: false
  DEL: false
  INS: false
  INV: false
  DUP: false
estimate_distance:
  BND: true
  DEL: true
  INS: true
  INV: true
  DUP: true
terms_relative: leuka?emia|blood|lymph|myelo|ha?ema|marrow|platel|thrombo|anemia|neutro
species: homo_sapiens
genome: GRCh37
version_snpeff: "87"
version_vep: 114
version_annotsv: v3.5
cache_snpeff: /doc/snpeff
cache_vep: /.vep
cache_annotsv: /doc/tool/annotator/annotsv
max_size_vep: 10000000
config_nanosv: /doc/nanosv/config.ini
bed_nanosv: /opt/nanosv/nanosv/bedfiles/hg19_genome_sample.bed
model_clair3: /doc/clair3/models/r941_prom_sup_g5014
bed_nvtr: /doc/sniffles/human_hs37d5.trf.chr.bed
model_svision: /doc/svision/svision-cnn-model.ckpt
```

</details>

### Execution Profile

<details>

<summary>Edit <code>workflow/profiles/default/config.yaml</code></summary>

```yaml
software-deployment-method:
  - apptainer
  - conda
apptainer-args: "--bind /"
conda-prefix: /.snakemake/envs/smk_sv
scheduler: greedy
rerun-trigger: mtime
printshellcmds: True
keep-incomplete: True
cores: all
resources:
  mem_mb: 500000  # 500GB
default-resources:
  mem_mb: 10000  # 10GB
set-threads:
  annotate_sv_snpeffnvep: 10
  annotsv: 10
  filter_sv_annotation: 10
  filter_sv: 10
  clair3: 10
  cutesv: 10
  debreak: 10
  nanosv: 10
  nanovar: 10
  severus: 10
  snpeff: 10
  sniffles: 10
  svision: 10
  minimap2: 10
  minimap2_index: 10
  vep: 10
set-resources:
  snpeff:
    mem_mb: 50000
  svision:
    runtime: 72h
```

</details>

### Sample Metadata

This workflow uses [**Portable Encapsulated Projects (PEP)**](https://pep.databio.org/) for sample management.

<details>

<summary>Edit <code>config/pep/config.yaml</code></summary>

```yaml
pep_version: 2.1.0
sample_table: samples.csv    # Path to the sample table (Required)
```

</details>

The sample table must include one mandatory column:

| **sample_name**                   |
| --------------------------------- |
| Unique identifier for each sample |

## Execution

### Local Execution

```shell
# Create environments
snakemake --conda-create-envs-only

# Run the workflow
snakemake
```

### Cluster Execution

If you want to run this pipeline on a cluster (*e.g.*, SLURM, or PBS), you should customise your own profile and place it into `~/.config/snakemake/`, and then run the pipeline with the profile:

```shell
snakemake --profile <your_profile_name>
```

Or run the pipeline with the profile you have set as an environment variable:

```shell
export SNAKEMAKE_PROFILE=<your_profile_name>
snakemake
```

You can refer to the profile I have been using (`workflow/profiles/mycluster`):

```shell
[ -d ~/.config/snakemake ] || mkdir -p ~/.config/snakemake
ln -s `pwd`/workflow/profiles/mycluster ~/.config/snakemake/
```

## License

Codes here are licensed under the [GNU General Public License v3](http://www.gnu.org/licenses/gpl-3.0.html).
