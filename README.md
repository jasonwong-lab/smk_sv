# A snakemake pipeline to call structural variants from ONT data

## Author

Minghao Jiang, <jiang01@icloud.com>

## Tools used

- SV callers
   - [cuteSV](https://github.com/tjiangHIT/cuteSV)
   - [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
   - [SVIM](https://github.com/eldariont/svim)
   - [SVision](https://github.com/xjtu-omics/SVision)
   - [Severus](https://github.com/KolmogorovLab/Severus)
   - [NanoSV](https://github.com/mroosmalen/nanosv)
   - [NanoVar](https://github.com/cytham/nanovar)
   - [Delly](https://github.com/dellytools/delly)
   - [Debreak](https://github.com/Maggi-Chen/DeBreak)
- Annotation tools
   - [AnnotSV](https://github.com/lgmgeo/AnnotSV)
   - [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)
   - [SnpEff](http://pcingola.github.io/SnpEff/snpeff/introduction/)
- Other tools
   - [Minimap2](https://github.com/lh3/minimap2)
   - [SAMtools](https://github.com/samtools/samtools)
   - [BCFtools](http://samtools.github.io/bcftools/bcftools.html)
   - [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
   - [vcf2maf](https://github.com/mskcc/vcf2maf)
   - [SnpSift](http://pcingola.github.io/SnpEff/snpsift/introduction/)
- R packages
   - [vroom](https://www.tidyverse.org/tags/vroom/)
   - [tibble](https://tibble.tidyverse.org/reference/tibble-package.html)
   - [glue](https://glue.tidyverse.org)
   - [dplyr](https://dplyr.tidyverse.org)
   - [tidyr](https://tidyr.tidyverse.org)
   - [purrr](https://purrr.tidyverse.org)

## Pipeline structure

```mermaid
%%{
  init: {
    'theme': 'base',
    'themeVariables': {
      'fontFamily': 'Comic Sans MS',
      'primaryColor': '#F4CE14',
      'primaryTextColor': '#FFFFFF',
      'lineColor': '#EE6983',
      'secondaryColor': '#00B8A9',
      'tertiaryColor': '#FFFFFF'
    }
  }
}%%

---
title: SV calling workflow
---

flowchart TD

  classDef myclass fill:#00B8A9, stroke-width:0px, padding:0px, margin:0px;
  classDef myclass2 fill:#A5DD9B, stroke-dasharray:5 5;

  fastq([FASTQ]) -- Minimap2 --> bam([BAM])
  bam([BAM]) -- cuteSV --> cutesv_vcf([VCF])
  bam([BAM]) -- Sniffles2 --> sniffles2_vcf([VCF])
  cutesv_vcf([VCF]) -- BCFtools --> cutesv_filtered_vcf([filtered VCF])
  sniffles2_vcf([VCF]) -- BCFtools --> sniffles2_filtered_vcf([filtered VCF])
  cutesv_filtered_vcf([filtered VCF]) --- survivor["SURVIVOR"]
  sniffles2_filtered_vcf([filtered VCF]) --- survivor["SURVIVOR"]
  survivor["SURVIVOR"] --> merged_vcf([merged VCF])
  merged_vcf([merged VCF]) -- VEP, AnnotSV, and SnpEff --> annotated_vcf([annotated VCF/TSV])
  annotated_vcf([annotated VCF/TSV]) --> somatic_vcf([somatic SVs])
  annotated_vcf([annotated VCF/TSV]) --> germline_vcf([germline SVs])

  bam([BAM]) -. "other callers (e.g. SVIM)" .-> other_vcfs([VCFs])
  other_vcfs([VCFs]) -. BCFtools .-> other_filtered_vcf([filtered VCFs])

  other_filtered_vcf([filtered VCFs]) -.- survivor["SURVIVOR"]

  survivor:::myclass
  other_vcfs:::myclass2
  other_filtered_vcf:::myclass2
```


## Usage

1. **Ensure you have clonned this repo and navigated to the directory `workflow`.**

   *Note: All steps below should be followed after you are in the `workflow` dir.*

2. **Build a `apptainer` sandbox from the Docker image** (using singularity as an example here):

   ```shell
   mkdir singularities
   singularity pull singularities/sv.sif docker://mhjiang97/sv:latest
   singularity build --sandbox singularities/sv singularities/sv.sif
   ```

   Or:

   ```shell
   mkdir singularities
   singularity build --sandbox singularities/sv docker://mhjiang97/sv:latest
   ```

   **Or build it from the def file** (You might need `--fakeroot` to build from the singularity def file successfully.):

   ```shell
   mkdir singularities
   singularity build --sandbox singularities/sv scripts/sv.def
   ```

   *Note: A Dockerfile is also provided in the directory `scripts`.*

3. ~~There are some tools not included in the built docker image, e.g., [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), [SnpEff](http://pcingola.github.io/SnpEff/), [vcf2maf](https://github.com/mskcc/vcf2maf), and some R packages, since one might have particular resources for annotating already. Thus, **you should install these missing tools manually**. (*Unfinished*)~~

   ~~*Note: R packages will be installed automatically when running the pipeline if needed*~~

   *They are included in the image now. However, the image size would be larger (~ 9GB).*

   *Note: When you prefer using a different version of VEP, please add `container: None` into the rule `annotate_sv`. Don't forget to make `vep` executable in your environment.*

4. **You should install [AnnotSV](https://github.com/lgmgeo/AnnotSV) by yourself**, as it's not included in the image due to its large annotation resources (~ 20GB) that cannot be specified elsewhere.

   *Note: Creating a lock file for each combination of sample/type_sv has been implemented. However, AnnotSV might still encounter errors since it doesn’t support processing multiple files within the same directory. To address this, you can add `threads: workflow.cores` to the rule `annotate_sv_annotsv` to ensure that only one instance of this rule runs at a time.*

5. **Modify the `../config/config.yaml`** to specify needed file paths.

   *Note: You must change the file paths specified in the config.*

6. **Modify the column `sample_name` of `../config/pep/samples.csv`.**

   *Note: Only `sample_name` in the table will be used.*

   *More information please see [Portable Encapsulated Projects (PEP)](https://pep.databio.org).*

7. **Modify the `profiles/default/config.yaml`** to:

   - bind directories you need in the container.
   - change the number of CPUs you prefer.
   - modify/add/delete other parameters of this snakemake pipeline.

8. **Run the whole pipeline**:

   ```shell
   snakemake
   ```

   If you want to run this pipeline on a cluster (e.g., SLURM, or PBS), you should customise your own profile and place it into `~/.config/snakemake/`, and then run the pipeline:

   ```shell
   snakemake --profile <your_profile_name>
   ```

   Or run the pipeline:

   ```shell
   export SNAKEMAKE_PROFILE=<your_profile_name>
   snakemake
   ```

   You can refer to the profile I have been using at `profiles/mycluster`, or turn to snakemake websites.

## License

Codes here are licensed under the [GNU General Public License v3](http://www.gnu.org/licenses/gpl-3.0.html).
