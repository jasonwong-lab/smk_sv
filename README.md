# A snakemake pipeline to call structural variants from ONT data

## Author

Minghao Jiang, <jiang01@icloud.com>

## Supported tools

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
   - [SnpEff](http://pcingola.github.io/SnpEff/)

## Pipeline structure

```mermaid
---
title: SV calling workflow
---
flowchart TD
		fastq([FASTQ]) -- minimap2 --> bam([BAM])
		bam([BAM]) -- cuteSV --> cutesv_vcf([VCF])
		bam([BAM]) -- Sniffles2 --> sniffles2_vcf([VCF])
		bam([BAM]) -- SVIM --> svim_vcf([VCF])
		cutesv_vcf([VCF]) -- bcftools --> cutesv_filtered_vcf([filtered VCF])
		sniffles2_vcf([VCF]) -- bcftools --> sniffles2_filtered_vcf([filtered VCF])
		svim_vcf([VCF]) -- bcftools --> svim_filtered_vcf([filtered VCF])
		cutesv_filtered_vcf([filtered VCF]) -- SURVIVOR --> merged_vcf([merged VCF])
		sniffles2_filtered_vcf([filtered VCF]) -- SURVIVOR --> merged_vcf([merged VCF])
		svim_filtered_vcf([filtered VCF]) -- SURVIVOR --> merged_vcf([merged VCF])
		merged_vcf([merged VCF]) -- VEP, AnnotSV, and SnpEff --> annotated_vcf([annotated VCF])
		annotated_vcf([annotated VCF]) -.-> somatic_vcf([somatic SVs])
		annotated_vcf([annotated VCF]) -.-> germline_vcf([germline SVs])
		fastq:::myclass; bam:::myclass; cutesv_vcf:::myclass; sniffles2_vcf:::myclass; svim_vcf:::myclass; cutesv_filtered_vcf:::myclass; sniffles2_filtered_vcf:::myclass; svim_filtered_vcf:::myclass; merged_vcf:::myclass; annotated_vcf:::myclass; somatic_vcf:::myclass; germline_vcf:::myclass
		classDef myclass fill: #A1DD70, stroke:#EE4E4E
```

## Usage

1. **Ensure you have clonned this repo and navigated to the directory `workflow`.**

   *Note: All steps below should be followed after you are in the `workflow` dir.*

2. **Build the `singularity` or `apptainer` image from the Docker Hub** (using singularity as an example here):

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

   **Or build the image from the def file:**

   ```shell
   mkdir singularities
   singularity build --sandbox singularities/sv scripts/sv.def
   ```

   You might need `--fakeroot` to run `singularity build` successfully.

3. ~~There are some tools not included in the built docker image, e.g., [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), [SnpEff](http://pcingola.github.io/SnpEff/), [vcf2maf](https://github.com/mskcc/vcf2maf), and some R packages, since one might have particular resources for annotating already. Thus, **you should install these missing tools manually**. (*Unfinished*)~~

   ~~*Note: R packages will be installed automatically when running the pipeline if needed*~~

   *They are included in the image now. However, the image size would be larger (~ 9GB).*

   *Note: When you prefer using a different version of VEP, please add `container: None` into the rule `annotate_sv`. Don't forget to make `vep` executable in your environment.*

4. **You should install [AnnotSV](https://github.com/lgmgeo/AnnotSV) by yourself**, as it's not included in the image due to its large annotation resources (~ 20GB) that cannot be specified elsewhere.

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