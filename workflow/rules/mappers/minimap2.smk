rule make_index_minimap2:
    conda:
        "../../envs/minimap2.yaml"
    input:
        fasta=config["fasta"],
    output:
        index_minimap2=protected(config["index_minimap2"]),
    threads: config["threads"]
    shell:
        "minimap2 -t {threads} -d {output.index_minimap2} {input.fasta}"


rule map_minimap2:
    conda:
        "../../envs/minimap2.yaml"
    input:
        index_minimap2=ancient(rules.make_index_minimap2.output.index_minimap2),
    output:
        bam=protected("minimap2/{sample}/{sample}.sorted.bam"),
        bai=protected("minimap2/{sample}/{sample}.sorted.bam.csi"),
    params:
        dir_data=config["dir_data"],
        suffix_fastq=config["suffix_fastq"],
    threads: math.floor(workflow.cores / NUMBER_SAMPLES)  # config["threads"]
    log:
        "logs/{sample}/map_minimap2.log",
    shell:
        """
        {{ minimap2 \\
        -ax map-ont -Y -t {threads} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:DNA\\tPL:ONT\" \\
        {input.index_minimap2} {params.dir_data}/{wildcards.sample}.{params.suffix_fastq} \\
            | samtools view -Sb -@ {threads} - \\
            | samtools sort -@ {threads} -o {output.bam} --write-index -

        echo -e "[INFO] minimap2 is done!"; }} \\
        1> {log} 2>&1
        """
