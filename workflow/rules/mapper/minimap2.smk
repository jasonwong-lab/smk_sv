rule minimap2_index:
    conda:
        "../../envs/minimap2.yaml"
    input:
        fasta=config["fasta"],
    output:
        index_minimap2=protected(config["index_minimap2"]),
    threads: 1
    log:
        "logs/minimap2_index.log",
    shell:
        "minimap2 -t {threads} -d {output.index_minimap2} {input.fasta}"


rule minimap2:
    conda:
        "../../envs/minimap2.yaml"
    input:
        fq=f"{config['dir_data']}/{{sample}}{config['suffix_fastq']}",
        index_minimap2=ancient(config["index_minimap2"]),
    output:
        bam=protected("minimap2/{sample}/{sample}.sorted.bam"),
        bai=protected("minimap2/{sample}/{sample}.sorted.bam.csi"),
    params:
        read_group="@RG\\tID:{sample}\\tSM:{sample}\\tLB:DNA\\tPL:ONT",
    threads: 1
    log:
        "logs/{sample}/minimap2.log",
    shell:
        """
        {{ minimap2 \\
            -ax map-ont \\
            -Y \\
            -t {threads} \\
            -R "{params.read_group}" \\
            {input.index_minimap2} {input.fq} \\
            | samtools sort \\
                -@ {threads} \\
                -o {output.bam} \\
                --write-index \\
                -; }} \\
        1> {log} 2>&1
        """
