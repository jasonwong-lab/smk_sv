rule make_index_minimap2:
    input:
        fasta=config["fasta"],
    output:
        index_minimap2=config["index_minimap2"],
    threads: config["threads"]
    shell:
        "minimap2 -t {threads} -d {output.index_minimap2} {input.fasta}"


rule map_minimap2:
    input:
        index_minimap2=ancient(rules.make_index_minimap2.output.index_minimap2),
    output:
        bam=protected("minimap2/{sample}/{sample}.sorted.bam"),
        bai=protected("minimap2/{sample}/{sample}.sorted.bam.bai"),
    params:
        dir_data=config["dir_data"],
        suffix_fastq=config["suffix_fastq"],
    threads: config["threads"]
    log:
        "logs/{sample}/map_minimap2.log",
    shell:
        """
        {{ minimap2 \\
        -ax map-ont -Y -t {threads} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:DNA\\tPL:ONT\" \\
        {input.index_minimap2} {params.dir_data}/{wildcards.sample}.{params.suffix_fastq} \\
            | samtools view -Sb -@ {threads} - \\
            | samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}; }} \\
        1> {log} 2>&1
        """
