rule faidx:
    input:
        fasta=config["fasta"],
    output:
        fai=f"{config['fasta']}.fai",
    log:
        "logs/faidx.log",
    shell:
        """
        samtools faidx {input.fasta} 1> {log} 2>&1
        """
