checkpoint svision:
    container:
        "docker://jiadongxjtu/svision:latest"
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        fasta=config["fasta"],
    output:
        touch(vcfs_svision),
        graphs=directory("svision/{sample}/chroms/graphs"),
        segments=temp(touch(directory("svision/{sample}/chroms/segments"))),
        predict_results=temp(
            touch(directory("svision/{sample}/chroms/predict_results"))
        ),
    params:
        dir="svision/{sample}/chroms",
        chroms=CHROMS,
        model=config["model_svision"],
        min_reads=config["min_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_size=config["min_size"],
        max_size=config["max_size"],
    threads: 1
    log:
        "logs/{sample}/svision.log",
    shell:
        """
        {{ for chrom in {params.chroms}; do
            lock_file="{params.dir}/${{chrom}}.lock"
            if [ -f ${{lock_file}} ]; then
                echo -e "[INFO] SVision already run for chromosome ${{chrom}}, skipping..."
                continue
            fi

            echo -e "[INFO] Running SVision on chromosome ${{chrom}}..."

            SVision \\
                -t {threads} \\
                -s {params.min_reads} \\
                --min_mapq {params.min_quality_mapping} \\
                --min_sv_size {params.min_size} \\
                --max_sv_size {params.max_size} \\
                --qname \\
                --graph \\
                --min_gt_depth {params.min_reads} \\
                -o {params.dir} \\
                -b {input.bam} \\
                -m {params.model} \\
                -g {input.fasta} \\
                -n {wildcards.sample}.${{chrom}} \\
                -c ${{chrom}}

            sleep 10
            touch ${{lock_file}}
            sleep 10
        done; }} \\
        1> {log} 2>&1
        """


rule format_svision:
    conda:
        "../../envs/bcftools.yaml"
    input:
        vcfs_svision,
    output:
        tab=temp("svision/{sample}/rename.tab"),
        vcf=protected("svision/{sample}/svision.vcf"),
    params:
        get_format_svision_parameters,
    log:
        "logs/{sample}/format_svision.log",
    shell:
        """
        {{ echo -e "{wildcards.sample}.{params[0][chrom_lead]}\\t{wildcards.sample}" > {output.tab}

        {{ grep '^#' {params[0][vcf_lead]}; cat {input} | grep -v '^#'; }} \\
            | awk '/^##INFO=<ID=GFA_L/ && !f {{print "##INFO=<ID=GFA_ID,Number=.,Type=String,Description=\\"GFA_ID\\">"; f=1}} 1' \\
            | awk 'BEGIN{{FS=OFS="\\t"}} /^#/ || $5 == "<CSV>" {{print; next}} {{split($8, a, ";"); for(i in a) {{if(a[i] ~ /^SVTYPE=/) {{split(a[i], b, "="); if(b[2] == "tDUP") b[2] = "DUP:TANDEM"; gsub("<SV>", "<"b[2]">", $5)}}}}}}1' \\
            | awk '/^##ALT/ && !f {{print "##ALT=<ID=INS,Description=\\"INS\\">\\n##ALT=<ID=INV,Description=\\"INV\\">\\n##ALT=<ID=DUP,Description=\\"DUP\\">\\n##ALT=<ID=DUP:TANDEM,Description=\\"DUP:TANDEM\\">\\n##ALT=<ID=DEL,Description=\\"DEL\\">"; f=1}} 1' \\
            | awk '/^##INFO=<ID=SUPPORT,/ {{sub("Type=String", "Type=Integer")}} 1' \\
            | bcftools reheader -s {output.tab} \\
            | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
            > {output.vcf}; }} \\
        1> {log} 2>&1
        """
