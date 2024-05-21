rule call_sv_svision:
    input:
        bam=ancient("minimap2/{sample}/{sample}.sorted.bam"),
        fasta=config["fasta"],
    output:
        vcf=protected("svision/{sample}/svision.vcf"),
        tab_rename=temp("svision/{sample}/rename.tab"),
    params:
        dir_out=directory("svision/{sample}"),
        model_svision=config["model_svision"],
        min_num_reads=config["min_num_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_length_sv=config["min_length_sv"],
    threads: config["threads"]
    log:
        "logs/{sample}/call_sv_svision.log",
    shell:
        """
        {{ echo -e "{wildcards.sample}.chr1\\t{wildcards.sample}" > {output.tab_rename}
        for chr in chr{{1..22}} chrX chrY; do
            echo -e "[INFO] Running SVision on chromosome ${{chr}}..."
            SVision \\
            -t {threads} \\
            -s {params.min_num_reads} \\
            --min_mapq {params.min_quality_mapping} \\
            --min_sv_size {params.min_length_sv} \\
            --max_sv_size 999999999999 \\
            --qname \\
            --graph \\
            --min_gt_depth {params.min_num_reads} \\
            -o {params.dir_out} \\
            -b {input.bam} \\
            -m {params.model_svision} \\
            -g {input.fasta} \\
            -n {wildcards.sample}.${{chr}} \\
            -c ${{chr}}
            sleep 10
        done

        vcf_svision_chr1={params.dir_out}/{wildcards.sample}.chr1.svision.s{params.min_num_reads}.graph.vcf
        {{ grep '^#' ${{vcf_svision_chr1}}; cat {params.dir_out}/{wildcards.sample}.*.vcf | grep -v '^#'; }} \\
            | awk '/^##INFO=<ID=GFA_L/ && !f {{print "##INFO=<ID=GFA_ID,Number=.,Type=String,Description=\\"GFA_ID\\">"; f=1}} 1' \\
            | awk 'BEGIN{{FS=OFS="\\t"}} /^#/ || $5 == "<CSV>" {{print; next}} {{split($8, a, ";"); for(i in a) {{if(a[i] ~ /^SVTYPE=/) {{split(a[i], b, "="); if(b[2] == "tDUP") b[2] = "DUP:TANDEM"; gsub("<SV>", "<"b[2]">", $5)}}}}}}1' \\
            | awk '/^##ALT/ && !f {{print "##ALT=<ID=INS,Description=\\"INS\\">\\n##ALT=<ID=INV,Description=\\"INV\\">\\n##ALT=<ID=DUP,Description=\\"DUP\\">\\n##ALT=<ID=DUP:TANDEM,Description=\\"DUP:TANDEM\\">\\n##ALT=<ID=DEL,Description=\\"DEL\\">"; f=1}} 1' \\
            | awk '/^##INFO=<ID=SUPPORT,/ {{sub("Type=String", "Type=Integer")}} 1' \\
            | bcftools reheader -s {output.tab_rename} \\
            > {output.vcf}

        mkdir {params.dir_out}/chrs
        mv {params.dir_out}/*chr*.vcf {params.dir_out}/*chr*.txt {params.dir_out}/*.log {params.dir_out}/chrs/

        echo -e "[INFO] SVision is done!"; }} \\
        1> {log} 2>&1
        """