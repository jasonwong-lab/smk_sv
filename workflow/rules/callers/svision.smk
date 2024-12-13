rule call_sv_svision:
    container:
        "docker://mhjiang97/smk_sv:latest"
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
    threads: math.floor(workflow.cores / NUMBER_SAMPLES / NUMBER_CALLERS)
    log:
        "logs/{sample}/call_sv_svision.log",
    shell:
        """
        {{ for chr in chr{{1..22}} chrX chrY; do
            lock_file="{params.dir_out}/${{chr}}.lock"
            if [ ! -f ${{lock_file}} ]; then
                echo -e "[INFO] Running SVision on chromosome ${{chr}}..."

                n_graphs=$(find {params.dir_out}/{wildcards.sample}/graphs/${{chr}}-* -type d | wc -l)
                if [ ${{n_graphs}} -gt 0]; then
                    rm -rf {params.dir_out}/{wildcards.sample}/graphs/${{chr}}-*/
                fi
                n_segments=$(find {params.dir_out}/{wildcards.sample}/segments/${{chr}}.*.bed | wc -l)
                if [ ${{n_segments}} -gt 0]; then
                    rm -rf {params.dir_out}/{wildcards.sample}/segments/${{chr}}.*.bed
                fi
                n_predict_results=$(find {params.dir_out}/{wildcards.sample}/predict_results/${{chr}}.predict.* | wc -l)
                if [ ${{n_predict_results}} -gt 0]; then
                    rm -rf {params.dir_out}/{wildcards.sample}/predict_results/${{chr}}.predict.*
                fi

                SVision -t {threads} -s {params.min_num_reads} --min_mapq {params.min_quality_mapping} --min_sv_size {params.min_length_sv} --max_sv_size 999999999999 --qname --graph --min_gt_depth {params.min_num_reads} -o {params.dir_out} -b {input.bam} -m {params.model_svision} -g {input.fasta} -n {wildcards.sample}.${{chr}} -c ${{chr}}

                sleep 10
                touch ${{lock_file}}
                sleep 10
            else
                echo -e "[INFO] SVision already run for chromosome ${{chr}}, skipping..."
            fi
        done

        vcf_svision_lead=$(find {params.dir_out}/{wildcards.sample}.chr*.svision.s{params.min_num_reads}.graph.vcf | head -1)
        chrom=$(basename $vcf_svision_lead | cut -d'.' -f2)
        echo -e "{wildcards.sample}.${{chrom}}\\t{wildcards.sample}" > {output.tab_rename}
        {{ grep '^#' ${{vcf_svision_lead}}; cat {params.dir_out}/{wildcards.sample}.chr*.svision.s{params.min_num_reads}.graph.vcf | grep -v '^#'; }} \\
            | awk '/^##INFO=<ID=GFA_L/ && !f {{print "##INFO=<ID=GFA_ID,Number=.,Type=String,Description=\\"GFA_ID\\">"; f=1}} 1' \\
            | awk 'BEGIN{{FS=OFS="\\t"}} /^#/ || $5 == "<CSV>" {{print; next}} {{split($8, a, ";"); for(i in a) {{if(a[i] ~ /^SVTYPE=/) {{split(a[i], b, "="); if(b[2] == "tDUP") b[2] = "DUP:TANDEM"; gsub("<SV>", "<"b[2]">", $5)}}}}}}1' \\
            | awk '/^##ALT/ && !f {{print "##ALT=<ID=INS,Description=\\"INS\\">\\n##ALT=<ID=INV,Description=\\"INV\\">\\n##ALT=<ID=DUP,Description=\\"DUP\\">\\n##ALT=<ID=DUP:TANDEM,Description=\\"DUP:TANDEM\\">\\n##ALT=<ID=DEL,Description=\\"DEL\\">"; f=1}} 1' \\
            | awk '/^##INFO=<ID=SUPPORT,/ {{sub("Type=String", "Type=Integer")}} 1' \\
            | bcftools reheader -s {output.tab_rename} \\
            | awk -F'\\t' -v OFS='\\t' '{{if ($0 ~ /^#/) {{print $0;}} else {{$3=$1"_"$2"_"$3; print $0}}}}' \\
        > {output.vcf}

        mkdir {params.dir_out}/chrs
        mv {params.dir_out}/*chr*.vcf {params.dir_out}/*chr*.txt {params.dir_out}/*.log {params.dir_out}/chrs/

        echo -e "[INFO] SVision is done!"; }} \\
        1> {log} 2>&1
        """
