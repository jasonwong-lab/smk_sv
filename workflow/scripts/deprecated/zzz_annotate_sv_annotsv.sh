#! /bin/bash

# """
# sleep $((RANDOM % 10 + 1))

# lockfile="annotsv.{wildcards.sample}.{wildcards.caller}.lock"
# while [ -f ${{lockfile}} ]; do
#     sleep 10
# done
# touch ${{lockfile}}

# input_annotsv={input.vcf}
# size_max=300
# size_input=$(du -m {input.vcf} | cut -f1)
# if {{ [ {wildcards.caller} == "cutesv" ] && [ {wildcards.type_sv} == "DEL" ]; }} \\
#     || {{ [ {wildcards.caller} == "debreak" ] && [ {wildcards.type_sv} == "DEL" ]; }} \\
#     && [ "${{size_input}}" -ge ${{size_max}} ]; then
#     awk 'BEGIN {{OFS=FS="\\t"}} !/^#/ {{$4 = "N"; $5 = "<DEL>"}} 1' "${{input_annotsv}}" > "${{input_annotsv%.*}}".4_5.vcf
#     input_annotsv="${{input_annotsv%.*}}".4_5.vcf
# fi
# AnnotSV -genomeBuild {params.genome} -SVinputFile ${{input_annotsv}} -outputFile {output.tsv} -SVminSize 1 -overwrite 1 1> {log} 2>&1

# sleep 5
# rm ${{lockfile}}
# """