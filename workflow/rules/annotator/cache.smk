rule download_snpeff_cache:
    conda:
        "../../envs/awscli.yaml"
    output:
        dir=directory(path_cache_snpeff),
    params:
        dir_src=f"{config['genome']}.{config['version_snpeff']}",
        dir_trgt=config["cache_snpeff"],
    log:
        "logs/download_snpeff_cache.log",
    shell:
        """
        {{ url_base="s3://annotation-cache/snpeff_cache"
        caches_aws=( $(aws s3 ls --no-sign-request ${{url_base}}/ | grep PRE | sed 's/^[[:space:]]*//' | cut -d ' ' -f 2 | sed 's/\\/$//' | tr '\n' ' ') )

        if [[ " ${{caches_aws[@]}} " =~ " {params.dir_src} " ]]; then
            aws s3 --no-sign-request sync ${{url_base}}/{params.dir_src} {params.dir_trgt}
        else
            snpEff download {params.dir_src} -dataDir {params.dir_trgt}
        fi; }} \\
        1> {log} 2>&1
        """


rule download_vep_cache:
    conda:
        "../../envs/awscli.yaml"
    output:
        dir=directory(path_cache_vep),
    params:
        dir_src=f"{config['version_vep']}_{config['genome']}",
        dir_trgt=f"{config['cache_vep']}/{config['species']}",
        cache=config["cache_vep"],
        species=config["species"],
        version=config["version_vep"],
        genome=config["genome"],
    log:
        "logs/download_vep_cache.log",
    shell:
        """
        {{ url_base="s3://annotation-cache/vep_cache"
        caches_aws=( $(aws s3 ls --no-sign-request ${{url_base}}/ | grep PRE | sed 's/^[[:space:]]*//' | cut -d ' ' -f 2 | sed 's/\\/$//' | tr '\n' ' ') )

        if [[ " ${{caches_aws[@]}} " =~ " {params.dir_src} " ]]; then
            aws s3 --no-sign-request sync ${{url_base}}/{params.dir_src} {params.dir_trgt}
        else
            vep_install \\
            --CACHEDIR {params.cache} \\
            --DESTDIR {params.cache} \\
            --PLUGINSDIR {params.cache}/Plugins/ \\
            --CACHE_VERSION {params.version} \\
            --SPECIES {params.species} \\
            --ASSEMBLY {params.genome} \\
            --PREFER_BIN --NO_UPDATE --AUTO cf
        fi; }} \\
        > {log} 2>&1
        """


rule download_annotsv_cache:
    conda:
        "../../envs/git.yaml"
    output:
        **get_annotsv_cache_outputs(),
        annotsv=temp(directory(f"{config['cache_annotsv']}/AnnotSV")),
    params:
        **get_annotsv_cache_parameters(),
        dir=config["cache_annotsv"],
        version=config["version_annotsv"],
    log:
        "logs/download_annotsv_cache.log",
    shell:
        """
        {{ git clone https://github.com/lgmgeo/AnnotSV.git {output.annotsv}
        cd {output.annotsv} || exit 1
        git checkout {params.version}

        make PREFIX=. install
        make PREFIX=. {params.arg_install}

        for dir in {params.dirs}; do
            mv ${{dir}} ../
        done; }} \\
        > {log} 2>&1
        """
