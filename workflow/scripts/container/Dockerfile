FROM jiadongxjtu/svision:latest

COPY delly /bin/delly
COPY duphold /bin/duphold

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/miniconda3/bin:/opt/ensembl-vep:/opt/snpEff/exec:/bin:/opt/bin:${PATH}
ENV BASH_ENV ~/.bashrc

RUN rm /'=3.6' && \
    apt-get update --fix-missing && \
    yes | apt-get upgrade && \
    apt-get install wget bzip2 make g++ libboost-graph-dev libcurl4-openssl-dev curl tabix perl libarchive-zip-perl libdbi-perl bioperl libmodule-build-perl libdbd-mysql-perl libset-intervaltree-perl libjson-perl libperlio-gzip-perl libtest-warnings-perl default-jdk -y && \
    apt install --no-install-recommends software-properties-common dirmngr -y && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/" && \
    apt install --no-install-recommends r-base -y && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh && \
    bash /opt/miniconda.sh -b -u -p /opt/miniconda3 && \
    rm -rf /opt/miniconda.sh && \
    /opt/miniconda3/bin/conda init bash && \
    . /root/.bashrc && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n clair3 python=3.9.0 -y && \
    conda install -n clair3 -c conda-forge pypy3.6 -y && \
    conda activate clair3 && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install mpmath==1.2.1 && \
    conda install -c conda-forge tensorflow-cpu==2.8.0 -y && \
    conda install -c conda-forge pytables -y && \
    pip install tensorflow-addons && \
    conda install -c anaconda pigz -y && \
    conda install -c anaconda cffi=1.14.4 -y && \
    conda install -c conda-forge parallel=20191122 zstd -y && \
    conda install -c conda-forge -c bioconda samtools=1.15.1 -y && \
    conda install -c conda-forge -c bioconda whatshap=1.7 -y && \
    conda install -c conda-forge xz zlib bzip2 -y && \
    conda install -c conda-forge automake -y && \
    conda install -c bioconda clair3 -y && \
    conda create -n severus severus -y && \
    conda create -n debreak debreak -y && \
    conda create -n nanosv nanosv -y && \
    conda create -n nanovar -c bioconda python=3.11 samtools bedtools minimap2 -y && \
    conda activate nanovar && \
    pip install nanovar && \
    pip install cloud-tpu-client && \
    cd /tmp && \
    curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf - && \
    mv /tmp/minimap2-2.28_x64-linux/minimap2 /bin/minimap2 && \
    conda create -n cutesv -c bioconda cutesv && \
    conda create -n sniffles sniffles=2.2 && \
    conda create -n svim -c bioconda svim && \
    echo "conda activate base" >> ~/.bashrc && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git /tmp/htslib && \
    git clone https://github.com/samtools/bcftools.git /tmp/bcftools && \
    cd /tmp/bcftools && \
    make && \
    mv /tmp/bcftools/bcftools /bin/bcftools && \
    git clone https://github.com/fritzsedlazeck/SURVIVOR.git /tmp/SURVIVOR && \
    cd /tmp/SURVIVOR/Debug && \
    make && \
    mv /tmp/SURVIVOR/Debug/SURVIVOR /bin/SURVIVOR && \
    curl -o /opt/snpEff_latest_core.zip 'https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip' && \
    cd /opt && \
    unzip /opt/snpEff_latest_core.zip && \
    rm /opt/snpEff_latest_core.zip && \
    git clone https://github.com/Ensembl/ensembl-vep.git /opt/ensembl-vep && \
    cd /opt/ensembl-vep && \
    /usr/bin/perl INSTALL.pl && \
    conda create -n vcf2maf -c bioconda vcf2maf && \
    conda clean -a && \
    apt-get clean && \
    rm -rf mskcc-vcf2maf.tar.gz && \
    rm -rf /bar/lib/apt/lists/* && \
    rm -rf /tmp/* && \
    rm -rf /opt/miniconda3/pkgs/* && \
    rm -rf /root/.cache/pip

RUN echo "conda activate clair3" >> ~/.bashrc && \
    mkdir /opt/bin && \
    ln -s /opt/miniconda3/envs/cutesv/bin/cuteSV /opt/bin/cuteSV && \
    ln -s /opt/miniconda3/envs/sniffles/bin/sniffles /opt/bin/sniffles && \
    ln -s /opt/miniconda3/envs/svim/bin/svim /opt/bin/svim && \
    cd /opt/bin && \
    ln -s /opt/miniconda3/envs/debreak/bin/debreak* ./ && \
    ln -s /opt/miniconda3/envs/nanosv/bin/NanoSV /opt/bin/NanoSV && \
    ln -s /opt/miniconda3/envs/nanovar/bin/nanovar /opt/bin/nanovar && \
    ln -s /opt/miniconda3/envs/clair3/bin/run_clair3.sh /opt/bin/run_clair3.sh && \
    ln -s /opt/miniconda3/envs/clair3/bin/scripts /opt/bin/scripts && \
    ln -s /opt/miniconda3/envs/severus/bin/severus /opt/bin/severus && \
    ln -s /opt/miniconda3/envs/vcf2maf/bin/vcf2maf.pl /opt/bin/vcf2maf.pl && \
    ln -s /opt/miniconda3/envs/clair3/bin/samtools /opt/bin/samtools && \
    ln -s /opt/miniconda3/envs/clair3/bin/clair3.py /opt/bin/clair3.py && \
    sed -i '1s/.*/\#\!\/opt\/miniconda3\/envs\/debreak\/bin\/python3/' /opt/bin/debreak && \
    sed -i '3c\SCRIPT_PATH=$(dirname $(/bin/readlink -f "$0"))' /opt/bin/run_clair3.sh && \
    sed -i '410s/^/# /' /opt/miniconda3/envs/vcf2maf/bin/vcf2maf.pl && \
    sed -i '1s/.*/\#\!\/opt\/miniconda3\/envs\/vcf2maf\/bin\/perl/' /opt/miniconda3/envs/vcf2maf/bin/vcf2maf.pl

RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt install -y g++-11 && \
    apt-get clean && \
    rm -rf /bar/lib/apt/lists/* && \
    rm -rf /tmp/*

ENV PATH ${PATH}:/opt/miniconda3/envs/clair3/bin
