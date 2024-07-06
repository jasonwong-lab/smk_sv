FROM jiadongxjtu/svision:latest

ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
ENV PATH /opt/miniconda3/bin:/opt/bin:/opt/samtools/bin:${PATH}:/opt/snpEff/exec:/bin:${PATH}
ENV BASH_ENV ~/.bashrc
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/lib

RUN rm /'=3.6' && \
    apt-get update --fix-missing && \
    yes | apt-get upgrade && \
    apt-get install libcurl4 openssl libssl-dev locales unzip wget bzip2 make g++ libboost-graph-dev libcurl4-openssl-dev curl tabix openjdk-17-jdk -y && \
    locale-gen en_US.UTF-8 && \
    ln -s /usr/lib/x86_64-linux-gnu/libcrypto.so /usr/lib/libcrypto.so.10 && \
    ln -s /usr/lib/x86_64-linux-gnu/libcurl.so.4 /usr/lib/libcurl.so.4 && \
    ldconfig && \
    apt install --no-install-recommends software-properties-common dirmngr -y && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh && \
    bash /opt/miniconda.sh -b -u -p /opt/miniconda3 && \
    rm -rf /opt/miniconda.sh && \
    /opt/miniconda3/bin/conda init bash && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git /tmp/htslib && \
    git clone https://github.com/samtools/bcftools.git /tmp/bcftools && \
    cd /tmp/bcftools && \
    make && \
    mv /tmp/bcftools/bcftools /bin/bcftools && \
    curl -o /opt/snpEff_latest_core.zip 'https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip' && \
    cd /opt && \
    unzip /opt/snpEff_latest_core.zip && \
    rm /opt/snpEff_latest_core.zip && \
    conda create -n vcf2maf -c bioconda vcf2maf==1.6.21 -y && \
    cd /tmp && \
    wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xvf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure --prefix=/opt/samtools && \
    make && \
    make install && \
    wget -O /bin/duphold https://github.com/brentp/duphold/releases/download/v0.2.3/duphold && \
    chmod +x /bin/duphold && \
    mkdir -p /opt/bin && \
    cd /opt/bin && \
    ln -s /opt/miniconda3/envs/vcf2maf/bin/vcf2maf.pl ./ && \
    ln -s /opt/miniconda3/envs/vcf2maf/bin/perl ./ && \
    sed -i '410s/^/# /' /opt/miniconda3/envs/vcf2maf/bin/vcf2maf.pl && \
    sed -i '1s/.*/\#\!\/opt\/miniconda3\/envs\/vcf2maf\/bin\/perl/' /opt/miniconda3/envs/vcf2maf/bin/vcf2maf.pl && \
    conda clean -a && \
    apt-get clean && \
    rm -rf /bar/lib/apt/lists/* && \
    rm -rf /tmp/* && \
    rm -rf /opt/miniconda3/pkgs/* && \
    rm -rf /root/.cache/pip
