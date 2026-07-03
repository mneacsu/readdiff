FROM python:3.12.7

RUN apt-get update && apt-get install -y \
    locales \
    perl \
    wget \
    zip \
    unzip \
    gzip \
    make \
    cmake \
    gcc \
    g++ \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libssl-dev \
    git \
    ca-certificates \
    bedtools \
    samtools \
    rna-star \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --trusted-host pypi.org --trusted-host files.pythonhosted.org \
    numpy \
    pandas \
    scipy \
    tqdm \
    scanpy \
    anndata \
    pydeseq2 \
    matplotlib \
    seaborn \
    umap-learn \
    tqdm \
    joypy \
    plotly \
    dash \
    dash-mantine-components \
    pyyaml \
    biopython \
    typer \
    snakemake \
    snakemake-executor-plugin-slurm

# HOMER installation
RUN mkdir -p /opt/homer && \
    cd /opt/homer && \
    wget http://homer.ucsd.edu/homer/configureHomer.pl && \
    chmod +x configureHomer.pl && \
    perl configureHomer.pl -install homer && \
    perl configureHomer.pl -install hg38 && \
    perl configureHomer.pl -install mm10

ENV HOMER_HOME=/opt/homer
ENV HOMER_DATA_DIR=/opt/homer/data
ENV PATH=/opt/homer/bin:$PATH

# NEEDLE installation
RUN git config --global http.sslVerify false && cd /opt && \
    git clone https://github.com/seqan/needle.git && \
    cd needle && git checkout needle-v1.0.3 && cd .. && \
    mkdir build-needle && cd build-needle && \
    cmake ../needle -DCMAKE_BUILD_TYPE=Release && \
    make
ENV PATH=/opt/build-needle/bin:$PATH

# readdiff installation
COPY . /opt/readdiff
RUN pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org --no-cache-dir /opt/readdiff

CMD ["/bin/bash"]

