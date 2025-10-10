# Dockerfile for VNTR_Pipeline

ARG image_name=ubuntu:20.04
FROM ${image_name}

LABEL maintainer="David Horner <david.horner@meduniwien.ac.at>"

SHELL ["/bin/bash", "-c"]
ENV PATH=/opt/bin:/opt/conda/bin:$PATH
ENV CONDA_AUTO_UPDATE_CONDA=false


# --- System dependencies ---
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        zlib1g-dev libbz2-dev liblzma-dev libboost-graph-dev \
        libcurl4-openssl-dev libssl-dev default-jdk pkg-config \
	cpanminus build-essential perl \
        xz-utils curl jq wget bzip2 make g++ git parallel && \
    rm -rf /var/lib/apt/lists/* && \
	cpanm Text::Fuzzy

# --- Miniconda installation ---
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    conda init bash

# --- Configure Conda ---
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda install -n base -c defaults "conda>=23.9.0" && \
    conda install -n base conda-libmamba-solver && \
    conda config --set solver libmamba && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# --- Create and configure 'clair3-env' + compile VNTR_Pipeline ---
WORKDIR /opt/bin

COPY /clair3 /opt/bin

RUN conda create -n clair3-env python=3.9 -y && \
    conda run -n clair3-env conda install -y -c conda-forge \
        pypy3.6 tensorflow-cpu=2.8.0 pytables \
        xz zlib bzip2 automake curl && \
    conda run -n clair3-env conda install -y -c bioconda \
        samtools=1.15.1 whatshap && \
    conda run -n clair3-env conda install -y -c anaconda \
        pigz cffi=1.14.4 && \
    conda run -n clair3-env pypy3 -m ensurepip && \
    conda run -n clair3-env pypy3 -m pip install mpmath==1.2.1 && \
    conda run -n clair3-env pip install tensorflow-addons && \
    conda clean -afy && \
    # Compile realign tools
    cd /opt/bin/preprocess/realign && \
    g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp && \
    g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp && \
    # Download clair3 models
    mkdir -p /opt/models && \
    wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz -P /opt/models && \
    tar -zxvf /opt/models/clair3_models.tar.gz -C /opt/models && \
    rm /opt/models/clair3_models.tar.gz && \
    # Compile full pipeline
    make -C /opt/bin PREFIX=/opt/conda/envs/clair3-env PYTHON=/opt/conda/envs/clair3-env/bin/python && \
    rm -rf /opt/bin/samtools-* /opt/bin/longphase-*

# --- Download rerio models ---
RUN git clone https://github.com/nanoporetech/rerio /opt/rerio && \
    python3 /opt/rerio/download_model.py --clair3

# --- Additional environments ---

WORKDIR /

# Python + MAFFT + R
COPY requirements.txt /
RUN conda create -n python-env python=3.8 -y && \
    conda run -n python-env pip install --no-cache-dir -r /requirements.txt && \
    conda run -n python-env pip install trviz && \
    conda run -n python-env pip install pandas && \
    conda run -n python-env pip install openpyxl && \
    conda run -n python-env conda install -y -c bioconda mafft && \
    conda run -n python-env conda install -y -c conda-forge \
        r-essentials r-base r-stringi r-openxlsx r-plotrix zlib && \
    conda run -n python-env conda install -y -c bioconda \
        bioconductor-biostrings bioconductor-shortread && \
    rm /requirements.txt

# Tool-specific envs
RUN conda create -n minimap-env minimap2 -y && \
    conda create -n samtools-env samtools -y && \
    conda create -n whatshap-env whatshap -y && \
    conda create -n bcftools-env bcftools -y

# Sniffles (latest available)
RUN version=$(conda search sniffles --full-name | awk 'END {print $2}') && \
    conda create -n sniffles-env sniffles="$version" -y

# --- Install Canu ---
RUN curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz -o canu.tar.xz && \
    tar -xJf canu.tar.xz && \
    rm canu.tar.xz && \
    git clone https://github.com/marbl/canu.git && \
    cd /canu/src && \
    make && \
    ln -s /canu/build/bin/canu /usr/bin/canu

# --- Copy external resources ---
COPY /test_data /test_data
COPY /references /references
COPY /scripts /scripts

RUN ln -s /scripts/bin/VNTR_pipeline.sh /usr/bin/VNTR_pipeline

# --- Final dummy init (optional) ---
RUN conda init && \
    bash

