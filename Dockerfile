FROM nfcore/base:1.9
LABEL authors="Drew Davies" \
      description="Docker image containing all software requirements for the nf-core/illuminavariantcalling pipeline"

# Install the conda environment
COPY environment.yml  /
RUN cat environment.yml
COPY CIN-2-QCM18-1264_S1_L001_R1_001.fastq.gz /
COPY CIN-2-QCM18-1264_S1_L001_R2_001.fastq.gz /
COPY ref.fa / 
COPY phix.fa / 
RUN conda update conda && \
    conda env create -f /environment.yml && \
    conda clean -a
RUN ls
#RUN wget -q https://get.nextflow.io | bash

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/variantcalling-1.0.0/bin:$PATH
RUN conda info --envs
#RUN /bin/bash -c "source activate variantcalling"

# Dump the details of the installed packages to a file for posterity
#RUN conda env export --name variantcalling > variantcalling.yml

# Dump the details of the installed packages to a file for posterity
#RUN conda env export --name nf-core-illuminavariantcalling-1.0dev > nf-core-illuminavariantcalling-1.0dev.yml
