FROM nfcore/base:1.9
LABEL authors="Drew Davies" \
      description="Docker image containing all software requirements for the nf-core/illuminavariantcalling pipeline"

COPY environment.yml  /
RUN cat environment.yml
COPY ./test/CIN-2-QCM18-1264_S1_L001_R1_001.fastq.gz /
COPY ./test/CIN-2-QCM18-1264_S1_L001_R2_001.fastq.gz /
RUN ls . -l
COPY ./test/ref.fa / 
COPY ./test/phix.fa / 
RUN ls . -l
RUN conda update conda && \
    conda env create -f /environment.yml && \
    conda clean -a

ENV PATH /opt/conda/envs/variantcalling-1.0.0/bin:$PATH
RUN conda info --envs
