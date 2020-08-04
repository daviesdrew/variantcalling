FROM nfcore/base:1.9
LABEL authors="Drew Davies" \
      description="Docker image for running the happymappy pipeline"

COPY environment.yml /
RUN conda update conda && \
    conda env create -f /environment.yml && \

ENV PATH /opt/conda/envs/variantcalling-1.0.0/bin:$PATH
