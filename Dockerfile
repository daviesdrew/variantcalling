FROM nfcore/base:1.9
LABEL authors="Drew Davies" \
      description="Docker image for running the happymappy pipeline"

COPY environment.yml /
RUN cat environment.yml
RUN ls . -la
RUN conda update conda && \
    conda env create -f /environment.yml && \
    conda clean -a

ENV PATH /opt/conda/envs/variantcalling-1.0.0/bin:$PATH
RUN conda info --envs
