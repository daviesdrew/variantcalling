FROM nfcore/base:1.9
LABEL authors="Drew Davies" \
      description="Docker image containing all software requirements for the nf-core/illuminavariantcalling pipeline"

# Install the conda environment
COPY environment.yml  /
RUN conda update conda && \
    conda env create -f /environment.yml && \
    conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/variantcalling/bin:$PATH


# Dump the details of the installed packages to a file for posterity
#RUN conda env export --name variantcalling > variantcalling.yml

# Dump the details of the installed packages to a file for posterity
#RUN conda env export --name nf-core-illuminavariantcalling-1.0dev > nf-core-illuminavariantcalling-1.0dev.yml
