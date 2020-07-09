FROM nfcore/base:1.7
LABEL authors="Phoenix Logan Lucy Li" \
      description="Docker image containing all requirements for nf-core/epigen pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-epigen-1.0dev/bin:$PATH

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron