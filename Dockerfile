FROM nfcore/base:1.7
LABEL authors="Phoenix Logan Lucy Li" \
      description="Docker image containing all requirements for nf-core/epigen pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-epigen-1.0dev/bin:$PATH

FROM ubuntu:18.10

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates software-properties-common gnupg2 gnupg1 \
      && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
      && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' \
      && apt-get install r-base 

#RUN R -e "devtools::install_github('lucymli/EpiGenR')"

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
