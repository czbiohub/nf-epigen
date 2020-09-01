FROM nfcore/base:1.10.2
LABEL authors="Phoenix Logan Lucy Li" \
      description="Docker image containing all requirements for nf-core/epigen pipeline"


RUN apt-get update && apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates software-properties-common gnupg2 gnupg1 dirmngr apt-transport-https \
      && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' \
      && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian buster-cran40/' \
      && apt update \
      && apt-get -y install r-base

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-epigen-1.0dev/bin:$PATH


RUN R -e "options(needs.promptUser = FALSE); devtools::install_github('lucymli/EpiGenR')"

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
