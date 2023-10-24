FROM rocker/tidyverse@sha256:0c62e932eeb3c5dd2529ed73274e670ffae169d7368073a59bc2282f28e99a6f
LABEL XenofilteR based on rocker/tidyverse 4.2.3

WORKDIR .

RUN mkdir ../opt/lib

RUN apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y \
 git \
 r-cran-futile.logger \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*


#packages from Bioconductor
COPY install_BioconductorLibs.R ../opt/lib
RUN Rscript ../opt/lib/install_BioconductorLibs.R

RUN R -e "install.packages('devtools', repos='http://cloud.r-project.org')"

RUN installGithub.r NKI-GCF/XenofilteR

ENV SHELL=bash

WORKDIR /usr/src/workdir






