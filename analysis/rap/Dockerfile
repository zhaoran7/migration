FROM rocker/r-ver:4.4.2
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libharfbuzz-dev \
    libfribidi-dev libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
    libgit2-dev libglpk-dev libgmp3-dev libmpfr-dev zlib1g-dev \
    gawk gzip unzip zip && rm -rf /var/lib/apt/lists/*
RUN R -q -e "install.packages(c('data.table','haven','tidyselect','survival', \
    'coxme','lme4','RANN','sf','openxlsx','jsonlite','xml2','zip','topr','cluster','ragg'), \
    repos='https://cloud.r-project.org')"
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 LANG=C.UTF-8 LC_ALL=C.UTF-8
WORKDIR /tmp
