FROM rocker/r-ver:4.2.1 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    g++ \
    git \
    nano \
    wget \
    zlib1g-dev \
    fftw3-dev \
    libgit2-dev \
    libssl-dev \
    libssh2-1-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libomp-dev \
    cmake \
    libxt6 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /Programme

ENV _R_SHLIB_STRIP_=true

RUN git clone --recursive https://github.com/microsoft/LightGBM \
 && cd LightGBM \
 && mkdir build \
 && cd build \
 && cmake .. \
 && make -j4

RUN install2.r --error --skipinstalled \
    shiny \
    git2r \
    data.table \
    plotrix \
    geosphere \
    DT \
    shinyFiles \
    shinydashboard \
    lubridate \
    wkb \
    shinybusy \
    shinythemes \
    shinyjs \
    lubridate \
    dplyr \
    aws.s3

COPY Rprofile.site /usr/local/lib/R/etc/

# Luna
RUN cp LightGBM/lib_lightgbm.so /usr/local/lib/ \
 && cp LightGBM/lib_lightgbm.so /usr/lib/ \
 && git clone https://github.com/remnrem/luna-base.git \
 && cd luna-base \ 
 && make -j 2 LGBM=1 LGBM_PATH=/Programme/LightGBM \
 && ln -s /Programme/luna-base/luna /usr/local/bin/luna \
 && ln -s /Programme/luna-base/destrat /usr/local/bin/destrat \
 && ln -s /Programme/luna-base/behead /usr/local/bin/behead \
 && ln -s /Programme/luna-base/fixrows /usr/local/bin/fixrows

## LunaR
RUN git clone https://github.com/remnrem/luna.git \
 && echo 'PKG_LIBS=include/libluna.a -L${LGBM_PATH} -lfftw3 -l_lightgbm' >> luna/src/Makevars \
 && LGBM=1 LGBM_PATH=/Programme/LightGBM/ EXTRA_PKG_LIBS="-L/Programme/LightGBM/lib -l_lightgbm" R CMD INSTALL luna

#------------------------------- Multi-stage build (keeps the image size down)-------------------------------------------
FROM rocker/r-ver:4.2.1
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    g++ \
    git \
    nano \
    wget \
    zlib1g-dev \
    fftw3-dev \
    libgit2-dev \
    libssl-dev \
    libssh2-1-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libomp-dev \
    cmake \
    libxt6 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /data
RUN mkdir /root/hypnoscope
COPY ui.R server.R helpers.R /root/hypnoscope/
COPY data /root/hypnoscope/data

ENV _R_SHLIB_STRIP_=true
COPY --from=builder /Programme/luna-base/luna /usr/local/bin/luna
COPY --from=builder /Programme/luna-base/destrat /usr/local/bin/destrat
COPY --from=builder /Programme/luna-base/behead /usr/local/bin/behead
COPY --from=builder /Programme/luna-base/fixrows /usr/local/bin/fixrows
COPY --from=builder /Programme/LightGBM/lib_lightgbm.so /usr/local/lib
COPY --from=builder /Programme/LightGBM/lib_lightgbm.so /usr/lib
COPY --from=builder /usr/local/lib/R /usr/local/lib/R

EXPOSE 3838
CMD ["R", "-q", "-e", "shiny::runApp('/root/hypnoscope')"]
