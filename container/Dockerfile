FROM r-base:4.4.0

RUN apt-get update && apt-get install -yq --no-install-recommends \
    build-essential \
    lsb-release \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    wget \
    libcairo2-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libxml2-dev \
    libxml2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('BiocManager')"

RUN R -e "BiocManager::install(c('optparse', 'scDblFinder', 'Seurat', 'data.table', 'AUCell', 'NMF'))"

COPY ./scripts /scripts