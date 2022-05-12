FROM alpine:3.15

RUN set -ex; \
    apk add --no-cache \
        autoconf \
        make \
        automake \
        zlib-dev bzip2-dev xz-dev libcurl openssl-dev \
        libtool \
        gcc g++ \
        cmake \
        wget; \
    rm -rf /var/lib/apt/lists/*

# Install htslib
ARG HTSLIB_VERSION='1.15.1'
ARG HTSCODECS_VERSION='1.2.2'
WORKDIR /opt
# hadolint ignore=DL3003,DL4006
RUN set -ex; \
    wget -qO- https://github.com/samtools/htslib/archive/refs/tags/${HTSLIB_VERSION}.tar.gz | tar xzf - -C $(pwd); \
    cd htslib-${HTSLIB_VERSION}; \
    wget -qO- https://github.com/samtools/htscodecs/archive/refs/tags/v${HTSCODECS_VERSION}.tar.gz | tar xzf - -C $(pwd); \
    rm -r htscodecs && mv htscodecs-${HTSCODECS_VERSION} htscodecs; \
    cd htscodecs; \
    autoreconf -i; \
    ./configure; \
    make; \
    cd ..; \
    autoreconf -i; \
    ./configure; \
    make install

# Install samtag
COPY . /opt/samtag
WORKDIR /opt/samtag/build
RUN set -ex; \
    cmake ..; \
    make install

ENV PATH="/opt/samtag:${PATH}"
ENTRYPOINT ["samtag"]
