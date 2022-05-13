FROM ubuntu:jammy

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Get dependencies
RUN RUN set -ex; \
    apt-get -y update; \
    apt-get -y install \
        build-essential \
        cmake \
        libhts-dev; \
    rm -rf /var/lib/apt/lists/*

# Install samtag
ARG architecture=broadwell
COPY . /opt/samtag
WORKDIR /opt/samtag/build
RUN set -ex; \
    cmake \
        -DCOMPILER_ARCHITECTURE=${architecture} \
        ..; \
    make install

ENV PATH="/opt/samtag:${PATH}"
ENTRYPOINT ["samtag"]
