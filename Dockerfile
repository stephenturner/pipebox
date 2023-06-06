FROM condaforge/mambaforge

COPY ./environment.yml /environment.yml

RUN mamba env update --file /environment.yml

RUN apt update && apt install -y vim gcc make zlib1g-dev

ARG VERSION_SEQTK="1.4"
RUN wget -q https://github.com/lh3/seqtk/archive/refs/tags/v${VERSION_SEQTK}.tar.gz && \
    tar xzf v${VERSION_SEQTK}.tar.gz && \
    cd seqtk-${VERSION_SEQTK} && \
    make && \
    make install

COPY ./src /src

ENTRYPOINT ["/bin/bash", "/src/pipebox.sh"]