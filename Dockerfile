FROM condaforge/mambaforge

# Copy conda environment yaml from the directory on the host to the root inside the container
COPY ./environment.yml /environment.yml

# Rather than creating a new environment, just update the base environment. 
# The base environment is automatically activated when the container is run.
RUN mamba env update --file /environment.yml

# Get some basic utilities 
RUN apt update && apt install -y vim gcc make zlib1g-dev

# Seqtk is available via conda, but let's demonstrate how to build something from source that isn't.
ARG VERSION_SEQTK="1.4"
RUN wget -q https://github.com/lh3/seqtk/archive/refs/tags/v${VERSION_SEQTK}.tar.gz && \
    tar xzf v${VERSION_SEQTK}.tar.gz && \
    cd seqtk-${VERSION_SEQTK} && \
    make && \
    make install

# Copy assets from the directory on the host to the root inside the container
COPY ./src /src

# The entrypoint into the container is the pipebox script
ENTRYPOINT ["/bin/bash", "/src/pipebox.sh"]