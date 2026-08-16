ARG JULIA_VERSION=1.12.6
FROM julia:${JULIA_VERSION}-bookworm

LABEL org.opencontainers.image.title="MOOSE"
LABEL org.opencontainers.image.description="Julia container for Mock Observation Of Synchrotron Emission"

ENV JULIA_DEPOT_PATH=/opt/julia-depot: \
    JULIA_PKG_PRECOMPILE_AUTO=0 \
    JULIA_NUM_THREADS=auto \
    PYTHONUNBUFFERED=1 \
    MOOSE_PROJECT=/app

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        fontconfig \
        git \
        libgl1 \
        libx11-6 \
        libxext6 \
        libxrender1 \
        libxt6 \
        python3 \
        tini \
    && rm -rf /var/lib/apt/lists/*

# Only Project.toml is copied: the Manifest is resolved at build time so the
# image tracks the compat bounds rather than a pinned resolution. Keeping this
# COPY before the sources preserves Docker layer caching for the dependencies.
COPY Project.toml ./

RUN for attempt in 1 2 3; do \
        julia --startup-file=no --project=/app -e 'using Pkg; Pkg.instantiate()' && break; \
        if [ "$attempt" = "3" ]; then exit 1; fi; \
        sleep 10; \
    done

COPY README.md setup.jl ./
COPY config ./config
COPY python ./python
COPY src ./src
COPY test ./test

RUN julia --startup-file=no --project=/app -e 'using Pkg; Pkg.precompile(); using Moose'
RUN python3 -m unittest discover -s test -p 'test_*.py'

RUN mkdir -p /data

VOLUME ["/data"]

ENTRYPOINT ["tini", "--"]

CMD ["julia", "--startup-file=no", "--project=/app", "-e", "using Moose; run_moose(help=true)"]
