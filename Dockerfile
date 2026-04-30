############################################################
# Dockerfile to build CRISPResso2 (multi-stage)
############################################################

# --- Build stage: resolve deps and install package ---
FROM ghcr.io/prefix-dev/pixi:0.43.0 AS build

USER root

RUN apt-get update \
  && apt-get install -y --no-install-recommends gcc g++ \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /CRISPResso2
COPY pixi.toml pyproject.toml ./
RUN pixi install

COPY . .
RUN pixi run -e default install

# --- Runtime stage: slim image with runtime environment only ---
FROM ubuntu:24.04

LABEL org.opencontainers.image.authors="support@edilytics.com"

# Install MS core fonts for plots
RUN echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
  && apt-get update \
  && apt-get install -y --no-install-recommends ttf-mscorefonts-installer \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && rm -rf /usr/share/zoneinfo

# Keep env path identical to build stage to avoid shebang/entrypoint path issues.
COPY --from=build /CRISPResso2/.pixi/envs/default /CRISPResso2/.pixi/envs/default

# Copy only minimal app files required at runtime.
COPY --from=build /CRISPResso2/CRISPResso2 /CRISPResso2/CRISPResso2
COPY --from=build /CRISPResso2/CRISPResso2.egg-info /CRISPResso2/CRISPResso2.egg-info
COPY --from=build /CRISPResso2/CRISPResso2_router.py /CRISPResso2/CRISPResso2_router.py
COPY --from=build /CRISPResso2/LICENSE.txt /CRISPResso2/LICENSE.txt

WORKDIR /CRISPResso2

ENV PATH="/CRISPResso2/.pixi/envs/default/bin:${PATH}"

# Verify
RUN CRISPResso -h \
  && CRISPRessoBatch -h \
  && CRISPRessoPooled -h \
  && CRISPRessoWGS -h \
  && CRISPRessoCompare -h

ENTRYPOINT ["python", "/CRISPResso2/CRISPResso2_router.py"]
