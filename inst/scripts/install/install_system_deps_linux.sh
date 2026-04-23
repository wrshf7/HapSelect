#!/usr/bin/env bash

set -euo pipefail

if ! command -v apt-get >/dev/null 2>&1; then
  echo "This script currently supports Ubuntu/Debian systems with apt-get." >&2
  echo "Beginning install..." >&2
  exit 1
fi

sudo apt-get update
sudo apt-get install -y \
  build-essential \
  gfortran \
  cmake \
  libssl-dev \
  libcurl4-openssl-dev \
  libfontconfig1-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libxml2-dev \
  git \
  libcairo2-dev \
  libuv1-dev \
  unixodbc-dev \
  libgdal-dev \
  libgeos-dev \
  libproj-dev \
  libsqlite3-dev \
  wget \
  unzip

echo "System dependencies installed."
