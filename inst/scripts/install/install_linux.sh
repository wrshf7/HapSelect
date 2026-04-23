#!/usr/bin/env bash

set -euo pipefail

# Resolve the package root from this script location.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Check that Rscript is available.
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is required but was not found. Please ensure R is installed before running this script." >&2
  exit 1
fi

# Check platform support.
if ! command -v apt-get >/dev/null 2>&1; then
  echo "This script currently supports Ubuntu/Debian systems with apt-get." >&2
  exit 1
fi

# Install system dependencies used by development tools and package builds.
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

# Download and install PLINK into the user's bin directory.
mkdir -p "$HOME/bin"
source "$HOME/.profile" 2>/dev/null || true
wget -P "$HOME/bin/" https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip
cd "$HOME/bin"
unzip -o plink_linux_x86_64_20250819.zip -d "$HOME/bin/"
rm -f plink_linux_x86_64_20250819.zip
export PATH="$PATH:$HOME/bin"
plink --version

# Install the R dependencies required by HapSelect.
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = 'https://cloud.r-project.org'); remotes::install_deps('$PACKAGE_DIR', dependencies = TRUE)"

# Install the genomicSimulation R dependency from GitHub.
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = 'https://cloud.r-project.org'); remotes::install_github('vllrs/genomicSimulation')"
