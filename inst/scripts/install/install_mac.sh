#!/usr/bin/env bash

set -euo pipefail

# Resolve the package root from this script location.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

PLINK_URL="https://s3.amazonaws.com/plink1-assets/plink_mac_20250819.zip"
PLINK_ARCHIVE="plink_mac_20250819.zip"

# Check platform support.
if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This script only supports macOS." >&2
  exit 1
fi

# Check that Rscript is available.
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is required but was not found. Please ensure R is installed before running this script." >&2
  exit 1
fi

# Check that the Xcode Command Line Tools are available for package compilation.
if ! xcode-select -p >/dev/null 2>&1; then
  echo "Xcode Command Line Tools are required but were not found." >&2
  echo "Run 'xcode-select --install', finish the installation, and then rerun this script." >&2
  exit 1
fi

# Check that Homebrew is available for dependency installs.
if ! command -v brew >/dev/null 2>&1; then
  echo "Homebrew is required for this script. Install it from https://brew.sh and rerun the script." >&2
  exit 1
fi

# Install system dependencies used by development tools and package builds.
brew update
brew install \
  cmake \
  openssl@3 \
  curl \
  fontconfig \
  harfbuzz \
  fribidi \
  freetype \
  libpng \
  libtiff \
  jpeg-turbo \
  libxml2 \
  git \
  cairo \
  libuv \
  unixodbc \
  gdal \
  geos \
  proj \
  sqlite \
  wget \
  unzip \
  gcc

# Download and install PLINK into the user's bin directory.
mkdir -p "$HOME/bin"
if [[ -f "$HOME/.zprofile" ]]; then
  # shellcheck disable=SC1090
  source "$HOME/.zprofile"
elif [[ -f "$HOME/.bash_profile" ]]; then
  # shellcheck disable=SC1090
  source "$HOME/.bash_profile"
fi
curl -L "$PLINK_URL" -o "$HOME/bin/$PLINK_ARCHIVE"
unzip -o "$HOME/bin/$PLINK_ARCHIVE" -d "$HOME/bin/"
rm -f "$HOME/bin/$PLINK_ARCHIVE"
chmod +x "$HOME/bin/plink" 2>/dev/null || true
export PATH="$HOME/bin:$PATH"

if ! plink --version; then
  if [[ "$(uname -m)" == "arm64" ]]; then
    echo "PLINK could not be executed. On Apple Silicon, install Rosetta 2 and rerun this script:" >&2
    echo "  softwareupdate --install-rosetta --agree-to-license" >&2
  fi
  exit 1
fi

# Install the R dependencies required by HapSelect.
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = 'https://cloud.r-project.org'); remotes::install_deps('$PACKAGE_DIR', dependencies = TRUE)"

# Install the genomicSimulation R dependency from GitHub.
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = 'https://cloud.r-project.org'); remotes::install_github('vllrs/genomicSimulation')"
