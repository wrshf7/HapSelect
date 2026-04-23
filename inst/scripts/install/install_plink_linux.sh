#!/usr/bin/env bash

set -euo pipefail

mkdir -p "$HOME/bin"
source "$HOME/.profile" 2>/dev/null || true

echo "Downloading PLINK"
wget -P "$HOME/bin/" https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip

echo "Extracting PLINK into $HOME/bin"
cd "$HOME/bin"
unzip -o plink_linux_x86_64_20250819.zip -d "$HOME/bin/"
rm -f plink_linux_x86_64_20250819.zip

export PATH="$PATH:$HOME/bin"

echo "Installed PLINK to $HOME/bin"
plink --version
