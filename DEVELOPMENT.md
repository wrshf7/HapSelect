# Development
Details for contributing to this package

## System Dependencies
Development in this package uses `devtools` which itself has several system level dependencies

### Linux (Ubuntu/Debian)
**The following dependencies are required to use `devtools`:**

```
# Core build tools
sudo apt install build-essential gfortran cmake

# SSL / curl (very common)
sudo apt install libssl-dev libcurl4-openssl-dev

# Font / graphics stack (devtools, ragg, systemfonts, etc.)
sudo apt install libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

# XML (xml2, rvest, many tidyverse deps)
sudo apt install libxml2-dev

# Git (devtools needs it)
sudo apt install git

# Cairo graphics (common for plotting packages)
sudo apt install libcairo2-dev

# fs package
sudo apt install libuv1-dev

# Database drivers (if you ever use DBI, odbc, etc.)
sudo apt install unixodbc-dev

# Spatial packages (sf, terra, etc.) - optional but common
sudo apt install libgdal-dev libgeos-dev libproj-dev libsqlite3-dev
```

### Windows
```TODO - Probably less painful than linux```


## PLINK
The PLINK package is an open source genome analysis tool set which can be used as an alternative to the FastStack implementation for operation like calculating LD, but requires additional installation steps.

### Linux (Ubuntu/Debian)
```
# Setup folders and download archive
mkdir -p ~/bin && source ~/.profile
wget -P ~/bin/ https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip

# Unzip the archive
cd ~/bin/
unzip plink_linux_x86_64_20250819.zip -d ~/bin/
rm plink_linux_x86_64_20250819.zip

# Path update
export PATH=$PATH:~/bin

# Test installation
plink --version
```

### Windows
```TODO```

## Package Dependencies
### Devtools
Devtools is used to install dependencies and run tests
Install devtools using `install.packages("devtools")`

### FastStack Dependencies
Next install the FastStack package level dependencies using devtools
```
devtools::install_deps(dependencies = TRUE)
```

## Running Tests
Now that dependencies are installed, run the tests:
```
devtools::test()
```

## Benchmarking LD Implementations
To compare the internal LD calculation against the PLINK-backed implementation on a synthetic dataset:
```
install.packages("pkgload")
Rscript inst/scripts/benchmark_ld.R --n_markers=500 --n_individuals=200 --n_chr=5 --missing_rate=0.02 --seed=1
```

This benchmark reports:
- `FastStack::pairwise_ld()` runtime
- `FastStack::plink_pairwise_ld()` runtime
