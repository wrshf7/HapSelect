# Development
Details for contributing to this package

## System Dependencies
Installs PLINK, GenomicSelection and all system dependencies

### Linux (Ubuntu/Debian)
`./inst/scripts/install/install_linux.sh`
### Windows
`Powershell -ExecutionPolicy bypass -File inst/scripts/install/install_windows.ps1`

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
- `HapSelect::pairwise_ld()` runtime
- `HapSelect::plink_pairwise_ld()` runtime
