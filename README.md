# HapSelect
Minimal R package for haplotype block workflows, localGEBV analysis, and parent selection for genomic selection.

## Contributing
For details on contributing to this package, see [the development guide](DEVELOPMENT.md).

## Install

From the package root in R:

```r
install.packages("devtools")
devtools::install(".")
```
Alternatively, download a .zip archive, decompress the .zip file, and run in R:

```r
install.packages("/path/to/extracted/folder/HapSelect/HapSelect", repos = NULL) 
```

## Run Example Workflow

```r
library(HapSelect)
source(system.file("examples", "example_workflow_minimal_comments.R", package = "HapSelect"))
```
