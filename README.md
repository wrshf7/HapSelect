# FastStack
Minimal R package for haplotype block workflows and localGEBV analysis.

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
install.packages("/path/to/extracted/folder/FastStack/", repos = NULL) 
```

## Run Example Workflow

```r
library(FastStack)
source(system.file("examples", "example_workflow_minimal_comments.R", package = "FastStack"))
```
