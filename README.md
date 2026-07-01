# HapSelect

HapSelect is an R package for haplotype-based genomic selection. It partitions the genome into linkage disequilibrium blocks (haploblocks), estimates per-block breeding value contributions (localGEBV), uses a genetic algorithm to select parents that maximise coverage of high-value haplotype alleles, and runs a basic simulation comparing genetic algorithm parents to truncation selection over time.

--- 

![HapSelect process overview](man/figures/overview_diagram.png)


---


## Documentation
The latest documentation for HapSelect is available on the **[HapSelect documentation site](https://wrshf7.github.io/HapSelect-Docs)**.

## Important Papers
<strong>localGEBV Method and Haploblock Formation:</strong><br>
[Shaffer et al. 2025. Local genomic estimates provide a powerful framework for haplotype discovery. bioRxiv (under review).](https://doi.org/10.1101/2025.08.28.672830)

<strong>Origin of the localGEBV Method and Parent Optimization with a Genetic Algorithm:</strong><br>
[Kemper et al. 2012. Long-term selection strategies for complex traits using high-density genetic markers. J Dairy Sci.](https://doi.org/10.3168/jds.2011-5289)

<strong>The First Implementation of Haploblocking with localGEBV:</strong><br>
[Voss-Fels et al. 2019. Breeding improves wheat productivity under contrasting agrochemical input levels. Nat Plants.](https://doi.org/10.1038/s41477-019-0445-5)

<strong>The Concept of the Ultimate Genotype:</strong><br>
[Hays et al. 2024. Potential approaches to create ultimate genotypes in crops and livestock. Nat Genet.](https://doi.org/10.1038/s41588-024-01942-0)

<strong>genomicSimulation Package Description and Simulation Comparing the Genetic Algorithm to Truncation Selection:</strong><br>
[Villiers et al. 2024. Evolutionary computing to assemble standing genetic diversity and achieve long-term genetic gain. Plant Genome.](https://doi.org/10.1002/tpg2.20467)


## Workflow at a glance

| Stage | What it does | Key function |
|-------|-------------|--------------|
| Pairwise LD | Compute r² between all marker pairs | `pairwise_ld()` / `plink_pairwise_ld()` |
| Haploblocking | Partition genome into LD-based haploblocks | `def_blocks()` |
| Genomic Prediction | Compute marker effects and prediction accuracy | `create_marker_effects_file()`, `cross_validation()` |
| LocalGEBV | Estimate per-block breeding value per individual | `compute_local_GEBV()` |
| Visualisation | Explore haploblock structure, localGEBV patterns | `plot_haploblocks()`, `plot_ld_decay()`, … |
| Parent selection | Optimise founder set using a genetic algorithm | `local_gebv_parent_selection()`, `haplotype_parent_selection()` |
| Basic Simulation | Compare GA and TS parent performance over time | `localGEBV_vs_TS_simulation()`, `Haplotype_vs_TS_simulation()` |

For full documentation, workflow guides, an in-depth installation guide, and parameter details, see the **[HapSelect documentation site](https://wrshf7.github.io/HapSelect-Docs)**.

## Install

```r
install.packages("devtools")
devtools::install("path/to/unzipped/HapSelect")
```

> **Note:** [genomicSimulation](https://github.com/vllrs/genomicSimulation) and [RTools 4.5](https://cran.r-project.org/) must be installed separately. PLINK 1.9 is optional but recommended for LD computation. See the [installation guide](https://wrshf7.github.io/HapSelect-Docs/installation) for details.

## Run Example Workflow

```r
library(HapSelect)
source(system.file("examples", "example_workflow_minimal_comments.R", package = "HapSelect"))
```

## Contributing

For details on contributing to this package, see [the development guide](DEVELOPMENT.md).

## Authors

Will Shaffer, Zane Carter, Victor Papin

The University of Queensland
