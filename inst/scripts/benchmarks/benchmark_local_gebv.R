pkg_name = "HapSelect"
has_installed_pkg = requireNamespace(pkg_name, quietly = TRUE)

if (has_installed_pkg) {
  suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
} else {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload to run this benchmark from the HapSelect source tree.")
  }
  pkgload::load_all(".", quiet = TRUE)
}
source(file.path("inst", "scripts", "benchmarks", "benchmark_support.R"))

# Settings ---------------------------------------------------------------------
defaults = list(
  n_markers = 500L,
  n_individuals = 200L,
  n_chr = 5L,
  n_blocks = 50L,
  missing_rate = 0.02,
  seed = 1L,
  n_reps = 3L,
  chunk_size = 100L
)

parse_args = function(defaults) {
  args = commandArgs(trailingOnly = TRUE)
  values = defaults

  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Arguments must be provided as --name=value")
    }

    key = sub("^--([^=]+)=.*$", "\\1", arg)
    value = sub("^--[^=]+=(.*)$", "\\1", arg)

    if (!key %in% names(values)) {
      stop("Unknown argument: ", key)
    }

    if (is.integer(values[[key]])) {
      values[[key]] = as.integer(value)
    } else if (is.numeric(values[[key]])) {
      values[[key]] = as.numeric(value)
    } else {
      values[[key]] = value
    }
  }

  values
}

cfg = parse_args(defaults)

n_markers     = cfg$n_markers
n_individuals = cfg$n_individuals
n_chr         = cfg$n_chr
n_blocks      = cfg$n_blocks
missing_rate  = cfg$missing_rate
seed          = cfg$seed
n_reps        = cfg$n_reps
chunk_size    = cfg$chunk_size

# Helpers ----------------------------------------------------------------------
simulate_marker_effects = function(geno, seed) {
  set.seed(seed)
  data.frame(
    SNP    = geno[, 1],
    Effect = rnorm(nrow(geno), mean = 0, sd = 0.1),
    stringsAsFactors = FALSE
  )
}

simulate_haploblocks = function(geno, n_blocks, seed) {
  set.seed(seed)

  chrom_markers = split(geno[, 1], geno[, 2])
  chrom_names = names(chrom_markers)
  chrom_sizes = lengths(chrom_markers)

  base_blocks = rep(floor(n_blocks / length(chrom_markers)), length(chrom_markers))
  if (sum(base_blocks) == 0L) {
    base_blocks[] = 1L
  }

  remainder = n_blocks - sum(base_blocks)
  if (remainder > 0L) {
    order_idx = order(chrom_sizes, decreasing = TRUE)
    base_blocks[order_idx[seq_len(remainder)]] = base_blocks[order_idx[seq_len(remainder)]] + 1L
  }
  base_blocks = pmax(base_blocks, 1L)
  names(base_blocks) = chrom_names

  block_rows = purrr::imap_dfr(chrom_markers, function(markers, chrom_name) {
    n_chr_blocks = min(length(markers), base_blocks[[chrom_name]])
    block_assign = sample(rep(seq_len(n_chr_blocks), length.out = length(markers)))
    block_strings = tapply(markers, block_assign, function(m) paste(m, collapse = ";"))

    data.frame(
      Chrom = geno[match(markers[1], geno[, 1]), 2],
      Block_ID = paste0("Chr", chrom_name, ":B", seq_along(block_strings)),
      Block = as.character(block_strings),
      stringsAsFactors = FALSE
    )
  })

  block_rows
}

make_row = function(method, times, result) {
  data.frame(
    Method         = method,
    Mean_s         = round(mean(times), 3),
    Min_s          = round(min(times),  3),
    Max_s          = round(max(times),  3),
    Blocks         = nrow(result$Haploblocks),
    Blocks_per_sec = round(nrow(result$Haploblocks) / mean(times)),
    stringsAsFactors = FALSE
  )
}

# Data -------------------------------------------------------------------------
geno           = simulate_genotypes(n_markers, n_individuals, n_chr, missing_rate, seed)
haploblocks_df = simulate_haploblocks(geno, n_blocks, seed)
marker_effects = simulate_marker_effects(geno, seed)

cat(
  "Markers:      ", n_markers,     "\n",
  "Individuals:  ", n_individuals, "\n",
  "Chromosomes:  ", n_chr,         "\n",
  "Blocks:       ", nrow(haploblocks_df), "\n",
  "Missing rate: ", missing_rate,  "\n",
  "Chunk size:   ", chunk_size,    "\n",
  "Reps:         ", n_reps,        "\n\n",
  sep = ""
)

# Run --------------------------------------------------------------------------
progressr::handlers("void")

cat("Benchmarking compute_local_GEBV (parallel = FALSE)\n")
benchmark_serial = time_reps(n_reps, function() {
  compute_local_GEBV(
    geno           = geno,
    marker_effects = marker_effects,
    haploblocks_df = haploblocks_df,
    set_missing_NA = TRUE,
    mean_adjust    = TRUE,
    parallel       = FALSE
  )
})
cat("  Elapsed (s): ", paste(round(benchmark_serial$times, 3), collapse = ", "),
    "  |  mean: ", round(mean(benchmark_serial$times), 3), "s\n", sep = "")

cat("Benchmarking compute_local_GEBV (parallel = TRUE)\n")
if (.Platform$OS.type == "windows" && !has_installed_pkg) {
  cat(
    "  Skipped: Windows parallel workers need an installed HapSelect package.\n",
    "  Run `Rscript -e \"install.packages('.', repos = NULL, type = 'source')\"` from the repo root,\n",
    "  then rerun this benchmark.\n",
    sep = ""
  )
  benchmark_parallel = NULL
} else {
  benchmark_parallel = tryCatch(
    time_reps(n_reps, function() {
      compute_local_GEBV(
        geno           = geno,
        marker_effects = marker_effects,
        haploblocks_df = haploblocks_df,
        set_missing_NA = TRUE,
        mean_adjust    = TRUE,
        parallel       = TRUE,
        chunk_size     = chunk_size
      )
    }),
    error = function(e) {
      cat("  Parallel benchmark failed: ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
}

if (!is.null(benchmark_parallel)) {
  cat("  Elapsed (s): ", paste(round(benchmark_parallel$times, 3), collapse = ", "),
      "  |  mean: ", round(mean(benchmark_parallel$times), 3), "s\n", sep = "")
}

# Summary ----------------------------------------------------------------------
rows = list(make_row("LocalGEBV (serial)", benchmark_serial$times, benchmark_serial$result))
if (!is.null(benchmark_parallel)) {
  rows[[length(rows) + 1L]] = make_row("LocalGEBV (parallel)", benchmark_parallel$times, benchmark_parallel$result)
}

summary_df = do.call(rbind, rows)
row.names(summary_df) = NULL

cat("\nBenchmark summary (", n_reps, " reps each)\n", sep = "")
print(summary_df, row.names = FALSE)

if (!is.null(benchmark_parallel)) {
  speedup = mean(benchmark_serial$times) / mean(benchmark_parallel$times)
  cat("\nParallel speedup: ", round(speedup, 2), "x\n", sep = "")
}
