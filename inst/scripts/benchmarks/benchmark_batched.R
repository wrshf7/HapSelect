# Run a sequence of benchmarks defined in a JSON config file and save results.
#
# Usage (from repo root):
#   Rscript inst/scripts/benchmarks/benchmark_batched.R [--config=FILE] [--output=FILE]
#
# Defaults:
#   --config  inst/scripts/benchmarks/benchmark_batched_example.json
#   --output  benchmark_results.json
#
# Config format: JSON array of job objects, each with:
#   "benchmark" : "ld" | "haploblocks" | "local_gebv"
#   "params"    : object of parameter overrides (all optional; missing keys use defaults)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop('Install jsonlite to run batch benchmarks: install.packages("jsonlite")')
}

if (!file.exists("DESCRIPTION")) {
  stop("Run benchmark_batched.R from the HapSelect repository root.")
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload to run this benchmark from the HapSelect source tree.")
}
pkgload::load_all(".", quiet = TRUE)

benchmarks_dir = file.path("inst", "scripts", "benchmarks")
source(file.path(benchmarks_dir, "benchmark_support.R"))

source(file.path(benchmarks_dir, "benchmark_ld.R"))
source(file.path(benchmarks_dir, "benchmark_haploblocks.R"))
source(file.path(benchmarks_dir, "benchmark_local_gebv.R"))

# Parse --config / --output args -----------------------------------------------
args        = commandArgs(trailingOnly = TRUE)
config_path = file.path(benchmarks_dir, "benchmark_batched_example.json")
output_path = "benchmark_results.json"

for (arg in args) {
  if (startsWith(arg, "--config=")) {
    config_path = sub("^--config=", "", arg)
  } else if (startsWith(arg, "--output=")) {
    output_path = sub("^--output=", "", arg)
  } else {
    stop("Unknown argument: ", arg, ". Supported: --config=FILE  --output=FILE")
  }
}

if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}

# Load and validate config -----------------------------------------------------
jobs = jsonlite::fromJSON(config_path, simplifyVector = FALSE)

if (!is.list(jobs) || length(jobs) == 0L) {
  stop("Config must be a non-empty JSON array of job objects.")
}

valid_benchmarks = c("ld", "haploblocks", "local_gebv")

for (i in seq_along(jobs)) {
  job = jobs[[i]]
  if (is.null(job$benchmark)) {
    stop("Job ", i, " is missing required field 'benchmark'.")
  }
  if (!job$benchmark %in% valid_benchmarks) {
    stop("Job ", i, ": unknown benchmark '", job$benchmark,
         "'. Valid options: ", paste(valid_benchmarks, collapse = ", "))
  }
}

# Run --------------------------------------------------------------------------
cat("Config:  ", config_path,  "\n")
cat("Output:  ", output_path,  "\n")
cat("Jobs:    ", length(jobs), "\n\n")

all_results = vector("list", length(jobs))

for (i in seq_along(jobs)) {
  job    = jobs[[i]]
  bench  = job$benchmark
  params = if (is.null(job$params)) list() else job$params

  cat(sprintf("[%d/%d] %s\n", i, length(jobs), bench))
  cat(strrep("-", 60), "\n", sep = "")

  all_results[[i]] = tryCatch(
    switch(
      bench,
      ld          = run_benchmark_ld(params),
      haploblocks = run_benchmark_haploblocks(params),
      local_gebv  = run_benchmark_local_gebv(params)
    ),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      list(
        benchmark = bench,
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        params    = params,
        error     = conditionMessage(e)
      )
    }
  )

  cat("\n")
}

# Save results -----------------------------------------------------------------
jsonlite::write_json(all_results, output_path, pretty = TRUE, auto_unbox = TRUE)
cat("Results written to: ", output_path, "\n", sep = "")
