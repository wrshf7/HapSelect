# Only load the package and support helpers when run directly (Rscript benchmark_local_gebv.R).
# When source()'d by benchmark_batched.R, sys.nframe() > 0 and these are already set up.
if (sys.nframe() == 0L) {
  pkg_name = "HapSelect"
  if (requireNamespace(pkg_name, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
  } else {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Install pkgload to run this benchmark from the HapSelect source tree.")
    }
    pkgload::load_all(".", quiet = TRUE)
  }
  source(file.path("inst", "scripts", "benchmarks", "benchmark_support.R"))
}

run_benchmark_local_gebv = function(params = list()) {
  params = coerce_params(params, list(
    n_markers     = 500L,
    n_individuals = 200L,
    n_chr         = 5L,
    n_blocks      = 50L,
    missing_rate  = 0.02,
    seed          = 1L,
    n_reps        = 3L,
    chunk_size    = 100L
  ))

  has_installed_pkg = requireNamespace("HapSelect", quietly = TRUE)

  n_markers     = params$n_markers
  n_individuals = params$n_individuals
  n_chr         = params$n_chr
  n_blocks      = params$n_blocks
  missing_rate  = params$missing_rate
  seed          = params$seed
  n_reps        = params$n_reps
  chunk_size    = params$chunk_size

  geno           = simulate_genotypes(n_markers, n_individuals, n_chr, missing_rate, seed)
  haploblocks_df = simulate_haploblocks(geno, n_blocks, seed)
  marker_effects = simulate_marker_effects(geno, seed)

  cat(
    "Markers:      ", n_markers,          "\n",
    "Individuals:  ", n_individuals,      "\n",
    "Chromosomes:  ", n_chr,              "\n",
    "Blocks:       ", nrow(haploblocks_df), "\n",
    "Missing rate: ", missing_rate,       "\n",
    "Chunk size:   ", chunk_size,         "\n",
    "Reps:         ", n_reps,             "\n\n",
    sep = ""
  )

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
  if (!has_installed_pkg) {
    cat(
      "  Skipped: parallel workers need an installed HapSelect package (multisession spawns fresh R processes).\n",
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

  make_row = function(method, times, result) {
    data.frame(
      Method         = method,
      Mean_s         = round(mean(times), 3),
      Min_s          = round(min(times),  3),
      Max_s          = round(max(times),  3),
      Blocks         = nrow(result$Haploblocks),
      Blocks_per_sec = round(nrow(result$Haploblocks) / mean(times), 3),
      stringsAsFactors = FALSE
    )
  }

  rows = list(make_row("LocalGEBV (serial)", benchmark_serial$times, benchmark_serial$result))
  if (!is.null(benchmark_parallel)) {
    rows[[length(rows) + 1L]] = make_row("LocalGEBV (parallel)", benchmark_parallel$times, benchmark_parallel$result)
  }

  summary_df = do.call(rbind, rows)
  row.names(summary_df) = NULL

  cat("\nBenchmark summary (", n_reps, " reps each)\n", sep = "")
  print(summary_df, row.names = FALSE)

  parallel_speedup = NULL
  if (!is.null(benchmark_parallel)) {
    parallel_speedup = round(mean(benchmark_serial$times) / mean(benchmark_parallel$times), 2)
    cat("\nParallel speedup: ", parallel_speedup, "x\n", sep = "")
  }

  list(
    benchmark        = "local_gebv",
    timestamp        = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    params           = params,
    results          = lapply(seq_len(nrow(summary_df)), function(i) as.list(summary_df[i, ])),
    parallel_speedup = parallel_speedup
  )
}

# Only execute when run directly, not when source()'d to load run_benchmark_local_gebv().
if (sys.nframe() == 0L) {
  run_benchmark_local_gebv(parse_args(list(
    n_markers     = 500L,
    n_individuals = 200L,
    n_chr         = 5L,
    n_blocks      = 50L,
    missing_rate  = 0.02,
    seed          = 1L,
    n_reps        = 3L,
    chunk_size    = 100L
  )))
}
