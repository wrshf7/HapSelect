# Only load the package and support helpers when run directly (Rscript benchmark_ld.R).
# When source()'d by benchmark_batched.R, sys.nframe() > 0 and these are already set up.
if (sys.nframe() == 0L) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload to run this benchmark from the HapSelect source tree.")
  }
  pkgload::load_all(".", quiet = TRUE)
  source(file.path("inst", "scripts", "benchmarks", "benchmark_support.R"))
}

run_benchmark_ld = function(params = list()) {
  params = coerce_params(params, list(
    n_markers     = 500L,
    n_individuals = 200L,
    n_chr         = 5L,
    missing_rate  = 0.02,
    seed          = 1L,
    n_reps        = 3L
  ))

  n_markers     = params$n_markers
  n_individuals = params$n_individuals
  n_chr         = params$n_chr
  missing_rate  = params$missing_rate
  seed          = params$seed
  n_reps        = params$n_reps

  geno = simulate_genotypes(n_markers, n_individuals, n_chr, missing_rate, seed)
  expected_pairs = sum(vapply(
    split(geno[[1]], geno[[2]]),
    function(x) choose(length(x), 2),
    numeric(1)
  ))

  work_dir = tempfile("hapselect_ld_benchmark_")
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE))
  dir.create(work_dir)

  text_prefix = file.path(work_dir, "synthetic_ld")
  bed_prefix  = file.path(work_dir, "synthetic_ld_bin")

  cat(
    "Markers:      ", n_markers,     "\n",
    "Individuals:  ", n_individuals, "\n",
    "Chromosomes:  ", n_chr,         "\n",
    "Missing rate: ", missing_rate,  "\n",
    "Reps:         ", n_reps,        "\n\n",
    sep = ""
  )

  cat("Writing PLINK PED/MAP files\n")
  write_plink_text_files(geno, text_prefix)

  cat("Converting to PLINK binary fileset\n")
  run_plink(c("--file", text_prefix, "--make-bed", "--out", bed_prefix))

  cat("\nBenchmarking HapSelect pairwise_ld (parallelize = FALSE)\n")
  benchmark_serial = time_reps(n_reps, function() pairwise_ld(geno, parallelize = FALSE))
  cat("  Elapsed (s): ", paste(round(benchmark_serial$times, 3), collapse = ", "),
      "  |  mean: ", round(mean(benchmark_serial$times), 3), "s\n", sep = "")

  cat("Benchmarking HapSelect pairwise_ld (parallelize = TRUE)\n")
  benchmark_parallel = tryCatch(
    time_reps(n_reps, function() pairwise_ld(geno, parallelize = TRUE)),
    error = function(e) {
      cat("  Skipped: parallel workers cannot load an uninstalled package.\n",
          "  Install HapSelect with install.packages('.', repos = NULL) and re-run.\n", sep = "")
      NULL
    }
  )
  if (!is.null(benchmark_parallel)) {
    cat("  Elapsed (s): ", paste(round(benchmark_parallel$times, 3), collapse = ", "),
        "  |  mean: ", round(mean(benchmark_parallel$times), 3), "s\n", sep = "")
  }

  cat("Benchmarking PLINK pairwise_ld\n")
  benchmark_plink = time_reps(n_reps, function() plink_pairwise_ld(
    prefix       = bed_prefix,
    ld_window    = 999999,
    ld_window_kb = 1000000,
    ld_window_r2 = 0
  ))
  cat("  Elapsed (s): ", paste(round(benchmark_plink$times, 3), collapse = ", "),
      "  |  mean: ", round(mean(benchmark_plink$times), 3), "s\n", sep = "")

  make_row = function(method, times, result) {
    data.frame(
      Method        = method,
      Mean_s        = round(mean(times), 3),
      Min_s         = round(min(times),  3),
      Max_s         = round(max(times),  3),
      Output_Pairs  = nrow(result),
      Pairs_per_sec = round(nrow(result) / mean(times)),
      stringsAsFactors = FALSE
    )
  }

  rows = list(make_row("HapSelect (serial)", benchmark_serial$times, benchmark_serial$result))
  if (!is.null(benchmark_parallel)) {
    rows[[length(rows) + 1]] = make_row("HapSelect (parallel)", benchmark_parallel$times, benchmark_parallel$result)
  }
  rows[[length(rows) + 1]] = make_row("PLINK", benchmark_plink$times, benchmark_plink$result)

  summary_df = do.call(rbind, rows)
  row.names(summary_df) = NULL

  cat("\nBenchmark summary (", n_reps, " reps each)\n", sep = "")
  print(summary_df, row.names = FALSE)
  cat("\nExpected pair count: ", expected_pairs, "\n", sep = "")

  list(
    benchmark      = "ld",
    timestamp      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    params         = params,
    results        = lapply(seq_len(nrow(summary_df)), function(i) as.list(summary_df[i, ])),
    expected_pairs = expected_pairs
  )
}

# Only execute when run directly, not when source()'d to load run_benchmark_ld().
if (sys.nframe() == 0L) {
  run_benchmark_ld(parse_args(list(
    n_markers     = 500L,
    n_individuals = 200L,
    n_chr         = 5L,
    missing_rate  = 0.02,
    seed          = 1L,
    n_reps        = 3L
  )))
}
