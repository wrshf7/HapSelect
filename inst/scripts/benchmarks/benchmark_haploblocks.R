# Only load the package and support helpers when run directly (Rscript benchmark_haploblocks.R).
# When source()'d by benchmark_batched.R, sys.nframe() > 0 and these are already set up.
if (sys.nframe() == 0L) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload to run this benchmark from the HapSelect source tree.")
  }
  pkgload::load_all(".", quiet = TRUE)
  source(file.path("inst", "scripts", "benchmarks", "benchmark_support.R"))
}

run_benchmark_haploblocks = function(params = list()) {
  params = coerce_params(params, list(
    threshold = 0.2,
    tolerance = 4L,
    tol_reset = TRUE,
    n_reps    = 3L
  ))

  threshold = params$threshold
  tolerance = params$tolerance
  tol_reset = params$tol_reset
  n_reps    = params$n_reps

  e = new.env(parent = emptyenv())
  load(file.path("data", "pairwise_ld.rda"), envir = e)
  load(file.path("data", "map.rda"),         envir = e)
  ld  = e$ld_pairs
  map = e$map

  configs = list(
    list(method = "flanking", start = "LD"),
    list(method = "flanking", start = "beginning"),
    list(method = "average",  start = "LD"),
    list(method = "average",  start = "beginning")
  )

  cat(
    "SNPs:        ", length(unique(c(ld$Name1, ld$Name2))), "\n",
    "LD pairs:    ", nrow(ld), "\n",
    "Chromosomes: ", length(unique(ld$Chrom)), "\n",
    "Threshold:   ", threshold, "  Tolerance: ", tolerance,
    "  Tol reset: ", tol_reset, "\n",
    "Reps:        ", n_reps, "\n\n",
    sep = ""
  )

  results = lapply(configs, function(cfg) {
    label = paste0("method=", cfg$method, "  start=", cfg$start)
    cat("Benchmarking: ", label, "\n", sep = "")

    benchmark = time_reps(n_reps, function() def_blocks(
      ld        = ld,
      map       = map,
      method    = cfg$method,
      threshold = threshold,
      tolerance = tolerance,
      tol_reset = tol_reset,
      start     = cfg$start,
      parallel  = FALSE
    ))
    n_blocks = sum(lengths(benchmark$result))

    cat("  Elapsed (s): ", paste(round(benchmark$times, 3), collapse = ", "),
        "  |  mean: ", round(mean(benchmark$times), 3), "s\n", sep = "")

    list(
      method         = cfg$method,
      start          = cfg$start,
      mean_s         = round(mean(benchmark$times), 3),
      min_s          = round(min(benchmark$times),  3),
      max_s          = round(max(benchmark$times),  3),
      total_blocks   = n_blocks,
      blocks_per_sec = round(n_blocks / mean(benchmark$times))
    )
  })

  summary_df = do.call(rbind, lapply(results, function(r) {
    data.frame(
      Method         = r$method,
      Start          = r$start,
      Mean_s         = r$mean_s,
      Min_s          = r$min_s,
      Max_s          = r$max_s,
      Total_Blocks   = r$total_blocks,
      Blocks_per_sec = r$blocks_per_sec,
      stringsAsFactors = FALSE
    )
  }))
  row.names(summary_df) = NULL

  cat("\nBenchmark summary (", n_reps, " reps each)\n", sep = "")
  print(summary_df, row.names = FALSE)

  list(
    benchmark = "haploblocks",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    params    = params,
    results   = results
  )
}

# Only execute when run directly, not when source()'d to load run_benchmark_haploblocks().
if (sys.nframe() == 0L) {
  run_benchmark_haploblocks(parse_args(list(
    threshold = 0.2,
    tolerance = 4L,
    tol_reset = TRUE,
    n_reps    = 3L
  )))
}
