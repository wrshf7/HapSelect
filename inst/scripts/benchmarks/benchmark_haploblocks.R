# Only load the package and support helpers when run directly (Rscript benchmark_haploblocks.R).
# When source()'d by benchmark_batched.R, sys.nframe() > 0 and these are already set up.
if (sys.nframe() == 0L) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload to run this benchmark from the HapSelect source tree.")
  }
  pkgload::load_all(".", quiet = TRUE)
  source(file.path("inst", "scripts", "benchmarks", "benchmark_support.R"))
}

# Benchmarks every blocking strategy def_blocks() dispatches on, over the same markers, so
# that runtimes are comparable and the block structures can be read side by side. The
# strategies are not interchangeable: the LD strategy needs a pairwise LD table and the
# window strategy needs only the map, so the inputs each one is given differ even though
# the markers do not. Every strategy runs at its own defaults, so this measures the
# methods as they ship rather than an attempt to equalise them.

# haploblocks_bench_defaults ---------------------------------------------------
# ld_strategy()'s own defaults for threshold, tolerance and tol_reset. window_strategy()
# has no default window, so window_snp is a common marker count and window_map is sized
# to give blocks of a comparable marker count on this dataset, whose median marker
# spacing is about 534 kb.
haploblocks_bench_defaults = list(
  threshold  = 0.7,
  tolerance  = 1L,
  tol_reset  = TRUE,
  window_snp = 10L,
  window_map = 5e6,
  n_reps     = 3L
)

# benchmark_haploblock_data ----------------------------------------------------
# Loads the bundled wheat genotypes and LD table, and derives the map from the genotypes.
#
# The map is derived rather than loaded because data/map.rda is a different dataset
# altogether - a maize map, whose markers do not appear in data/geno.rda or
# data/pairwise_ld.rda at all. perform_ld_blocking() orders each chromosome's markers by
# looking their positions up in the map, so a map that does not contain them leaves the
# order to chance.
benchmark_haploblock_data = function() {
  e = new.env(parent = emptyenv())
  load(file.path("data", "geno.rda"),        envir = e)
  load(file.path("data", "pairwise_ld.rda"), envir = e)

  list(map = order_map(e$geno[, 1:3], verbose = FALSE), ld = e$ld_pairs)
}

run_benchmark_haploblocks = function(params = list()) {
  params = coerce_params(params, haploblocks_bench_defaults)

  n_reps = params$n_reps
  data   = benchmark_haploblock_data()
  map    = data$map
  ld     = data$ld

  # Strategies are built once, outside the timed calls: ld_strategy() copies the LD table
  # into the strategy object, and that cost belongs to neither method's blocking time.
  ld_config = function(method, start) {
    list(strategy = "ld", config = paste0(method, ", start=", start),
         object = ld_strategy(ld, method = method, threshold = params$threshold,
                              tolerance = params$tolerance, tol_reset = params$tol_reset,
                              start = start))
  }
  configs = list(
    ld_config("flanking", "LD"),
    ld_config("flanking", "beginning"),
    ld_config("average",  "LD"),
    ld_config("average",  "beginning"),
    list(strategy = "window", config = paste0("window_snp, ", params$window_snp, " markers"),
         object = window_strategy(params$window_snp, method = "window_snp")),
    list(strategy = "window", config = paste0("window_map, ", params$window_map / 1e6, " Mb"),
         object = window_strategy(params$window_map, method = "window_map"))
  )

  cat(
    "Markers:     ", nrow(map), "\n",
    "LD pairs:    ", nrow(ld), "\n",
    "Chromosomes: ", length(unique(map$Chromosome)), "\n",
    "Reps:        ", n_reps, "\n",
    "\nEach strategy at its own defaults:\n",
    "  ld    : threshold ", params$threshold, ", tolerance ", params$tolerance,
    ", tol_reset ", params$tol_reset, "\n",
    "  window: ", params$window_snp, " markers / ", params$window_map / 1e6, " Mb\n\n",
    sep = ""
  )

  results = lapply(configs, function(cfg) {
    cat("Benchmarking: ", cfg$strategy, " - ", cfg$config, "\n", sep = "")

    benchmark = time_reps(n_reps, function() {
      suppressMessages(def_blocks(map = map, strategy = cfg$object))
    })

    # Block structure, measured once on the last result rather than inside the timing
    block_df = block_obj_to_df(benchmark$result, map)
    summary  = block_summary(block_df)

    cat("  Elapsed (s): ", paste(round(benchmark$times, 3), collapse = ", "),
        "  |  mean: ", round(mean(benchmark$times), 3), "s\n", sep = "")

    list(
      strategy           = cfg$strategy,
      config             = cfg$config,
      mean_s             = round(mean(benchmark$times), 3),
      min_s              = round(min(benchmark$times),  3),
      max_s              = round(max(benchmark$times),  3),
      total_blocks       = nrow(block_df),
      markers_placed     = sum(block_df$Num_SNP),
      mean_snp_per_block = round(summary$Mean_SNP_per_Block, 2),
      max_snp_per_block  = summary$Max_SNP_per_Block,
      pct_singletons     = round(summary$Percent_Singleton_Blocks, 1),
      blocks_per_sec     = round(nrow(block_df) / mean(benchmark$times))
    )
  })

  summary_df = do.call(rbind, lapply(results, function(r) {
    data.frame(
      Strategy       = r$strategy,
      Config         = r$config,
      Mean_s         = r$mean_s,
      Min_s          = r$min_s,
      Max_s          = r$max_s,
      Blocks         = r$total_blocks,
      Markers        = r$markers_placed,
      Mean_SNP       = r$mean_snp_per_block,
      Max_SNP        = r$max_snp_per_block,
      Pct_Singleton  = r$pct_singletons,
      Blocks_per_sec = r$blocks_per_sec,
      stringsAsFactors = FALSE
    )
  }))
  row.names(summary_df) = NULL

  # Both strategies should place every marker, so a shortfall is worth calling out
  dropped = nrow(map) - summary_df$Markers
  notes = vapply(which(dropped > 0), function(i) {
    paste0("Markers not placed by ", summary_df$Strategy[i], " (", summary_df$Config[i],
           "): ", dropped[i], " of ", nrow(map))
  }, character(1))

  print_benchmark_table(
    summary_df,
    title = paste0("Benchmark summary (", n_reps, " reps each, ", nrow(map), " markers)"),
    group = summary_df$Strategy,
    notes = notes
  )

  list(
    benchmark = "haploblocks",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    params    = params,
    results   = results
  )
}

# Only execute when run directly, not when source()'d to load run_benchmark_haploblocks().
if (sys.nframe() == 0L) {
  # invisible(): the summary is already printed, and the returned list is for
  # benchmark_batched.R to serialise, not for reading in the terminal
  invisible(run_benchmark_haploblocks(parse_args(haploblocks_bench_defaults)))
}
