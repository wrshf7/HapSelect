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
# strategies are not interchangeable: the LD strategy needs a pairwise LD table, the graph
# strategy needs the genotypes, and the window strategy needs neither, so the inputs each
# one is given differ even though the markers do not.
#
# Speed alone is misleading here, because the methods do not produce comparable partitions.
# The window and LD strategies place every marker, while the graph strategy discards markers
# that fall inside a block's span without enough LD to join it. The summary therefore reports
# markers placed alongside the timings.
#
# Every strategy runs at its own defaults, so this measures the methods as they ship rather
# than an attempt to equalise them. Do not add a shared r^2 threshold here: ld_strategy()'s
# threshold and graph_strategy()'s theta_core are not equivalent knobs. Measured on this
# dataset, sweeping theta_core alone from 0.8 to 0.9 changes the graph blocking by two
# blocks out of 1700, because theta_extend and theta_bridge stay at 0.20 and reassemble
# whatever a stricter core threshold breaks apart. Setting theta_core to match an LD
# threshold therefore compares a strict setting against a permissive one, which is how an
# earlier version of this script came to report the graph strategy producing smaller blocks
# than it really does. To compare at equal stringency, move theta_core, theta_extend and
# theta_bridge together.

# benchmark_haploblock_data ----------------------------------------------------
# Loads the bundled wheat genotypes and LD table, and derives the map from the genotypes.
#
# The map is derived rather than loaded because data/map.rda is a different dataset
# altogether - a maize map, whose markers do not appear in data/geno.rda or
# data/pairwise_ld.rda at all. Taking the map from the genotypes keeps all three inputs on
# the same markers, which the graph strategy requires and the LD strategy silently needs:
# perform_ld_blocking() orders each chromosome's markers by looking their positions up in
# the map, so a map that does not contain them leaves the order to chance.
benchmark_haploblock_data = function() {
  e = new.env(parent = emptyenv())
  load(file.path("data", "geno.rda"),        envir = e)
  load(file.path("data", "pairwise_ld.rda"), envir = e)

  geno = e$geno
  ld   = e$ld_pairs
  map  = order_map(geno[, 1:3], verbose = FALSE)

  # order_map() transforms the three map columns and leaves the dosage columns behind,
  # so order_geno() applies the same ordering to the genotypes: same markers, same rows,
  # and the map's chromosome numbering, which check_ld_matrix() needs.
  geno = order_geno(geno, map)

  list(geno = geno, map = map, ld = ld)
}

run_benchmark_haploblocks = function(params = list()) {
  # ld_strategy()'s own defaults for threshold, tolerance and tol_reset. window_strategy()
  # has no default window, so window_snp is a common marker count and window_map is sized
  # to give blocks of a comparable marker count on this dataset, whose median marker
  # spacing is about 534 kb. The graph strategy takes no parameters here: it runs at
  # graph_strategy()'s defaults, which are the values the method was validated at.
  params = coerce_params(params, list(
    threshold   = 0.7,
    tolerance   = 1L,
    tol_reset   = TRUE,
    window_snp  = 10L,
    window_map  = 5e6,
    n_reps      = 3L
  ))

  n_reps = params$n_reps
  data   = benchmark_haploblock_data()
  geno   = data$geno
  map    = data$map
  ld     = data$ld

  # Strategies are built once, outside the timed calls: ld_strategy() copies the LD table
  # into the strategy object, and that cost belongs to neither method's blocking time.
  configs = list(
    list(strategy = "ld", config = "flanking, start=LD",
         object = ld_strategy(ld, method = "flanking", threshold = params$threshold,
                              tolerance = params$tolerance, tol_reset = params$tol_reset,
                              start = "LD")),
    list(strategy = "ld", config = "flanking, start=beginning",
         object = ld_strategy(ld, method = "flanking", threshold = params$threshold,
                              tolerance = params$tolerance, tol_reset = params$tol_reset,
                              start = "beginning")),
    list(strategy = "ld", config = "average, start=LD",
         object = ld_strategy(ld, method = "average", threshold = params$threshold,
                              tolerance = params$tolerance, tol_reset = params$tol_reset,
                              start = "LD")),
    list(strategy = "ld", config = "average, start=beginning",
         object = ld_strategy(ld, method = "average", threshold = params$threshold,
                              tolerance = params$tolerance, tol_reset = params$tol_reset,
                              start = "beginning")),
    list(strategy = "window", config = paste0("window_snp, ", params$window_snp, " markers"),
         object = window_strategy(params$window_snp, method = "window_snp")),
    list(strategy = "window", config = paste0("window_map, ", params$window_map / 1e6, " Mb"),
         object = window_strategy(params$window_map, method = "window_map")),
    list(strategy = "graph", config = "defaults",
         object = graph_strategy())
  )

  graph_defaults = graph_strategy()

  cat(
    "Markers:     ", nrow(map), "\n",
    "Individuals: ", ncol(geno) - 3L, "\n",
    "LD pairs:    ", nrow(ld), "\n",
    "Chromosomes: ", length(unique(map$Chromosome)), "\n",
    "Reps:        ", n_reps, "\n",
    "\nEach strategy at its own defaults:\n",
    "  ld    : threshold ", params$threshold, ", tolerance ", params$tolerance,
    ", tol_reset ", params$tol_reset, "\n",
    "  window: ", params$window_snp, " markers / ", params$window_map / 1e6, " Mb\n",
    "  graph : theta_core ", graph_defaults$theta_core,
    ", theta_extend ", graph_defaults$theta_extend,
    ", theta_bridge ", graph_defaults$theta_bridge,
    ", theta_refill ", graph_defaults$theta_refill, "\n",
    "          window_ld ", graph_defaults$window_ld,
    ", max_gap_snps ", graph_defaults$max_gap_snps,
    ", max_gap_markers ", graph_defaults$max_gap_markers,
    ", min_block_snps ", graph_defaults$min_block_snps, "\n\n",
    sep = ""
  )

  results = lapply(configs, function(cfg) {
    cat("Benchmarking: ", cfg$strategy, " - ", cfg$config, "\n", sep = "")

    # Only the graph strategy reads the genotypes; def_blocks() warns if the others are
    # handed them, since they would be ignored
    needs_geno = inherits(cfg$object, "graph_strategy")
    benchmark = time_reps(n_reps, function() {
      suppressMessages(def_blocks(
        strategy = cfg$object,
        map      = map,
        geno     = if (needs_geno) geno else NULL
      ))
    })

    # Block structure, measured once on the last result rather than inside the timing
    block_df = block_obj_to_df(benchmark$result, map)
    summary  = block_summary(block_df)

    cat("  Elapsed (s): ", paste(round(benchmark$times, 3), collapse = ", "),
        "  |  mean: ", round(mean(benchmark$times), 3), "s\n", sep = "")

    list(
      strategy         = cfg$strategy,
      config           = cfg$config,
      mean_s           = round(mean(benchmark$times), 3),
      min_s            = round(min(benchmark$times),  3),
      max_s            = round(max(benchmark$times),  3),
      total_blocks     = nrow(block_df),
      markers_placed   = sum(block_df$Num_SNP),
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

  cat("\nBenchmark summary (", n_reps, " reps each, ", nrow(map), " markers)\n", sep = "")
  print(summary_df, row.names = FALSE)

  dropped = nrow(map) - summary_df$Markers
  if (any(dropped > 0)) {
    cat("\nMarkers not placed in any block:\n")
    for (i in which(dropped > 0)) {
      cat("  ", summary_df$Strategy[i], " (", summary_df$Config[i], "): ",
          dropped[i], " of ", nrow(map), "\n", sep = "")
    }
  }

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
  invisible(run_benchmark_haploblocks(parse_args(list(
    threshold  = 0.7,
    tolerance  = 1L,
    tol_reset  = TRUE,
    window_snp = 10L,
    window_map = 5e6,
    n_reps     = 3L
  ))))
}
