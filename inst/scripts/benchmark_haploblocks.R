if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload to run this benchmark from the HapSelect source tree.")
}
pkgload::load_all(".", quiet = TRUE)

# Data -------------------------------------------------------------------------
e = new.env(parent = emptyenv())
load(file.path("data", "pairwise_ld.rda"), envir = e)
load(file.path("data", "map.rda"),         envir = e)
ld  = e$ld_pairs
map = e$map

# Settings ---------------------------------------------------------------------
threshold = 0.2
tolerance = 4L
tol_reset = TRUE
n_reps    = 3L

configs = list(
  list(method = "flanking", start = "LD"),
  list(method = "flanking", start = "beginning"),
  list(method = "average",  start = "LD"),
  list(method = "average",  start = "beginning")
)

# Run --------------------------------------------------------------------------
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

  times    = numeric(n_reps)
  n_blocks = NA_integer_

  for (i in seq_len(n_reps)) {
    t = system.time({
      blocks = def_blocks(
        ld        = ld,
        map       = map,
        method    = cfg$method,
        threshold = threshold,
        tolerance = tolerance,
        tol_reset = tol_reset,
        start     = cfg$start,
        parallel  = FALSE
      )
    })
    times[i] = t[["elapsed"]]
    if (i == 1L) n_blocks = sum(lengths(blocks))
  }

  cat("  Elapsed (s): ", paste(round(times, 3), collapse = ", "),
      "  |  mean: ", round(mean(times), 3), "s\n", sep = "")

  data.frame(
    Method         = cfg$method,
    Start          = cfg$start,
    Mean_s         = round(mean(times), 3),
    Min_s          = round(min(times),  3),
    Max_s          = round(max(times),  3),
    Total_Blocks   = n_blocks,
    Blocks_per_sec = round(n_blocks / mean(times)),
    stringsAsFactors = FALSE
  )
})

summary_df = do.call(rbind, results)
row.names(summary_df) = NULL

cat("\nBenchmark summary (", n_reps, " reps each)\n", sep = "")
print(summary_df, row.names = FALSE)
