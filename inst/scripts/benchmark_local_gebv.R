if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload to run this benchmark from the FastStack source tree.")
}
pkgload::load_all(".", quiet = TRUE)
source(file.path("inst", "scripts", "benchmark_support.R"))

# Settings ---------------------------------------------------------------------
n_markers     = 500L
n_individuals = 200L
n_chr         = 5L
n_blocks      = 50L   # number of haploblocks to simulate
missing_rate  = 0.02
seed          = 1L
n_reps        = 3L

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
    markers = geno[, 1]
    # assign each marker to a random block (keep within chromosome for realism)
    block_assign = sample(rep(seq_len(n_blocks), length.out = length(markers)))

    blocks = tapply(markers, block_assign, function(m) paste(m, collapse = ";"))

    data.frame(
        Block_ID = paste0("B", names(blocks)),
        Block    = as.character(blocks),
        stringsAsFactors = FALSE
    )
}

make_row = function(method, times, n_blocks) {
  data.frame(
    Method         = method,
    Mean_s         = round(mean(times), 3),
    Min_s          = round(min(times),  3),
    Max_s          = round(max(times),  3),
    Blocks         = n_blocks,
    Blocks_per_sec = round(n_blocks / mean(times))
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
  "Blocks:       ", n_blocks,      "\n",
  "Missing rate: ", missing_rate,  "\n",
  "Reps:         ", n_reps,        "\n\n",
  sep = ""
)

# Run --------------------------------------------------------------------------
progressr::handlers("void")

cat("Benchmarking compute_local_GEBV\n")
r_gebv = time_reps(n_reps, function() {
  compute_local_GEBV(
    geno           = geno,
    marker_effects = marker_effects,
    haploblocks_df = haploblocks_df,
    set_missing_NA = TRUE,
    center         = TRUE
  )
})
cat("  Elapsed (s): ", paste(round(r_gebv$times, 3), collapse = ", "),
    "  |  mean: ", round(mean(r_gebv$times), 3), "s\n", sep = "")

# Summary ----------------------------------------------------------------------
summary_df = make_row("LocalGEBV (Serial)", r_gebv$times, n_blocks)
row.names(summary_df) = NULL

cat("\nBenchmark summary (", n_reps, " reps each)\n", sep = "")
print(summary_df, row.names = FALSE)