if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload to run this benchmark from the FastStack source tree.")
}
pkgload::load_all(".", quiet = TRUE)

# Settings ---------------------------------------------------------------------
n_markers    = 500L
n_individuals = 200L
n_chr        = 5L
missing_rate = 0.02
seed         = 1L
n_reps       = 3L

# Helpers ----------------------------------------------------------------------
simulate_genotypes = function(n_markers, n_individuals, n_chr, missing_rate, seed) {
  set.seed(seed)

  chrom_sizes = rep(floor(n_markers / n_chr), n_chr)
  chrom_sizes[seq_len(n_markers %% n_chr)] = chrom_sizes[seq_len(n_markers %% n_chr)] + 1L

  geno = matrix(
    sample(c(0L, 1L, 2L), n_markers * n_individuals, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
    nrow = n_markers,
    ncol = n_individuals
  )

  if (missing_rate > 0) {
    geno[sample(length(geno), floor(length(geno) * missing_rate))] = NA_integer_
  }

  out = data.frame(
    marker = paste0("m", seq_len(n_markers)),
    chrom  = rep(seq_len(n_chr), times = chrom_sizes),
    pos    = ave(seq_len(n_markers), rep(seq_len(n_chr), times = chrom_sizes), FUN = seq_along) * 100L,
    geno,
    stringsAsFactors = FALSE
  )
  colnames(out)[4:ncol(out)] = paste0("ind", seq_len(n_individuals))
  out
}

write_plink_text_files = function(geno, prefix) {
  utils::write.table(
    data.frame(CHR = geno[[2]], SNP = geno[[1]], CM = 0, BP = geno[[3]], stringsAsFactors = FALSE),
    file      = paste0(prefix, ".map"),
    quote     = FALSE,
    sep       = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

  dosage_to_calls = function(x) {
    a1 = ifelse(is.na(x), "0", ifelse(x == 2L, "G", "A"))
    a2 = ifelse(is.na(x), "0", ifelse(x == 0L, "A", "G"))
    c(rbind(a1, a2))
  }

  geno_matrix = as.matrix(geno[, -(1:3)])
  ped = do.call(rbind, lapply(seq_len(ncol(geno_matrix)), function(i) {
    c(paste0("F", i), paste0("I", i), "0", "0", "0", "-9", dosage_to_calls(geno_matrix[, i]))
  }))

  utils::write.table(
    ped,
    file      = paste0(prefix, ".ped"),
    quote     = FALSE,
    sep       = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}

run_plink = function(args) {
  status = attr(system2("plink", args = args, stdout = TRUE, stderr = TRUE), "status")
  if (!is.null(status) && status != 0) stop("PLINK command failed.")
  invisible(NULL)
}

time_reps = function(n_reps, expr_fn) {
  times = numeric(n_reps)
  result = NULL
  for (i in seq_len(n_reps)) {
    t = system.time({ result = expr_fn() })
    times[i] = t[["elapsed"]]
  }
  list(times = times, result = result)
}

# Data -------------------------------------------------------------------------
geno = simulate_genotypes(n_markers, n_individuals, n_chr, missing_rate, seed)
expected_pairs = sum(vapply(
  split(geno[[1]], geno[[2]]),
  function(x) choose(length(x), 2),
  numeric(1)
))

work_dir = tempfile("faststack_ld_benchmark_")
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

# Run --------------------------------------------------------------------------
cat("Writing PLINK PED/MAP files\n")
write_plink_text_files(geno, text_prefix)

cat("Converting to PLINK binary fileset\n")
run_plink(c("--file", text_prefix, "--make-bed", "--out", bed_prefix))

cat("\nBenchmarking FastStack pairwise_ld (parallelize = FALSE)\n")
r_serial = time_reps(n_reps, function() pairwise_ld(geno, parallelize = FALSE))
cat("  Elapsed (s): ", paste(round(r_serial$times, 3), collapse = ", "),
    "  |  mean: ", round(mean(r_serial$times), 3), "s\n", sep = "")

cat("Benchmarking FastStack pairwise_ld (parallelize = TRUE)\n")
r_parallel = tryCatch(
  time_reps(n_reps, function() pairwise_ld(geno, parallelize = TRUE)),
  error = function(e) {
    cat("  Skipped: parallel workers cannot load an uninstalled package.\n",
        "  Install FastStack with install.packages('.', repos = NULL) and re-run.\n", sep = "")
    NULL
  }
)
if (!is.null(r_parallel)) {
  cat("  Elapsed (s): ", paste(round(r_parallel$times, 3), collapse = ", "),
      "  |  mean: ", round(mean(r_parallel$times), 3), "s\n", sep = "")
}

cat("Benchmarking PLINK pairwise_ld\n")
r_plink = time_reps(n_reps, function() plink_pairwise_ld(
  prefix       = bed_prefix,
  ld_window    = 999999,
  ld_window_kb = 1000000,
  ld_window_r2 = 0
))
cat("  Elapsed (s): ", paste(round(r_plink$times, 3), collapse = ", "),
    "  |  mean: ", round(mean(r_plink$times), 3), "s\n", sep = "")

# Summary ----------------------------------------------------------------------
make_row = function(method, times, result) {
  data.frame(
    Method         = method,
    Mean_s         = round(mean(times), 3),
    Min_s          = round(min(times),  3),
    Max_s          = round(max(times),  3),
    Output_Pairs   = nrow(result),
    Pairs_per_sec  = round(nrow(result) / mean(times)),
    stringsAsFactors = FALSE
  )
}

rows = list(make_row("FastStack (serial)", r_serial$times, r_serial$result))
if (!is.null(r_parallel)) {
  rows[[length(rows) + 1]] = make_row("FastStack (parallel)", r_parallel$times, r_parallel$result)
}
rows[[length(rows) + 1]] = make_row("PLINK", r_plink$times, r_plink$result)

summary_df = do.call(rbind, rows)
row.names(summary_df) = NULL

cat("\nBenchmark summary (", n_reps, " reps each)\n", sep = "")
print(summary_df, row.names = FALSE)
cat("\nExpected pair count: ", expected_pairs, "\n", sep = "")
