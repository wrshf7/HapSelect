# Only load the package and support helpers when run directly (Rscript benchmark_ld.R).
# When source()'d by benchmark_batched.R, sys.nframe() > 0 and these are already set up.
if (sys.nframe() == 0L) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload to run this benchmark from the HapSelect source tree.")
  }
  pkgload::load_all(".", quiet = TRUE)
  source(file.path("inst", "scripts", "benchmarks", "benchmark_support.R"))
}

# Times the R and compiled pairwise LD implementations, serially and in parallel, against
# PLINK, on the same simulated genotypes. Each method runs in two modes: every
# within-chromosome pair, and a forward window of markers. Results are only checked for
# sanity, not compared across methods; tests/testthat/test-plink-vs-hapstack-ld.R covers
# agreement.
#
# Parameters, all overridable as --name=value or from a batched config:
#   window        : forward marker window used by the windowed mode
#   include_r     : time the R pairwise_ld(). It builds a data frame per marker pair, so
#                   turn this off for large marker counts
#   include_plink : time PLINK --r2. Skipped with a message when PLINK is not found
#
# Parallel runs need an installed HapSelect: multisession workers are fresh R processes
# that load the installed package rather than this source tree, so they run whatever
# version was last installed. Reinstall after changing the C++ before trusting them.

# ld_bench_defaults ------------------------------------------------------------
ld_bench_defaults = list(
  n_markers     = 500L,
  n_individuals = 200L,
  n_chr         = 5L,
  missing_rate  = 0.02,
  seed          = 1L,
  n_reps        = 3L,
  window        = 20L,
  include_r     = TRUE,
  include_plink = TRUE
)

# ld_expected_pairs ------------------------------------------------------------
# The number of within-chromosome pairs a run should report, before any pair with an
# undefined r^2 is dropped.
#
# chrom  : chromosome of each marker
# window : forward marker window, or NULL for every pair
ld_expected_pairs = function(chrom, window) {
  sizes = as.numeric(table(chrom))
  sum(vapply(sizes, function(n) {
    if (is.null(window) || window >= n - 1) return(choose(n, 2))
    window * n - window * (window + 1) / 2
  }, numeric(1)))
}

# ld_sanity --------------------------------------------------------------------
# Checks one LD table on its own terms, returning "ok" or the first check it fails.
#
# ld             : LD table as returned by pairwise_ld() or plink_pairwise_ld()
# expected_pairs : the most pairs the run could report
# window         : forward marker window the run used, or NULL for every pair
ld_sanity = function(ld, expected_pairs, window) {
  required = c("Chrom", "Locus1", "Locus2", "Name1", "Name2", "LD")
  missing  = setdiff(required, names(ld))
  if (length(missing) > 0) return(paste("missing", paste(missing, collapse = ", ")))
  if (nrow(ld) == 0) return("no pairs")
  if (nrow(ld) > expected_pairs) return("too many pairs")
  if (anyNA(ld$LD)) return("NA r^2")
  if (any(ld$LD < 0 | ld$LD > 1 + 1e-9)) return("r^2 outside [0, 1]")
  if (!is.null(window) && any(abs(ld$Locus2 - ld$Locus1) > window)) return("pair outside window")
  "ok"
}

run_benchmark_ld = function(params = list()) {
  params = coerce_params(params, ld_bench_defaults)

  n_markers     = params$n_markers
  n_individuals = params$n_individuals
  n_chr         = params$n_chr
  missing_rate  = params$missing_rate
  seed          = params$seed
  n_reps        = params$n_reps

  if (params$window < 1L) stop("window must be a positive marker count.")

  # NULL is pairwise_ld()'s every-pair window
  modes = list("all pairs" = NULL, "windowed" = params$window)

  has_installed_pkg = requireNamespace("HapSelect", quietly = TRUE)
  has_plink = params$include_plink &&
    !inherits(tryCatch(find_plink(), error = function(e) e), "error")

  geno = simulate_genotypes(n_markers, n_individuals, n_chr, missing_rate, seed)
  expected_pairs = vapply(modes, function(w) ld_expected_pairs(geno[[2]], w), numeric(1))

  cat(
    "Markers:      ", n_markers,     "\n",
    "Individuals:  ", n_individuals, "\n",
    "Chromosomes:  ", n_chr,         "\n",
    "Missing rate: ", missing_rate,  "\n",
    "Window:       ", params$window, " markers (windowed mode)\n",
    "Reps:         ", n_reps,        "\n\n",
    sep = ""
  )

  # PLINK reads a binary fileset, written once and shared by both modes. --chr-set, or
  # PLINK reads chromosomes past 22 as human X, Y, XY and MT
  if (params$include_plink && !has_plink) {
    cat("PLINK executable not found, so PLINK is skipped.\n\n")
  } else if (has_plink) {
    work_dir = tempfile("hapselect_ld_benchmark_")
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
    dir.create(work_dir)

    text_prefix = file.path(work_dir, "synthetic_ld")
    bed_prefix  = file.path(work_dir, "synthetic_ld_bin")
    chr_set     = c("--chr-set", as.character(n_chr))

    cat("Preparing PLINK binary fileset\n\n")
    write_plink_text_files(geno, text_prefix)
    run_plink(c("--file", text_prefix, "--make-bed", "--out", bed_prefix, chr_set))
  }

  if (!has_installed_pkg) {
    cat(
      "Parallel runs skipped: workers need an installed HapSelect package (multisession spawns fresh R processes).\n",
      "Run `Rscript -e \"install.packages('.', repos = NULL, type = 'source')\"` from the repo root,\n",
      "then rerun this benchmark.\n\n",
      sep = ""
    )
  }

  # One entry per method: implementation, execution, and a function of the window that
  # computes the LD table. Parallel runs sit behind the installed-package check
  methods = list()
  add_method = function(implementation, execution, fn) {
    methods[[length(methods) + 1L]] <<- list(implementation = implementation,
                                             execution = execution, fn = fn)
  }

  # Each parallel worker attaches HapSelect and prints its startup banner, once per rep,
  # partly as messages and partly to stdout
  quietly = function(expr) {
    utils::capture.output(result <- suppressMessages(expr))
    result
  }

  if (params$include_r) {
    add_method("R", "serial", function(w) pairwise_ld(geno, parallelize = FALSE, window = w))
    if (has_installed_pkg) {
      add_method("R", "parallel", function(w) quietly(pairwise_ld(geno, parallelize = TRUE, window = w)))
    }
  }
  add_method("C++", "serial", function(w) pairwise_ld_c(geno, parallelize = FALSE, window = w))
  if (has_installed_pkg) {
    add_method("C++", "parallel", function(w) quietly(pairwise_ld_c(geno, parallelize = TRUE, window = w)))
  }
  if (has_plink) {
    # PLINK counts the index marker in its window, so w markers apart is --ld-window w + 1
    add_method("PLINK", "serial", function(w) plink_pairwise_ld(
      prefix       = bed_prefix,
      ld_window    = if (is.null(w)) 999999 else w + 1L,
      ld_window_kb = 1000000,
      ld_window_r2 = 0,
      extra_args   = chr_set
    ))
  }

  rows = list()
  for (mode in names(modes)) {
    window = modes[[mode]]
    baseline_mean = NULL

    for (m in methods) {
      cat("Benchmarking ", mode, ": ", m$implementation, " (", m$execution, ")\n", sep = "")

      benchmark = tryCatch(
        time_reps(n_reps, function() m$fn(window)),
        error = function(e) {
          cat("  Failed: ", conditionMessage(e), "\n", sep = "")
          NULL
        }
      )
      if (is.null(benchmark)) next

      cat("  Elapsed (s): ", paste(round(benchmark$times, 3), collapse = ", "),
          "  |  mean: ", round(mean(benchmark$times), 3), "s\n", sep = "")

      # The first method timed in a mode is the baseline: R serial, or C++ serial when
      # the R implementation is skipped
      mean_s = mean(benchmark$times)
      if (is.null(baseline_mean)) baseline_mean = mean_s

      rows[[length(rows) + 1L]] = data.frame(
        Mode           = mode,
        Implementation = m$implementation,
        Execution      = m$execution,
        Mean_s         = round(mean_s, 3),
        Min_s          = round(min(benchmark$times), 3),
        Max_s          = round(max(benchmark$times), 3),
        Pairs          = nrow(benchmark$result),
        Pairs_per_sec  = round(nrow(benchmark$result) / mean_s),
        Speedup        = round(baseline_mean / mean_s, 2),
        Sane           = ld_sanity(benchmark$result, expected_pairs[[mode]], window),
        stringsAsFactors = FALSE
      )
    }
    cat("\n")
  }

  summary_df = do.call(rbind, rows)
  row.names(summary_df) = NULL

  notes = c(
    paste0("Expected pairs: ",
           paste0(names(expected_pairs), " ", format(expected_pairs, big.mark = ","),
                  collapse = ", "),
           " (pairs with an undefined r^2 are dropped, so a run may report fewer)"),
    paste0("Speedup is relative to the first method in each mode: ",
           if (params$include_r) "R" else "C++", " serial")
  )
  if (has_plink) {
    notes = c(notes, "PLINK times include reading the binary fileset and parsing its .ld output")
  }

  print_benchmark_table(
    summary_df,
    title = paste0("Benchmark summary (", n_reps, " reps each)"),
    group = summary_df$Mode,
    notes = notes
  )

  list(
    benchmark      = "ld",
    timestamp      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    params         = params,
    results        = lapply(seq_len(nrow(summary_df)), function(i) as.list(summary_df[i, ])),
    expected_pairs = as.list(expected_pairs)
  )
}

# Only execute when run directly, not when source()'d to load run_benchmark_ld().
if (sys.nframe() == 0L) {
  # invisible(): the summary is already printed, and the returned list is for
  # benchmark_batched.R to serialise, not for reading in the terminal
  invisible(run_benchmark_ld(parse_args(ld_bench_defaults)))
}
