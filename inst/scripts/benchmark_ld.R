parse_args = function(args){
  defaults = list(
    n_markers = 500L,
    n_individuals = 200L,
    n_chr = 5L,
    missing_rate = 0.02,
    seed = 1L,
    work_dir = NULL
  )

  if(length(args) == 0){
    return(defaults)
  }

  parsed = lapply(args, function(arg){
    if(!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)){
      stop("Arguments must use the form --name=value")
    }

    pieces = strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    setNames(list(pieces[2]), pieces[1])
  })

  parsed = do.call(c, parsed)
  unknown = setdiff(names(parsed), names(defaults))
  if(length(unknown) > 0){
    stop("Unknown argument: ", paste(unknown, collapse = ", "))
  }

  defaults[names(parsed)] = parsed
  defaults$n_markers = as.integer(defaults$n_markers)
  defaults$n_individuals = as.integer(defaults$n_individuals)
  defaults$n_chr = as.integer(defaults$n_chr)
  defaults$missing_rate = as.numeric(defaults$missing_rate)
  defaults$seed = as.integer(defaults$seed)
  defaults
}

simulate_genotypes = function(n_markers, n_individuals, n_chr, missing_rate, seed){
  set.seed(seed)

  chrom_sizes = rep(floor(n_markers / n_chr), n_chr)
  chrom_sizes[seq_len(n_markers %% n_chr)] = chrom_sizes[seq_len(n_markers %% n_chr)] + 1L

  geno = matrix(
    sample(c(0L, 1L, 2L), n_markers * n_individuals, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
    nrow = n_markers,
    ncol = n_individuals
  )

  if(missing_rate > 0){
    geno[sample(length(geno), floor(length(geno) * missing_rate))] = NA_integer_
  }

  out = data.frame(
    marker = paste0("m", seq_len(n_markers)),
    chrom = rep(seq_len(n_chr), times = chrom_sizes),
    pos = ave(seq_len(n_markers), rep(seq_len(n_chr), times = chrom_sizes), FUN = seq_along) * 100L,
    geno,
    stringsAsFactors = FALSE
  )

  colnames(out)[4:ncol(out)] = paste0("ind", seq_len(n_individuals))
  out
}


write_plink_text_files = function(geno, prefix){
  utils::write.table(
    data.frame(CHR = geno[[2]], SNP = geno[[1]], CM = 0, BP = geno[[3]], stringsAsFactors = FALSE),
    file = paste0(prefix, ".map"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

  dosage_to_calls = function(x){
    a1 = ifelse(is.na(x), "0", ifelse(x == 2L, "G", "A"))
    a2 = ifelse(is.na(x), "0", ifelse(x == 0L, "A", "G"))
    c(rbind(a1, a2))
  }

  geno_matrix = as.matrix(geno[, -(1:3)])
  ped = do.call(rbind, lapply(seq_len(ncol(geno_matrix)), function(i){
    c(paste0("F", i), paste0("I", i), "0", "0", "0", "-9", dosage_to_calls(geno_matrix[, i]))
  }))

  utils::write.table(
    ped,
    file = paste0(prefix, ".ped"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}


run_plink = function(args){
  status = attr(system2("plink", args = args, stdout = TRUE, stderr = TRUE), "status")
  if(!is.null(status) && status != 0){
    stop("PLINK command failed.")
  }
  invisible(NULL)
}


benchmark_ld = function(n_markers = 500L, n_individuals = 200L, n_chr = 5L,
                        missing_rate = 0.02, seed = 1L, work_dir = NULL){
  if(!requireNamespace("pkgload", quietly = TRUE)){
    stop("Install the pkgload package to run this benchmark from the FastStack source tree.")
  }
  pkgload::load_all(".", quiet = TRUE)

  cleanup = is.null(work_dir)
  if(cleanup){
    work_dir = tempfile("faststack_ld_benchmark_")
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  geno = simulate_genotypes(n_markers, n_individuals, n_chr, missing_rate, seed)
  text_prefix = file.path(work_dir, "synthetic_ld")
  bed_prefix = file.path(work_dir, "synthetic_ld_bin")
  expected_pairs = sum(vapply(split(geno[[1]], geno[[2]]), function(x) choose(length(x), 2), numeric(1)))

  cat(
    "Generating synthetic genotype data\n",
    "Markers: ", n_markers, "\n",
    "Individuals: ", n_individuals, "\n",
    "Chromosomes: ", n_chr, "\n",
    "Missing rate: ", missing_rate, "\n",
    "Work directory: ", work_dir, "\n\n",
    sep = ""
  )

  cat("Writing PLINK PED/MAP files\n")
  write_plink_text_files(geno, text_prefix)

  cat("Converting PED/MAP to PLINK binary fileset\n")
  run_plink(c("--file", text_prefix, "--make-bed", "--out", bed_prefix))

  cat("Running FastStack pairwise_ld\n")
  faststack_time = system.time({
    faststack_ld = pairwise_ld(geno, parallelize = FALSE)
  })

  cat("Running PLINK-backed LD calculation\n")
  plink_time = system.time({
    plink_ld = plink_pairwise_ld(
      prefix = bed_prefix,
      ld_window = 999999,
      ld_window_kb = 1000000,
      ld_window_r2 = 0
    )
  })

  summary_df = data.frame(
    Method = c("FastStack pairwise_ld", "PLINK pairwise_ld"),
    Elapsed = c(faststack_time[["elapsed"]], plink_time[["elapsed"]]),
    Output_Rows = c(nrow(faststack_ld), nrow(plink_ld)),
    stringsAsFactors = FALSE
  )

  cat("\nBenchmark summary\n")
  print(summary_df, row.names = FALSE)
  cat(
    "\nExpected pair count: ", expected_pairs, "\n",
    "FastStack pair count: ", nrow(faststack_ld), "\n",
    "PLINK pair count: ", nrow(plink_ld), "\n",
    sep = ""
  )

  invisible(summary_df)
}


main = function(){
  do.call(benchmark_ld, parse_args(commandArgs(trailingOnly = TRUE)))
}


if(sys.nframe() == 0){
  main()
}
