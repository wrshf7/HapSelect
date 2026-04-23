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
  status = attr(call_plink(args, stdout = TRUE, stderr = TRUE), "status")
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
