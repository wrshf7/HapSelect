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

# Parse --name=value command-line arguments, coercing to match types in defaults.
parse_args = function(defaults) {
  args   = commandArgs(trailingOnly = TRUE)
  values = defaults

  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Arguments must be provided as --name=value")
    }

    key   = sub("^--([^=]+)=.*$", "\\1", arg)
    value = sub("^--[^=]+=(.*)$", "\\1", arg)

    if (!key %in% names(values)) {
      stop("Unknown argument: ", key)
    }

    if (is.integer(values[[key]])) {
      values[[key]] = as.integer(value)
    } else if (is.numeric(values[[key]])) {
      values[[key]] = as.numeric(value)
    } else if (is.logical(values[[key]])) {
      values[[key]] = as.logical(value)
    } else {
      values[[key]] = value
    }
  }

  values
}

# Merge params over defaults and coerce each value to match the default's type.
coerce_params = function(params, defaults) {
  params = modifyList(defaults, params)
  for (nm in names(defaults)) {
    if (is.integer(defaults[[nm]])) {
      params[[nm]] = as.integer(params[[nm]])
    } else if (is.numeric(defaults[[nm]])) {
      params[[nm]] = as.numeric(params[[nm]])
    } else if (is.logical(defaults[[nm]])) {
      params[[nm]] = as.logical(params[[nm]])
    }
  }
  params
}

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

  chrom_markers = split(geno[, 1], geno[, 2])
  chrom_names   = names(chrom_markers)
  chrom_sizes   = lengths(chrom_markers)

  base_blocks = rep(floor(n_blocks / length(chrom_markers)), length(chrom_markers))
  if (sum(base_blocks) == 0L) {
    base_blocks[] = 1L
  }

  remainder = n_blocks - sum(base_blocks)
  if (remainder > 0L) {
    order_idx = order(chrom_sizes, decreasing = TRUE)
    base_blocks[order_idx[seq_len(remainder)]] = base_blocks[order_idx[seq_len(remainder)]] + 1L
  }
  base_blocks = pmax(base_blocks, 1L)
  names(base_blocks) = chrom_names

  block_rows = purrr::imap_dfr(chrom_markers, function(markers, chrom_name) {
    n_chr_blocks = min(length(markers), base_blocks[[chrom_name]])
    block_assign = sample(rep(seq_len(n_chr_blocks), length.out = length(markers)))
    block_strings = tapply(markers, block_assign, function(m) paste(m, collapse = ";"))

    data.frame(
      Chrom    = geno[match(markers[1], geno[, 1]), 2],
      Block_ID = paste0("Chr", chrom_name, ":B", seq_along(block_strings)),
      Block    = as.character(block_strings),
      stringsAsFactors = FALSE
    )
  })

  block_rows
}
