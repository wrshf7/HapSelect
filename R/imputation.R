########################################
###### Beagle-Based Imputation #########
########################################

# nocov start

##### Impute (and phase) an existing VCF with Beagle #####
# Beagle always phases, and phasing fills in missing genotypes as it goes, so a single run both
# phases and imputes. Its impute= argument only controls markers that are absent from the input
# but present in a reference panel, so it has no effect unless ref= is supplied.
# vcf_in    : path to the input VCF (.vcf or .vcf.gz) containing the genotypes to impute
# out_prefix: output file prefix; Beagle writes <out_prefix>.vcf.gz (and <out_prefix>.log)
# ref       : optional reference panel VCF, for reference-based imputation
# map       : optional PLINK-format genetic map, passed through to Beagle's map= argument
# extra_args: additional Beagle key=value arguments, e.g. c("ne=100", "window=40")
# Returns the path to the imputed output VCF.
beagle_impute = function(vcf_in, out_prefix, ref = NULL, map = NULL, extra_args = character()) {
  report_step("beagle_impute")

  beagle_impute_run(vcf_in, out_prefix, ref = ref, map = map, extra_args = extra_args)
}

##### Run Beagle over a VCF, following its progress #####
# Shared by beagle_impute() and run_beagle_on_geno(). Reports no step itself: each exported
# function reports its own, so a call shows one step whichever input it starts from.
#
# markers: data frame of the input's Chromosome and Position, in file order, which Beagle's
#          progress is measured against. Read from vcf_in when not given; run_beagle_on_geno()
#          passes the genotypes' own, having just written the file from them.
#
# The remaining arguments are documented on beagle_impute().
beagle_impute_run = function(vcf_in, out_prefix, ref = NULL, map = NULL, extra_args = character(),
                             markers = NULL) {
  # Check the beagle input VCF exists
  if (!file.exists(vcf_in)) {
    stop("Beagle input VCF not found: ", vcf_in)
  }

  # Construct the args
  args = c(paste0("gt=", vcf_in), paste0("out=", out_prefix))

  # If reference panel is provided
  if (!is.null(ref)) {
    # If it was not loadable, stop
    if (!file.exists(ref)) {
      stop("Beagle reference panel not found: ", ref)
    }
    args = c(args, paste0("ref=", ref))
  }
  
  # If the plink map was provided
  if (!is.null(map)) {
    # If it was not loadable, stop
    if (!file.exists(map)) {
      stop("Beagle genetic map not found: ", map)
    }
    args = c(args, paste0("map=", map))
  }

  # Run the beagle command, ticking as each of its windows completes
  args = c(args, extra_args)
  if (is.null(markers)) {
    report_status("Scanning markers")
    markers = read_vcf_positions(vcf_in)
  }
  report_status("Running Beagle")
  with_progress({
    p = progressor(steps = nrow(markers))
    run_beagle_command(args, on_line = beagle_progress_parser(markers, p))
  })

  # Return the vcf outpute path
  paste0(out_prefix, ".vcf.gz")
}

##### Run Beagle over a genotype data frame #####
# Shared by beagle_impute_geno() and beagle_phase_geno(), which differ only in how they read the
# result back. Writes geno out as a temporary VCF, runs Beagle over it, and hands the output VCF
# to reader. Reports no step itself, as each caller reports its own.
#
# geno  : data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#         cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
# reader: function taking the output VCF path, e.g. read_vcf_geno or read_vcf_phased
#
# ref, map and extra_args are passed straight through to beagle_impute(), documented there.
run_beagle_on_geno = function(geno, reader, ref = NULL, map = NULL, extra_args = character()) {
  # Check the geno dataframe is valid
  if (!is.data.frame(geno) || ncol(geno) < 4) {
    stop("geno must be a data frame with columns: marker, chromosome, position, and at least one genotype column.")
  }

  # Beagle only handles diploid genotypes, so stop before running it on anything else
  dosages = as.matrix(geno[, -(1:3), drop = FALSE])
  if (any(!(dosages[!is.na(dosages)] %in% c(0, 1, 2)))) {
    stop("Beagle only supports diploid genotypes, so geno dosages must be 0, 1, 2, or NA.")
  }

  # Create a temporary VCF file
  in_vcf = tempfile("hapselect_beagle_in_", fileext = ".vcf")
  out_prefix = tempfile("hapselect_beagle_out_")
  on.exit(unlink(c(in_vcf, paste0(out_prefix, c(".vcf.gz", ".log"))), force = TRUE), add = TRUE)
  report_status("Writing VCF")
  write_vcf_geno(geno, in_vcf)

  # Run Beagle and read the resulting vcf file back in the requested format
  out_vcf = beagle_impute_run(in_vcf, out_prefix, ref = ref, map = map, extra_args = extra_args,
                              markers = data.frame(Chromosome = geno[[2]], Position = geno[[3]]))
  reader(out_vcf)
}

##### Impute missing genotypes in a genotype data frame using Beagle #####
# Returns the same data frame layout as geno (marker, chromosome, position, then one dosage
# column per individual, in the original sample order and names).
#
# Note that this discards the phasing Beagle produced, since a dosage cannot express which
# haplotype an allele sits on. Use beagle_phase_geno() for the haplotype pipeline.
#
# geno      : data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#             cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
# ref       : optional reference panel VCF, for reference-based imputation
# map       : optional PLINK-format genetic map, passed through to Beagle's map= argument
# extra_args: additional Beagle key=value arguments, e.g. c("ne=100", "window=40")
# prefer_ds : TRUE to read Beagle's DS field back instead of hard-calling its GT, keeping the
#             imputation's uncertainty as a fractional dosage. Off by default because it changes
#             the return type: the whole numbers this function otherwise promises become floats,
#             which the LD and blocking functions accept but which PLINK's text format, for one,
#             cannot represent at all.
#
#             Beagle only writes DS when it is imputing against a reference panel. A run without
#             ref= produces a GT-only VCF whatever flags it is given, so prefer_ds has nothing to
#             act on there and warns rather than quietly returning hard calls.
beagle_impute_geno = function(geno, ref = NULL, map = NULL, extra_args = character(),
                              prefer_ds = FALSE) {
  report_step("beagle_impute")

  if (isTRUE(prefer_ds) && is.null(ref)) {
    warning("prefer_ds = TRUE has no effect without a reference panel: Beagle writes the DS ",
            "field only when imputing against ref=, so a gt=-only run produces a GT-only VCF ",
            "and the dosages returned here will be whole-number hard calls.")
  }

  # A closure rather than another argument through run_beagle_on_geno(), whose other caller
  # reads haplotypes and has no dosage field to choose.
  reader = function(path) read_vcf_geno(path, prefer_ds = prefer_ds)
  run_beagle_on_geno(geno, reader, ref = ref, map = map, extra_args = extra_args)
}

##### Phase (and impute) a genotype data frame using Beagle #####
# Same Beagle run as beagle_impute_geno(), but keeps the phasing: returns the geno_phased layout
# of marker, chromosome, position, then <individual>_1 and <individual>_2 columns of 0/1 allele
# presence, which is what compute_haplotype_effects() expects.
#
# geno      : data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#             cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
# ref       : optional reference panel VCF, for reference-based imputation
# map       : optional PLINK-format genetic map, passed through to Beagle's map= argument
# extra_args: additional Beagle key=value arguments, e.g. c("ne=100", "window=40")
beagle_phase_geno = function(geno, ref = NULL, map = NULL, extra_args = character()) {
  report_step("beagle_phase")

  run_beagle_on_geno(geno, read_vcf_phased, ref = ref, map = map, extra_args = extra_args)
}
# nocov end


##### Follow Beagle's progress from its output #####
# Beagle works through the genome in overlapping windows, announcing each one with a line like
#
#   Window 3 [2:257460-40222928]
#
# naming its chromosome and position range. It never says how many windows there will be, but
# the markers are known, so progress is measured in markers: once the next window is announced,
# or Beagle prints its closing statistics, every marker up to the end of the previous window is
# done. Every other line is ignored.
#
# markers: data frame of Chromosome and Position, in the order Beagle reads them
# p      : progressor with one step per marker
# Returns a function taking one line of output, for run_beagle_command()'s on_line.
beagle_progress_parser = function(markers, p) {
  chrom = as.character(markers[[1]])
  pos = as.numeric(markers[[2]])
  # Row indices of each chromosome's markers, in file order, which Beagle requires be by position
  by_chrom = split(seq_along(chrom), factor(chrom, levels = unique(chrom)))

  done = 0L
  pending_end = NULL

  # Move the progressor to the last marker at or before a window's end
  advance_to_window_end = function(window_end) {
    rows = by_chrom[[window_end$chrom]]
    if (is.null(rows)) return(invisible(NULL))
    k = findInterval(window_end$pos, pos[rows])
    if (k > 0) advance(rows[[k]])
  }

  # Ticks forward only, since overlapping windows step back over markers already counted
  advance = function(to) {
    if (to > done) {
      p(amount = to - done)
      done <<- to
    }
  }

  function(line) {
    # The chromosome is matched greedily, so a name containing a colon still parses
    window = regmatches(line, regexec("^Window [0-9]+ \\[(.+):([0-9]+)-([0-9]+)\\]", line))[[1]]
    if (length(window) == 4) {
      if (!is.null(pending_end)) advance_to_window_end(pending_end)
      pending_end <<- list(chrom = window[[2]], pos = as.numeric(window[[4]]))
    } else if (startsWith(line, "Cumulative Statistics")) {
      advance(length(chrom))
      pending_end <<- NULL
    }
    invisible(NULL)
  }
}
