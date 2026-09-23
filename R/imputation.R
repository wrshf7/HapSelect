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

  # Run the beagle command
  args = c(args, extra_args)
  run_beagle_command(args)

  # Return the vcf outpute path
  paste0(out_prefix, ".vcf.gz")
}

##### Run Beagle over a genotype data frame #####
# Shared by beagle_impute_geno() and beagle_phase_geno(), which differ only in how they read the
# result back. Writes geno out as a temporary VCF, runs Beagle over it, and hands the output VCF
# to reader.
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
  write_vcf_geno(geno, in_vcf)

  # Run Beagle and read the resulting vcf file back in the requested format
  out_vcf = beagle_impute(in_vcf, out_prefix, ref = ref, map = map, extra_args = extra_args)
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
beagle_impute_geno = function(geno, ref = NULL, map = NULL, extra_args = character()) {
  run_beagle_on_geno(geno, read_vcf_geno, ref = ref, map = map, extra_args = extra_args)
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
  run_beagle_on_geno(geno, read_vcf_phased, ref = ref, map = map, extra_args = extra_args)
}
# nocov end
