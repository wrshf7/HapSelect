########################################
###### Beagle-Based Imputation #########
########################################

##### Write a genotype data frame out as a minimal VCF for Beagle #####
# geno: data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#       cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
# Dosages are written as unphased calls (0/0, 0/1, 1/1) with missing values as ./.,
# using placeholder REF/ALT alleles (A/G) since dosage alone does not carry allele identity.
write_vcf_geno = function(geno, path) {
  # Use the geno column names as VCF sample IDs, falling back to generated names if missing
  sample_names = colnames(geno)[-(1:3)]
  if (is.null(sample_names) || any(!nzchar(sample_names))) {
    sample_names = paste0("S", seq_len(ncol(geno) - 3))
  }

  # Extract the dosage matrix, stripping out marker/chrom/position
  geno_matrix = as.matrix(geno[, -(1:3), drop = FALSE])

  # Beagle is built for diploid genomes, so dosages must be 0, 1, 2, or NA.
  # Check every dosage is 0, 1, 2, or NA before encoding it as a genotype call
  if (any(!(geno_matrix[!is.na(geno_matrix)] %in% c(0, 1, 2)))) {
    stop("Genotype dosages must be 0, 1, 2, or NA.")
  }

  # Look up the unphased VCF genotype call for each dosage, missing values become ./.
  gt_calls = c("0" = "0/0", "1" = "0/1", "2" = "1/1")
  gt = matrix(gt_calls[as.character(geno_matrix)], nrow = nrow(geno_matrix))
  gt[is.na(geno_matrix)] = "./."

  # Build the VCF body: one row per marker, with placeholder REF/ALT alleles
  body = data.frame(
    geno[[2]], geno[[3]], geno[[1]], "A", "G", ".", "PASS", ".", "GT", gt,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  colnames(body) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_names)

  # Open the output file for writing
  con = file(path, "w")
  on.exit(close(con))

  # Write the minimal VCF header Beagle expects
  writeLines(
    c("##fileformat=VCFv4.2",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"),
    con
  )

  # Write the body rows, tab-delimited, with the #CHROM line as the column header
  utils::write.table(body, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

##### Every diploid biallelic GT call mapped to its ALT-allele dosage #####
# The decoding half of the gt_calls table in write_vcf_geno(). Beagle returns phased calls while
# write_vcf_geno() writes unphased ones, so each genotype appears under both separators.
dosage_lookup = c("0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2,
                  "0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)

# Missing genotype calls, which become NA rather than an error
missing_calls = c("./.", ".|.")

##### Convert a vector of VCF GT calls to ALT-allele dosages #####
# gt: GT field values, one per sample, e.g. c("0/0", "0|1", "./.")
# Missing calls become NA. Anything else the lookup does not cover - a multi-allelic call such as
# 1/2, or a non-diploid one such as 0/0/1/1 - is an error rather than a silent NA, since a
# bagle only supports a diploid/biallelic genotype.
gt_to_dosage = function(gt) {
  dosage = dosage_lookup[gt]

  # A missing call is expected and becomes NA; any other unmapped call is not
  unknown = is.na(dosage) & !(gt %in% missing_calls)
  if (any(unknown)) {
    stop("Unrecognised GT call(s): ", paste(unique(gt[unknown]), collapse = ", "),
         ". HapSelect dosages require diploid biallelic genotypes.")
  }

  unname(dosage)
}

##### Convert a vector of phased VCF GT calls to haplotype allele columns #####
# gt: phased GT values, one per sample, e.g. c("0|0", "0|1")
# Returns an integer matrix of 0/1 allele presence with one column per haplotype. Beagle output
# is always phased and complete, so an unphased ("0/1") or missing ("./.") call means something
# upstream went wrong and is an error rather than an NA.
gt_to_haplotypes = function(gt) {
  alleles = strsplit(gt, "|", fixed = TRUE)
  if (any(lengths(alleles) != 2)) {
    stop("Expected phased diploid GT calls separated by \"|\", e.g. 0|1.")
  }
  matrix(as.integer(unlist(alleles)), ncol = 2, byrow = TRUE)
}

##### Parse a (optionally gzipped) VCF into map columns and per-sample GT strings #####
# The shared front half of read_vcf_geno() and read_vcf_phased(): everything up to and including
# pulling the GT field out of each sample column, which is the same whichever format is wanted.
# Returns a list of:
#   map - data frame of SNP / Chromosome / Position
#   gt  - named list of GT character vectors, one per sample, in VCF column order
read_vcf_gt = function(path) {
  # Check the file exists
  if (!file.exists(path)) {
    stop("VCF file not found: ", path)
  }

  # Read the VCF file
  con = gzfile(path, "rt")
  lines = readLines(con)
  close(con)

  # Find the #CHROM header line, it names every column including the samples
  header_idx = which(startsWith(lines, "#CHROM"))
  if (length(header_idx) != 1) {
    stop("VCF file must contain exactly one #CHROM header line: ", path)
  }

  # Extract the column names, then drop the header/meta lines to leave just the data rows
  col_names = strsplit(sub("^#", "", lines[header_idx]), "\t")[[1]]
  data_lines = lines[-seq_len(header_idx)]
  data_lines = data_lines[nzchar(data_lines)]

  # Everything after the fixed VCF columns is a sample
  sample_cols = setdiff(col_names, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))

  # If there are no records, hand back an empty map and empty GT vectors
  if (length(data_lines) == 0) {
    return(list(
      map = data.frame(SNP = character(), Chromosome = integer(), Position = numeric(),
                       stringsAsFactors = FALSE, check.names = FALSE),
      gt = stats::setNames(rep(list(character()), length(sample_cols)), sample_cols)
    ))
  }

  # Split every line into its tab-delimited fields
  tab = do.call(rbind, strsplit(data_lines, "\t", fixed = TRUE))
  colnames(tab) = col_names

  # Work out which position in FORMAT holds GT, per record
  format_fields = strsplit(tab[, "FORMAT"], ":", fixed = TRUE)
  gt_index = vapply(format_fields, function(f) match("GT", f), integer(1))
  if (any(is.na(gt_index))) {
    stop("VCF FORMAT field does not include GT for one or more records: ", path)
  }

  # Pull the GT field out of every sample column
  gt = lapply(sample_cols, function(s) {
    sample_fields = strsplit(tab[, s], ":", fixed = TRUE)
    vapply(seq_along(sample_fields), function(i) sample_fields[[i]][gt_index[i]], character(1))
  })
  names(gt) = sample_cols

  list(
    map = data.frame(
      SNP = tab[, "ID"],
      Chromosome = utils::type.convert(tab[, "CHROM"], as.is = TRUE),
      Position = as.numeric(tab[, "POS"]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    gt = gt
  )
}

##### Read a (optionally gzipped) VCF back into a HapSelect genotype data frame #####
# Converts the GT (genotype) field of every sample column back to a dosage (0 / 1 / 2 / NA), counting
# ALT alleles; phased ("|") and unphased ("/") genotypes are both accepted.
# Returns a geno object, so columns are SNP / Chromosome / Position, matching order_map().
# LD and haploblock tables use Chrom instead - do not align this to those.
read_vcf_geno = function(path) {
  parsed = read_vcf_gt(path)

  # One dosage column per sample, in the original sample order and names
  dosage_cols = lapply(parsed$gt, gt_to_dosage)

  result = data.frame(parsed$map, dosage_cols, stringsAsFactors = FALSE, check.names = FALSE)
  row.names(result) = NULL

  return(result)
}

##### Read a (optionally gzipped) VCF back into a HapSelect phased genotype data frame #####
# Splits each phased GT call into its two haplotypes instead of summing them to a dosage, giving
# the geno_phased layout: SNP / Chromosome / Position, then <sample>_1 and <sample>_2 columns of
# 0/1 allele presence. This is the format compute_haplotype_effects() expects.
read_vcf_phased = function(path) {
  parsed = read_vcf_gt(path)

  # Two columns per sample, kept adjacent so the order is S1_1, S1_2, S2_1, S2_2, ...
  hap_cols = list()
  for (s in names(parsed$gt)) {
    haplotypes = gt_to_haplotypes(parsed$gt[[s]])
    hap_cols[[paste0(s, "_1")]] = haplotypes[, 1]
    hap_cols[[paste0(s, "_2")]] = haplotypes[, 2]
  }

  result = data.frame(parsed$map, hap_cols, stringsAsFactors = FALSE, check.names = FALSE)
  row.names(result) = NULL

  return(result)
}

########################################
##### Beagle Imputation Functions #####
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
run_beagle_on_geno = function(geno, reader, ref = NULL, map = NULL, extra_args = character()) {
  # Check the geno dataframe is valid
  if (!is.data.frame(geno) || ncol(geno) < 4) {
    stop("geno must be a data frame with columns: marker, chromosome, position, and at least one genotype column.")
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
# geno: data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#       cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
beagle_impute_geno = function(geno, ref = NULL, map = NULL, extra_args = character()) {
  run_beagle_on_geno(geno, read_vcf_geno, ref = ref, map = map, extra_args = extra_args)
}

##### Phase (and impute) a genotype data frame using Beagle #####
# Same Beagle run as beagle_impute_geno(), but keeps the phasing: returns the geno_phased layout
# of marker, chromosome, position, then <individual>_1 and <individual>_2 columns of 0/1 allele
# presence, which is what compute_haplotype_effects() expects.
#
# geno: data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#       cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
beagle_phase_geno = function(geno, ref = NULL, map = NULL, extra_args = character()) {
  run_beagle_on_geno(geno, read_vcf_phased, ref = ref, map = map, extra_args = extra_args)
}
# nocov end
