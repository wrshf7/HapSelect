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

##### Read a (optionally gzipped) VCF back into a HapSelect genotype data frame #####
# Converts the GT (genotype) field of every sample column back to a dosage (0 / 1 / 2 / NA), counting
# ALT alleles; phased ("|") and unphased ("/") genotypes are both accepted.
read_vcf_geno = function(path) {
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

  # If there are no records, return an empty geno data frame with the right columns
  if (length(data_lines) == 0) {
    result = data.frame(Marker = character(), Chrom = integer(), Position = numeric(),
                         stringsAsFactors = FALSE, check.names = FALSE)
    result[sample_cols] = numeric()
    return(result)
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

  # Convert a vector of GT calls (phased or unphased) to dosages, ./. becomes NA
  gt_to_dosage = function(gt) {
    alleles = strsplit(gsub("|", "/", gt, fixed = TRUE), "/", fixed = TRUE)
    vapply(alleles, function(a) if (any(a == ".")) NA_real_ else sum(as.integer(a)), numeric(1))
  }

  # Pull out and convert the GT field for every sample column
  dosage_cols = lapply(sample_cols, function(s) {
    sample_fields = strsplit(tab[, s], ":", fixed = TRUE)
    gt = vapply(seq_along(sample_fields), function(i) sample_fields[[i]][gt_index[i]], character(1))
    gt_to_dosage(gt)
  })
  names(dosage_cols) = sample_cols

  # Assemble the final geno-shaped data frame
  result = data.frame(
    Marker = tab[, "ID"],
    Chrom = utils::type.convert(tab[, "CHROM"], as.is = TRUE),
    Position = as.numeric(tab[, "POS"]),
    dosage_cols,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  row.names(result) = NULL

  return(result)
}

########################################
##### Beagle Imputation Functions #####
########################################

# nocov start

##### Impute (and phase) an existing VCF with Beagle #####
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

##### Impute missing genotypes in a genotype data frame using Beagle #####
# Writes geno out as a temporary VCF, imputes it with Beagle, and reads the result back into
# the same data frame layout (marker, chromosome, position, then one dosage column per
# individual, in the original sample order and names).
#
# geno: data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#       cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
beagle_impute_geno = function(geno, ref = NULL, map = NULL, extra_args = character()) {
  # Check the geno dataframe is valid
  if (!is.data.frame(geno) || ncol(geno) < 4) {
    stop("geno must be a data frame with columns: marker, chromosome, position, and at least one genotype column.")
  }

  # Create a temporary VCF file
  in_vcf = tempfile("hapselect_beagle_in_", fileext = ".vcf")
  out_prefix = tempfile("hapselect_beagle_out_")
  on.exit(unlink(c(in_vcf, paste0(out_prefix, c(".vcf.gz", ".log"))), force = TRUE), add = TRUE)
  write_vcf_geno(geno, in_vcf)

  # Perform imputation and read the resulting vcf file
  out_vcf = beagle_impute(in_vcf, out_prefix, ref = ref, map = map, extra_args = extra_args)
  read_vcf_geno(out_vcf)
}
# nocov end
