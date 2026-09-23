########################################
########## VCF Read and Write ##########
########################################
#
# Reading and writing VCF files, and converting VCF GT calls to and from the
# HapSelect geno (dosage) and geno_phased (haplotype) layouts. Beagle is the
# reason these exist - see imputation.R - but nothing here is specific to it,
# so any part of the package that needs genotypes in or out of a VCF uses
# these rather than parsing a VCF itself.

##### Write a genotype data frame out as a minimal VCF for Beagle #####
# geno: data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#       cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
# path: file path to write the VCF to
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
# See haplotype_lookup below for the phase-preserving counterpart.
dosage_lookup = c("0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2,
                  "0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)

# Missing genotype calls, which become NA rather than an error. A bare "." is the
# VCF spelling for a genotype that was not called at all, as opposed to "./." for
# a diploid one whose alleles are both unknown; both mean the same thing here.
missing_calls = c("./.", ".|.", ".")

##### Convert a vector of VCF GT calls to ALT-allele dosages #####
# gt: GT field values, one per sample, e.g. c("0/0", "0|1", "./.")
# Missing calls become NA. Anything else the lookup does not cover - a multi-allelic call such as
# 1/2, or a non-diploid one such as 0/0/1/1 - is an error rather than a silent NA, since
# Beagle only supports diploid biallelic genotypes.
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

##### Every phased biallelic GT call mapped to its two haplotype alleles #####
# The haplotype counterpart of dosage_lookup. Only phased ("|") calls appear, since an unphased
# call does not say which allele sits on which haplotype, and no missing form appears, since a
# haplotype column cannot hold an NA the way a dosage can.
haplotype_lookup = list("0|0" = c(0L, 0L), "0|1" = c(0L, 1L),
                        "1|0" = c(1L, 0L), "1|1" = c(1L, 1L))

##### Convert a vector of phased VCF GT calls to haplotype allele columns #####
# gt: phased GT values, one per sample, e.g. c("0|0", "0|1")
# Returns an integer matrix of 0/1 allele presence with one column per haplotype. Beagle output
# is always phased and complete, so an unphased ("0|1" written as "0/1") or missing (".|.") call
# means something upstream went wrong and is an error rather than an NA.
gt_to_haplotypes = function(gt) {
  haplotypes = haplotype_lookup[gt]

  # Anything the table does not cover: unphased, missing, multi-allelic (0|2) or non-diploid
  unknown = vapply(haplotypes, is.null, logical(1))
  if (any(unknown)) {
    stop("Unrecognised GT call(s): ", paste(unique(gt[unknown]), collapse = ", "),
         ". Expected phased diploid GT calls separated by \"|\", e.g. 0|1.")
  }

  matrix(as.integer(unlist(haplotypes)), ncol = 2, byrow = TRUE)
}

##### Parse a (optionally gzipped) VCF into its map, sample columns and FORMAT layout #####
# The shared front half of every reader below: everything up to but not including pulling a
# particular FORMAT field out of the sample columns, which is the same whichever field is wanted.
# path: VCF to read, optionally gzipped
# Returns a list of:
#   map         - data frame of SNP / Chromosome / Position
#   tab         - character matrix of data rows, columns named by the #CHROM header, or NULL
#                 when the file holds no records
#   sample_cols - sample column names, in VCF column order
#   format      - per-record FORMAT field names, as a list of character vectors
read_vcf_records = function(path) {
  # Check the file exists
  if (!file.exists(path)) {
    stop("VCF file not found: ", path)
  }

  # Read the VCF file
  con = gzfile(path, "rt")
  lines = readLines(con)
  close(con)

  # Line endings need no handling here: readLines() accepts LF, CRLF or CR as the
  # terminator whatever mode the connection was opened in, so a VCF written on
  # Windows arrives with its carriage returns already removed.

  # Find the #CHROM header line, it names every column including the samples
  header_idx = which(startsWith(lines, "#CHROM"))
  if (length(header_idx) != 1) {
    stop("VCF file must contain exactly one #CHROM header line: ", path)
  }

  # Extract the column names, then drop the header/meta lines to leave just the data rows
  col_names = strsplit(sub("^#", "", lines[header_idx]), "\t")[[1]]

  # VCF is tab-delimited by spec, so a single column means the file is delimited some other way
  if (length(col_names) < 2) {
    stop("VCF #CHROM header is not tab-delimited, so no columns could be read: ", path,
         "\nVCF requires tab separators between columns. Check the file has not been re-saved ",
         "with spaces or commas, for example by a spreadsheet editor.")
  }

  data_lines = lines[-seq_len(header_idx)]
  data_lines = data_lines[nzchar(data_lines)]

  # Everything after the fixed VCF columns is a sample
  sample_cols = setdiff(col_names, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))

  # If there are no records, hand back an empty map and no rows to pull fields from
  if (length(data_lines) == 0) {
    return(list(
      map = data.frame(SNP = character(), Chromosome = integer(), Position = numeric(),
                       stringsAsFactors = FALSE, check.names = FALSE),
      tab = NULL,
      sample_cols = sample_cols,
      format = list()
    ))
  }

  # Split every line into its tab-delimited fields
  tab = do.call(rbind, strsplit(data_lines, "\t", fixed = TRUE))
  colnames(tab) = col_names

  list(
    map = data.frame(
      SNP = tab[, "ID"],
      Chromosome = utils::type.convert(tab[, "CHROM"], as.is = TRUE),
      Position = as.numeric(tab[, "POS"]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    tab = tab,
    sample_cols = sample_cols,
    format = strsplit(tab[, "FORMAT"], ":", fixed = TRUE),
    path = path
  )
}

##### Pull one FORMAT field out of every sample column of a parsed VCF #####
# FORMAT can differ from record to record, so the field's position is looked up per record rather
# than once for the file.
# records : output of read_vcf_records()
# field   : FORMAT field name, e.g. "GT" or "DS"
# required: TRUE to stop when a record's FORMAT does not list the field, FALSE to return NULL.
#           An optional field is all or nothing - one record missing it means the whole file
#           falls back, since a half-filled column of dosages is not something to guess at.
# Returns a named list of character vectors, one per sample, in VCF column order, or NULL.
read_vcf_field = function(records, field, required = TRUE) {
  # No records means no values to pull, but the samples are still known
  if (is.null(records$tab)) {
    return(stats::setNames(rep(list(character()), length(records$sample_cols)),
                           records$sample_cols))
  }

  # Where the field sits within each record's colon-separated FORMAT
  field_index = vapply(records$format, function(f) match(field, f), integer(1))
  if (any(is.na(field_index))) {
    if (!required) return(NULL)
    stop("VCF FORMAT field does not include ", field, " for one or more records: ", records$path)
  }

  values = lapply(records$sample_cols, function(s) {
    sample_fields = strsplit(records$tab[, s], ":", fixed = TRUE)
    vapply(seq_along(sample_fields), function(i) sample_fields[[i]][field_index[i]], character(1))
  })
  stats::setNames(values, records$sample_cols)
}

##### Parse a (optionally gzipped) VCF into map columns and per-sample GT strings #####
# The GT-field form of read_vcf_records(), kept because GT is what read_vcf_phased() and the
# default read_vcf_geno() path both want.
# path: VCF to read, optionally gzipped
# Returns a list of:
#   map - data frame of SNP / Chromosome / Position
#   gt  - named list of GT character vectors, one per sample, in VCF column order
read_vcf_gt = function(path) {
  records = read_vcf_records(path)
  list(map = records$map, gt = read_vcf_field(records, "GT"))
}

##### Convert a vector of VCF DS values to ALT-allele dosages #####
# DS is the imputed ALT allele dosage: a float on the same 0 to 2 scale as a GT-derived dosage, but
# carrying the imputation's uncertainty rather than rounding it away, so 0.83 stays 0.83 instead of
# becoming the hard call 1. A missing value (".") becomes NA, as it does for GT.
# ds: DS field values, one per sample, e.g. c("0.00", "0.83", ".")
ds_to_dosage = function(ds) {
  dosage = suppressWarnings(as.numeric(ds))

  # A missing value is expected and becomes NA; anything else that will not parse is not
  unknown = is.na(dosage) & !(ds %in% c(".", "") | is.na(ds))
  if (any(unknown)) {
    stop("Unrecognised DS value(s): ", paste(unique(ds[unknown]), collapse = ", "),
         ". HapSelect dosages require numeric ALT allele dosages.")
  }

  # A DS outside 0 to 2 is not a diploid dosage, and usually means the field was misread
  out_of_range = !is.na(dosage) & (dosage < 0 | dosage > 2)
  if (any(out_of_range)) {
    stop("DS value(s) outside the 0 to 2 dosage range: ",
         paste(unique(dosage[out_of_range]), collapse = ", "),
         ". HapSelect dosages require diploid ALT allele dosages.")
  }

  dosage
}

##### Read a (optionally gzipped) VCF back into a HapSelect genotype data frame #####
# Converts the GT (genotype) field of every sample column back to a dosage (0 / 1 / 2 / NA), counting
# ALT alleles; phased ("|") and unphased ("/") genotypes are both accepted.
# Returns a geno object, so columns are SNP / Chromosome / Position, matching order_map().
# LD and haploblock tables use Chrom instead - do not align this to those.
# path     : VCF to read, optionally gzipped
# prefer_ds: TRUE to read dosages from the DS field instead of GT when every record carries one,
#            keeping fractional imputed dosages rather than the rounded hard calls. Imputation
#            software such as Beagle writes DS alongside GT. When any record lacks DS the whole
#            file falls back to GT, so a file is read one way or the other, never half and half.
#            The dosage columns are then floats rather than whole numbers, which every LD and
#            blocking function here accepts but which downstream code expecting 0 / 1 / 2 may not.
read_vcf_geno = function(path, prefer_ds = FALSE) {
  records = read_vcf_records(path)

  # DS when asked for and available, GT otherwise. field_values is one character vector per
  # sample either way, so only the converter differs.
  field_values = NULL
  if (prefer_ds) {
    field_values = read_vcf_field(records, "DS", required = FALSE)
    converter = ds_to_dosage
  }
  if (is.null(field_values)) {
    field_values = read_vcf_field(records, "GT")
    converter = gt_to_dosage
  }

  # One dosage column per sample, in the original sample order and names.
  # Any converter error is combined with the file it was read from.
  dosage_cols = tryCatch(
    lapply(field_values, converter),
    error = function(e) stop("Cannot read ", path, " as dosages: ", conditionMessage(e), call. = FALSE)
  )

  result = data.frame(records$map, dosage_cols, stringsAsFactors = FALSE, check.names = FALSE)
  row.names(result) = NULL

  return(result)
}

##### Read a (optionally gzipped) VCF back into a HapSelect phased genotype data frame #####
# Splits each phased GT call into its two haplotypes instead of summing them to a dosage, giving
# the geno_phased layout: SNP / Chromosome / Position, then <sample>_1 and <sample>_2 columns of
# 0/1 allele presence. This is the format compute_haplotype_effects() expects.
# path: VCF to read, optionally gzipped
read_vcf_phased = function(path) {
  parsed = read_vcf_gt(path)

  # Two columns per sample, kept adjacent so the order is S1_1, S1_2, S2_1, S2_2, ...
  hap_cols = list()
  for (s in names(parsed$gt)) {
    # Any gt_to_haplotypes() error is combined with the file, the sample and the fix
    haplotypes = tryCatch(
      gt_to_haplotypes(parsed$gt[[s]]),
      error = function(e) stop(
        "Cannot read ", path, " as phased haplotypes, in sample ", s, ": ", conditionMessage(e),
        "\nread_vcf_phased() needs a phased VCF such as Beagle output. Use beagle_phase_geno() ",
        "to phase, or read_vcf_geno() to read this file as dosages.",
        call. = FALSE
      )
    )
    hap_cols[[paste0(s, "_1")]] = haplotypes[, 1]
    hap_cols[[paste0(s, "_2")]] = haplotypes[, 2]
  }

  result = data.frame(parsed$map, hap_cols, stringsAsFactors = FALSE, check.names = FALSE)
  row.names(result) = NULL

  return(result)
}
