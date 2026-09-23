########################################
########## VCF Read and Write ##########
########################################
#
# Reading and writing VCF files, and converting VCF GT calls to and from the
# HapSelect geno (dosage) and geno_phased (haplotype) layouts. Any part of the
# package that needs genotypes in or out of a VCF uses these rather than parsing
# a VCF itself.

##### Check a ploidy argument #####
# Every reader and writer here takes the ploidy as declared rather than inferring it from the
# file, so that a call with the wrong number of alleles is caught instead of read as something else.
# ploidy: number of allele copies per genotype, a single whole number of at least 1
# Returns ploidy as an integer.
check_ploidy = function(ploidy) {
  if (!is.numeric(ploidy) || length(ploidy) != 1 || !is.finite(ploidy) ||
      ploidy < 1 || ploidy != round(ploidy)) {
    stop("ploidy must be a single whole number of at least 1, e.g. 2L for diploid or 4L for tetraploid.")
  }
  as.integer(ploidy)
}

##### Write a genotype data frame out as a minimal VCF #####
# geno  : data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#         cols 4+ = dosage values (0 to ploidy, or NA) per individual
# path  : file path to write the VCF to
# ploidy: allele copies per genotype, 2L for diploid
# Each dosage d is written as an unphased call of ploidy - d REF alleles followed by d ALT alleles,
# so a tetraploid dosage of 1 is 0/0/0/1, and a missing value as one "." per allele (./././.).
# An unphased call records only how many copies carry each allele, not which copy carries which,
# so it holds exactly what a dosage does and nothing is lost by always writing REF first.
# Placeholder REF/ALT alleles (A/G) are used since dosage alone does not carry allele identity.
write_vcf_geno = function(geno, path, ploidy = 2L) {
  ploidy = check_ploidy(ploidy)

  # Use the geno column names as VCF sample IDs, falling back to generated names if missing
  sample_names = colnames(geno)[-(1:3)]
  if (is.null(sample_names) || any(!nzchar(sample_names))) {
    sample_names = paste0("S", seq_len(ncol(geno) - 3))
  }

  # Extract the dosage matrix, stripping out marker/chrom/position
  geno_matrix = as.matrix(geno[, -(1:3), drop = FALSE])

  # A GT call can only express a whole number of ALT copies between none and all of them
  if (any(!(geno_matrix[!is.na(geno_matrix)] %in% 0:ploidy))) {
    stop("Genotype dosages must be whole numbers from 0 to ", ploidy, ", or NA, for ploidy ",
         ploidy, ".")
  }

  # Look up the unphased VCF genotype call for each dosage, built once for every possible dosage
  gt_calls = vapply(0:ploidy, function(d) {
    paste(rep(c("0", "1"), c(ploidy - d, d)), collapse = "/")
  }, character(1))
  names(gt_calls) = 0:ploidy
  gt = matrix(gt_calls[as.character(geno_matrix)], nrow = nrow(geno_matrix))
  gt[is.na(geno_matrix)] = paste(rep(".", ploidy), collapse = "/")

  # Build the VCF body: one row per marker, with placeholder REF/ALT alleles
  body = data.frame(
    geno[[2]], geno[[3]], geno[[1]], "A", "G", ".", "PASS", ".", "GT", gt,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  colnames(body) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_names)

  # Open the output file for writing
  con = file(path, "w")
  on.exit(close(con))

  # Write a minimal VCF header: the file format and the GT field
  writeLines(
    c("##fileformat=VCFv4.2",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"),
    con
  )

  # Write the body rows, tab-delimited, with the #CHROM line as the column header
  utils::write.table(body, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

##### Convert a vector of VCF GT calls to ALT-allele dosages #####
# gt    : GT field values, one per sample, e.g. c("0/0", "0|1", "./.")
# ploidy: allele copies per genotype, 2L for diploid
# Each call is split into its alleles on "/" or "|", so phased, unphased and mixed calls are all
# read, and its dosage is the number of ALT (1) alleles: 0/0/1/1 is 2. A call with any unknown
# allele - ".", "./.", or a partly missing 0/. - becomes NA, since its count is not known.
# Anything else is an error rather than a silent NA: an allele other than 0, 1 or . is a second
# ALT allele (1/2), which one ALT count cannot express, and a call whose allele count differs
# from ploidy means the file is not the ploidy declared. A bare "." is exempt from the count, being
# the VCF spelling for a genotype that was not called at all.
# Only the distinct calls are parsed, then matched back, since a file holds only a handful of them.
gt_to_dosage = function(gt, ploidy = 2L) {
  ploidy = check_ploidy(ploidy)

  calls = unique(gt)
  alleles = strsplit(calls, "[/|]")

  # A NA or empty call splits into nothing usable, so it is caught here too
  recognised = vapply(alleles, function(a) {
    length(a) > 0 && !anyNA(a) && all(a %in% c("0", "1", "."))
  }, logical(1))
  if (any(!recognised)) {
    stop("Unrecognised GT call(s): ", paste(calls[!recognised], collapse = ", "),
         ". Alleles must be 0 (REF), 1 (ALT) or . (missing); multi-allelic calls are not supported.")
  }

  # Check ploidy count
  wrong_count = lengths(alleles) != ploidy & calls != "."
  if (any(wrong_count)) {
    stop("GT call(s) with the wrong number of alleles for ploidy ", ploidy, ": ",
         paste(calls[wrong_count], collapse = ", "), ". Set ploidy to match the VCF.")
  }

  dosage = vapply(alleles, function(a) {
    if (any(a == ".")) NA_real_ else sum(a == "1")
  }, numeric(1))

  dosage[match(gt, calls)]
}

##### Convert a vector of phased VCF GT calls to haplotype allele columns #####
# gt    : phased GT values, one per sample, e.g. c("0|0", "0|1")
# ploidy: allele copies per genotype, 2L for diploid
# Returns an integer matrix of 0/1 allele presence with one column per haplotype, ploidy columns
# in all. The haplotype layout needs every allele placed on a haplotype, so each call must be fully
# phased ("|" between every allele) with exactly ploidy alleles, each 0 or 1. An unphased ("0/1")
# or missing (".|.") call is an error rather than an NA, since a haplotype column cannot hold one.
# A haploid call ("0" or "1") has a single allele, so there is nothing to phase.
# Only the distinct calls are parsed, then matched back, as in gt_to_dosage().
gt_to_haplotypes = function(gt, ploidy = 2L) {
  ploidy = check_ploidy(ploidy)

  calls = unique(gt)
  alleles = strsplit(calls, "|", fixed = TRUE)

  # "/" is not split on, so an unphased call leaves an allele such as "0/1" that fails here
  recognised = vapply(alleles, function(a) {
    length(a) == ploidy && !anyNA(a) && all(a %in% c("0", "1"))
  }, logical(1))
  if (any(!recognised)) {
    example = paste(c(rep("0", ploidy - 1), "1"), collapse = "|")
    stop("Unrecognised GT call(s): ", paste(calls[!recognised], collapse = ", "),
         ". Expected phased calls of ", ploidy, " alleles, each 0 or 1, separated by \"|\", e.g. ",
         example, ".")
  }

  haplotypes = matrix(as.integer(unlist(alleles)), ncol = ploidy, byrow = TRUE)
  haplotypes[match(gt, calls), , drop = FALSE]
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
      sample_cols = sample_cols
    ))
  }

  # Split every line into its tab-delimited fields
  fields = strsplit(data_lines, "\t", fixed = TRUE)

  # rbind() would pad a short record by reusing its own values, so check every record first
  if (any(lengths(fields) != length(col_names))) {
    stop("VCF records do not all have the ", length(col_names), " tab-separated columns the ",
         "#CHROM header names: ", path, call. = FALSE)
  }
  tab = do.call(rbind, fields)
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
    path = path
  )
}

##### Pull one FORMAT field out of every sample column of a parsed VCF #####
# Each record's FORMAT (e.g. "GT" or "GT:DS:GP") names the colon-separated values in its sample
# cells, in order. FORMAT can differ from record to record, but records that share a FORMAT string
# hold the field at the same position, so each distinct FORMAT is handled once, over all of its
# records and samples together, rather than cell by cell. A file usually has a single FORMAT.
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

  format_strings = records$tab[, "FORMAT"]
  # Unique format strings
  formats = unique(format_strings)
  # Each distinct FORMAT split into its field names, e.g. "GT:DS" -> c("GT", "DS")
  format_fields = strsplit(formats, ":", fixed = TRUE)
  # Position of the requested field within each distinct FORMAT, NA where it is not listed
  field_index = vapply(format_fields, function(f) match(field, f), integer(1))

  # A FORMAT that does not list the field: an error, or NULL when the field is optional
  if (any(is.na(field_index))) {
    if (!required) return(NULL)
    stop("VCF FORMAT field does not include ", field, " for one or more records: ", records$path)
  }

  # How many values each distinct FORMAT holds, so a single-field FORMAT can skip the split
  n_fields = lengths(format_fields)

  # The records under each distinct FORMAT, in file order
  format_rows = split(seq_along(format_strings),
                      factor(match(format_strings, formats), levels = seq_along(formats)))

  # For each sample column, pull out the field, giving one character vector per sample
  values = lapply(records$sample_cols, function(s) {
    column = records$tab[, s]

    # Fast path: one FORMAT for the whole file, so the column is taken in one call
    if (length(formats) == 1) return(extract_format_value(column, field_index, n_fields))

    # Otherwise, one call per FORMAT group, writing each group's values back into its records' positions
    out = rep(NA_character_, length(column))
    for (g in seq_along(formats)) {
      rows = format_rows[[g]]
      out[rows] = extract_format_value(column[rows], field_index[g], n_fields[g])
    }
    out
  })
  stats::setNames(values, records$sample_cols)
}

##### Pull the k-th colon-separated value out of a set of VCF sample cells #####
# One vectorised call over every cell, rather than a split and an R function call per cell.
# cells    : character vector of sample cells that share one FORMAT
# k        : position of the wanted value within that FORMAT
# n_fields : number of values that FORMAT holds
# Returns the values in the same order as cells, NA where a cell stops before position k.
extract_format_value = function(cells, k, n_fields) {
  # A FORMAT of a single field means the cell is the value, so there is nothing to split
  if (n_fields == 1) return(cells)

  # Skip k - 1 values and their colons, keep the next value, drop the rest
  pattern = paste0("^(?:[^:]*:){", k - 1, "}([^:]*).*$")
  out = sub(pattern, "\\1", cells, perl = TRUE)

  # sub() returns a cell it cannot match unchanged. Past the first value, a cell that matched has
  # lost at least one colon, so an unchanged cell is one that stops before position k
  if (k > 1) out[out == cells] = NA_character_

  out
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
# DS is the imputed ALT allele dosage: a float on the same scale as a GT-derived dosage, but
# carrying the imputation's uncertainty rather than rounding it away, so 0.83 stays 0.83 instead of
# becoming the hard call 1. A missing value (".") becomes NA, as it does for GT.
# ds    : DS field values, one per sample, e.g. c("0.00", "0.83", ".")
# ploidy: allele copies per genotype, 2L for diploid, which bounds the dosage from above
ds_to_dosage = function(ds, ploidy = 2L) {
  ploidy = check_ploidy(ploidy)
  dosage = suppressWarnings(as.numeric(ds))

  # A missing value is expected and becomes NA; anything else that will not parse is not
  unknown = is.na(dosage) & !(ds %in% c(".", "") | is.na(ds))
  if (any(unknown)) {
    stop("Unrecognised DS value(s): ", paste(unique(ds[unknown]), collapse = ", "),
         ". DS must hold numeric ALT allele dosages.")
  }

  # A dosage counts ALT copies, so it lies between none and all of them. A DS outside that range
  # usually means the field was misread, or that the file is not the ploidy declared
  out_of_range = !is.na(dosage) & (dosage < 0 | dosage > ploidy)
  if (any(out_of_range)) {
    stop("DS value(s) outside the 0 to ", ploidy, " dosage range for ploidy ", ploidy, ": ",
         paste(unique(dosage[out_of_range]), collapse = ", "), ". Set ploidy to match the VCF.")
  }

  dosage
}

##### Read a (optionally gzipped) VCF back into a HapSelect genotype data frame #####
# Converts the GT (genotype) field of every sample column to a dosage, counting ALT alleles, with
# missing calls as NA; phased ("|") and unphased ("/") genotypes are both accepted. See
# gt_to_dosage() for which calls are read and which are errors.
# Returns a geno object, so columns are SNP / Chromosome / Position, matching order_map().
# LD and haploblock tables use Chrom instead - do not align this to those.
# path     : VCF to read, optionally gzipped
# prefer_ds: TRUE to read dosages from the DS field instead of GT when every record carries one,
#            keeping fractional imputed dosages rather than the rounded hard calls. Imputation
#            software commonly writes DS alongside GT. When any record lacks DS the whole file
#            falls back to GT, so a file is read one way or the other, never half and half.
#            The dosage columns are then floats rather than whole numbers, which every LD and
#            blocking function here accepts but which downstream code expecting whole numbers
#            may not. Only the DS range is checked against ploidy; GT is not read as well.
# ploidy   : allele copies per genotype, 2L for diploid. Every GT call must carry this many
#            alleles, and dosages run from 0 to ploidy
read_vcf_geno = function(path, prefer_ds = FALSE, ploidy = 2L) {
  ploidy = check_ploidy(ploidy)
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
    lapply(field_values, converter, ploidy = ploidy),
    error = function(e) stop("Cannot read ", path, " as dosages: ", conditionMessage(e), call. = FALSE)
  )

  result = data.frame(records$map, dosage_cols, stringsAsFactors = FALSE, check.names = FALSE)
  row.names(result) = NULL

  return(result)
}

##### Read a (optionally gzipped) VCF back into a HapSelect phased genotype data frame #####
# Splits each phased GT call into its haplotypes instead of summing them to a dosage, giving the
# geno_phased layout: SNP / Chromosome / Position, then <sample>_1 to <sample>_<ploidy> columns
# of 0/1 allele presence. 
# path  : VCF to read, optionally gzipped
# ploidy: allele copies per genotype, 2L for diploid. Every GT call must carry this many alleles
read_vcf_phased = function(path, ploidy = 2L) {
  ploidy = check_ploidy(ploidy)
  parsed = read_vcf_gt(path)

  # ploidy columns per sample, kept adjacent so the order is S1_1, S1_2, S2_1, S2_2, ...
  hap_cols = list()
  for (s in names(parsed$gt)) {
    # Any gt_to_haplotypes() error is combined with the file, the sample and the fix
    haplotypes = tryCatch(
      gt_to_haplotypes(parsed$gt[[s]], ploidy = ploidy),
      error = function(e) stop(
        "Cannot read ", path, " as phased haplotypes, in sample ", s, ": ", conditionMessage(e),
        "\nread_vcf_phased() needs a VCF whose calls are all phased. Phase the genotypes first, ",
        "or use read_vcf_geno() to read this file as dosages.",
        call. = FALSE
      )
    )
    for (k in seq_len(ploidy)) {
      hap_cols[[paste0(s, "_", k)]] = haplotypes[, k]
    }
  }

  result = data.frame(parsed$map, hap_cols, stringsAsFactors = FALSE, check.names = FALSE)
  row.names(result) = NULL

  return(result)
}
