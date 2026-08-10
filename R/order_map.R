#####################################
##### Map file Check and Order ######
#####################################

####natural_key() - sort key that compares digits as numbers####
# String sorting works one character at a time, so "chr10" sorts before "chr2".
# Padding every number to the same width lines the digits up and fixes that:
#
#   "chr2" -> "chr000000000002",  "chr10" -> "chr000000000010"
#
# Rather than measure each number, we add more zeros than any number could need and
# then trim back to `width`. Labels with no digits ("chrX") pass through untouched
# and sort after the numbers.
#
# x     : character vector
# width : digits every number is padded to. Longer numbers are left unpadded.
natural_key = function(x, width = 12){

  # Pad every run of digits, leaving each number at least `width` long. "[0-9]+" is
  # greedy, so a run is padded once rather than digit by digit.
  x = gsub("([0-9]+)", paste0(strrep("0", width - 1), "\\1"), x)

  # Trim to the last `width` digits. "0*" sits outside the capture group, so the
  # zeros it greedily eats - exactly the surplus added above - are the ones dropped.
  return(gsub(paste0("0*([0-9]{", width, "})"), "\\1", x))
}

####number_chromosomes() - map a chromosome column onto integers####
# Every distinct label is treated as a distinct chromosome.
#
# Labels that are already plain numbers keep their value, which leaves a numeric
# map untouched. Every other label is assigned a number above the highest of
# them, in natural order so that related names stay in the order their digits
# imply. Anything ambiguous emits a warning.
#
# Examples:
#   "1", "2", "10"                    -> 1, 2, 10   (kept, never renumbered)
#   "1", "2", "10", "X"               -> 1, 2, 10, 11
#   "chr1", "chr2", "chr10"           -> 1, 2, 3
#   "11", "12", "scaffold_12"         -> 11, 12, 13
#   NA, ""                            -> NA
#
# x       : chromosome column (numeric, character or factor)
# verbose : report the label to integer mapping
number_chromosomes = function(x, verbose = TRUE){

  # The column is already a set of numerics, so there is nothing to do.
  if(is.numeric(x)) return(x)

  # Convert to character and treat empty strings as missing. One entry per row.
  label_per_row = trimws(as.character(x))
  label_per_row[!is.na(label_per_row) & label_per_row == ""] = NA

  # Create a vector of distinct labels, duplicates removed and missing values excluded.
  distinct_labels = unique(label_per_row[!is.na(label_per_row)])

  # There is nothing to do if there are no valid labels, so return the column as-is.
  if(!length(distinct_labels)) return(as.numeric(label_per_row))

  # Sort the labels in natural order so that "chr2" comes before "chr10".
  sorted_labels = distinct_labels[order(natural_key(distinct_labels), method = "radix")]

  # Whether a label is a plain number is determined by a regular expression that matches only digits.
  is_number = grepl("^[0-9]+$", sorted_labels)

  # Create an empty numeric vector of the same length as the labels, initialized with NA.
  number_per_label = rep(NA_real_, length(sorted_labels))

  # Fill in the numbers for labels that are plain numbers
  number_per_label[is_number] = as.numeric(sorted_labels[is_number])
  
  # Assign numbers to labels that are not plain numbers, starting from one above the maximum of the existing numbers.
  max_number = max(c(0, number_per_label), na.rm = TRUE)
  number_per_label[!is_number] = max_number + seq_len(sum(!is_number))

  # Warning: Similar spelling of a chromosome label
  # Finds labels that are very close in spelling (ignoring case and punctuation) and warns the user that they were numbered separately.
  similar = split(sorted_labels, tolower(gsub("[^[:alnum:]]", "", sorted_labels)))
  similar = similar[lengths(similar) > 1]
  if(length(similar)){
    warning("These chromosome labels differ only in case or punctuation and were ",
            "numbered separately. You should merge these yourself if they are the same ",
            "chromosome: ",
            paste(vapply(similar, paste, character(1), collapse = " / "), collapse = ", "), ".")
  }

  # Warning: Labels that look like missing values
  bad_label_markers = c("na", "n/a", "nan", "null", ".", "-", "?", "-9", "unknown")
  bad_labels = sorted_labels[tolower(sorted_labels) %in% bad_label_markers]
  if(length(bad_labels)){
    warning("These chromosome labels look like missing values but were numbered ",
            "as real chromosomes: \n", paste(bad_labels, collapse = ", "),
            ". Set them to NA if they are missing.")
  }

  # If verbose is TRUE, display a message showing the mapping of chromosome labels to integers.
  if(verbose){
    shown = order(number_per_label)[seq_len(length(sorted_labels))]
    message("Chromosome labels were mapped to integers:\n",
            paste0("  ", format(sorted_labels[shown]), " -> ", number_per_label[shown], collapse = "\n"))
  }

  # Return the numeric vector corresponding to the original input, using the mapping from labels to numbers.
  return(unname(setNames(number_per_label, sorted_labels)[label_per_row]))
}

#####check file structure######
# map     : map file, with SNP ID, chromosome and position in columns 1 to 3
# verbose : passed to number_chromosomes() to report the chromosome numbering
check_file = function(map, verbose = FALSE){
  # Check file structure - make sure it's a data frame with at least 3 columns (SNP, chrom, pos)
  if(!is.data.frame(map) || ncol(map) < 3){
    stop("map must be a data frame with at least 3 columns: SNP ID (column 1), chromosome (column 2, numeric), and position (column 3, numeric).")
  }

  # Make sure SNP ID are characters - if not make them characters and give a warning
  if(!is.character(map[,1])){
    map[,1] = as.character(map[,1])
    warning("SNP ID were not characters - coercing to characters. For proper function, ensure they are characters in other files and they match this output.")
  }

  # Make sure chromosomes are numeric
  if(!is.numeric(map[,2])){
    # If chromosomes are not numeric, call number_chromosomes() to convert them to integers
    map[,2] = number_chromosomes(map[,2], verbose = verbose)
  }

  # Check that positions are numeric. If not make them numeric and give a warning
  if(!is.numeric(map[,3])){
    map[,3] = as.numeric(as.character(map[,3]))
    warning("Positions were not numeric - attempting to coerce to numeric. Check the output is correct. For proper function, ensure they are numeric in other files and they match this output.")
  }

  return(map)
}

####order_map() - order the entire map file by chromosome, then by position####
#provide the whole map file with columns "SNP" (name of the snp), "chrom", and "pos"
#
# map     : map file, with SNP ID, chromosome and position in columns 1 to 3
# verbose : report the chromosome numbering when the chromosome column is not
#           already numeric
order_map = function(map, verbose = FALSE){
  #check the files
  map = check_file(map, verbose = verbose)
  colnames(map)[1:3] = c("SNP", "Chromosome", "Position")

  #the chromosome column is numeric by this point, so a single sort by chromosome
  #then position orders the whole map. Markers with a missing chromosome sort last.
  map = map[order(map[,2], map[,3]), ]
  rownames(map) = NULL

  return(map)
}


