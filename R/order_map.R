#####################################
##### Map file Check and Order ######
#####################################

####order_chromo() - function to order the map file for an individual chromosome (called inside of overall map order function)####
# Given the map of a chromosome (take whole map file and split by chromo), order it based on marker position
order_chromo = function(chromo){
  #use the position column to order
  chromo = chromo[order(chromo[,3]),]
  return(chromo)
}

####natural_key() - sort key that compares digits as numbers####
# Left-pads every run of digits so that ordinary string sorting puts "chr2"
# before "chr10" rather than after it.
#
# x : character vector
natural_key = function(x){
  runs = gregexpr("[0-9]+", x)
  regmatches(x, runs) = lapply(regmatches(x, runs), function(digits)
    paste0(strrep("0", pmax(0, 12 - nchar(digits))), digits))
  return(x)
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

  # Convert to character and treat empty strings as missing.
  raw = trimws(as.character(x))
  raw[!is.na(raw) & raw == ""] = NA

  # Sort the distinct labels in natural order and assign each a number. 
  labels = unique(raw[!is.na(raw)])

  # There is nothing to do if there are no valid labels, so return the raw column as-is.
  if(!length(labels)) return(as.numeric(raw))

  # Sort the labels in natural order so that "chr2" comes before "chr10".
  labels = labels[order(natural_key(labels), method = "radix")]
  
  # Whether a label is a plain number is determined by a regular expression that matches only digits.
  is_number = grepl("^[0-9]+$", labels)

  # Create an empty numeric vector of the same length as the labels, initialized with NA.
  number = rep(NA_real_, length(labels))
  # Fill in the numbers for labels that are plain numbers
  number[is_number] = as.numeric(labels[is_number])
  # Assign numbers to labels that are not plain numbers, starting from one above the maximum of the existing numbers.
  max_number = max(c(0, number), na.rm = TRUE)
  number[!is_number] = max_number + seq_len(sum(!is_number))

  # Warning: Similar spelling of a chromosome label
  # Finds labels that are very close in spelling (ignoring case and punctuation) and warns the user that they were numbered separately.
  similar = split(labels, tolower(gsub("[^[:alnum:]]", "", labels)))
  similar = similar[lengths(similar) > 1]
  if(length(similar)){
    warning("These chromosome labels differ only in case or punctuation and were ",
            "numbered separately - merge them yourself if they are the same ",
            "chromosome: ",
            paste(vapply(similar, paste, character(1), collapse = " / "), collapse = ", "), ".")
  }

  # Warning: Labels that look like missing values
  bad_label_markers = c("na", "n/a", "nan", "null", ".", "-", "?", "-9", "unknown")
  bad_labels = labels[tolower(labels) %in% bad_label_markers]
  if(length(bad_labels)){
    warning("These chromosome labels look like missing values but were numbered ",
            "as real chromosomes: ", paste(bad_labels, collapse = ", "),
            ". Set them to NA if they are missing.")
  }

  # If verbose is TRUE, display a message showing the mapping of chromosome labels to integers.
  if(verbose){
    shown = order(number)[seq_len(length(labels))]
    message("Chromosome labels were mapped to integers",
            " (apply the same mapping to your genotype and marker effect files):\n",
            paste0("  ", format(labels[shown]), " -> ", number[shown], collapse = "\n"))
  }

  # Return the numeric vector corresponding to the original input, using the mapping from labels to numbers.
  return(unname(setNames(number, labels)[raw]))
}

#####check file structure######
# map     : map file, with SNP ID, chromosome and position in columns 1 to 3
# verbose : passed to number_chromosomes() to report the chromosome numbering
check_file = function(map, verbose = TRUE){
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

####order_map() - order the entire map file, which calls the order_chromo() function####
#provide the whole map file with columns "SNP" (name of the snp), "chrom", and "pos"
#
# map     : map file, with SNP ID, chromosome and position in columns 1 to 3
# verbose : report the chromosome numbering when the chromosome column is not
#           already numeric
order_map = function(map, verbose = FALSE){
  #check the files
  map = check_file(map, verbose = verbose)
  colnames(map)[1:3] = c("SNP", "Chromosome", "Position")

  #create a progress bar - might not be needed as it's so fast
  handlers("txtprogressbar")

  #split the map file up by chromosome and order within chromosome
  map_split = split(map, map[,2])

  #call the progress bar and while it's active do the ordering
  with_progress({
    #define how many times the progress bar should update - can also use `along = list_name` which will
    #automatically define the length based on a list
    p = progressor(steps = length(map_split))

    #iterate over the chromosomes sequentially and combine individual data frames into rows from the list (dfr part of map)
    ordered_map = map_dfr(map_split, function(chromo){

      #call the map ordering function on the map of the chromosome
      chromo = order_chromo(chromo)

      #iterate the progress bar to indicate a step has been completed
      p()

      #return the ordered chromosome map, which will be row binded
      return(chromo)
    })
  })
}


