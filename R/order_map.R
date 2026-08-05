#####################################
##### Map file Check and Order ######
#####################################

####order_chromo() - function to order the map file for an individual chromosome (called inside of overall map order function)####
#Given the map of a chromosome (take whole map file and split by chromo), order it based on marker position
order_chromo = function(chromo){
  #use the position column to order
  chromo = chromo[order(chromo[,3]),]
  return(chromo)
}

####parse_chromosome() - convert a chromosome column to numbers####
# Chromosomes that already carry a number are kept. Labels with no numeric are assigned one.
#
# Examples:
#   "1", "10"                  -> 1, 10      (kept as-is, never renumbered)
#   "chr2", "ch9", "Chr03"     -> 2, 9, 3    (leading non-numeric prefix stripped)
#   "X", "MT", "chrX", "1A"    -> Assigned numbers above the highest real
#                                 chromosome, in sorted label order, so they are
#                                 deterministic and cannot collide with a real one
#   NA, ""                     -> NA
#
# x : chromosome column (numeric, character or factor)
parse_chromosome = function(x){

  # If the column is already numeric, return it as-is.
  if(is.numeric(x)) return(x)

  raw = trimws(as.character(x))
  # Set empty strings to NA
  raw[which(raw == "")] = NA

  # A chromosome column holds only a handful of distinct labels, so each one is
  # parsed once and matched back rather than being parsed row by row.
  labels = unique(raw)

  # Drop any leading non-digit characters, then read the number that remains
  digits = as.character(0:9)
  numbers = vapply(labels, function(label){
    # If the label is NA, return NA.
    if(is.na(label)) return(NA_real_)

    # Split into characters and find the first digit.
    chars = strsplit(label, "", fixed = TRUE)[[1]]
    first_match = match(TRUE, chars %in% digits)

    # If there is no digit, return NA.
    if(is.na(first_match)) return(NA_real_)

    # Finally, extract the substring starting at the first digit and convert to numeric.
    suppressWarnings(as.numeric(substring(label, first_match)))
  }, numeric(1), USE.NAMES = FALSE)

  # Parse the original vector by matching the raw labels to the unique labels and indexing into the parsed numbers.
  parsed = numbers[match(raw, labels)]

  # The labels carrying no number at all, e.g. X, Y, MT, chrX, scaffold_a
  unparsed = is.na(parsed) & !is.na(raw)

  message = paste("Chromosomes were not numeric - parsed to numeric (e.g. 'chr10' -> 10).")

  if(any(unparsed)){
    labels   = sort(unique(raw[unparsed]))
    highest  = if(all(is.na(parsed))) 0 else max(parsed, na.rm = TRUE)
    assigned = setNames(highest + seq_along(labels), labels)

    parsed[unparsed] = assigned[raw[unparsed]]

    message = paste0(message, " Labels with no number of their own were assigned one: ",
                     paste(names(assigned), assigned, sep = " -> ", collapse = ", "), ".")
  }

  # If any labels were assigned a number, issue a warning with the message.
  warning(message)
  return(parsed)
}


#####check file structure######
check_file = function(map){
  # Check file structure - make sure it's a data frame with at least 3 columns (SNP, chrom, pos)
  if(!is.data.frame(map) || ncol(map) < 3){
    stop("map must be a data frame with at least 3 columns: SNP ID (column 1), chromosome (column 2, numeric), and position (column 3, numeric).")
  }

  #make sure SNP ID are characters - if not make them characters and give a warning
  if(!is.character(map[,1])){
    map[,1] = as.character(map[,1])
    warning("SNP ID were not characters - coercing to characters. For proper function, ensure they are characters in other files and they match this output.")
  }

  #check chromosomes are numeric
  if(!is.numeric(map[,2])){
    map[,2] = parse_chromosome(map[,2])
  }

  #check positions are numeric
  if(!is.numeric(map[,3])){
    map[,3] = as.numeric(as.character(map[,3]))
    warning("Positions were not numeric - attempting to coerce to numeric. Check the output is correct. For proper function, ensure they are numeric in other files and they match this output.")
  }

  return(map)
}

####order_map() - order the entire map file, which calls the order_chromo() function####
#provide the whole map file with columns "SNP" (name of the snp), "chrom", and "pos"
order_map = function(map){
  #check the files
  map = check_file(map)
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


