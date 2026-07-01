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

#####check file structure######
check_file = function(map){
  # Check file structure - make sure it's a data frame with at least 3 columns (SNP, chrom, pos)
  if(!is.data.frame(map) || ncol(map) < 3){
    stop("map must be a data frame with at least 3 columns: SNP ID (column 1), chromosome (column 2, numeric), and position (column 3, numeric).")
  }

  #make sure SNP ID are characters - if not make them characters and give a warning
  if(is.numeric(map[,1])){
    map[,1] = as.character(map[,1])
    warning("SNP ID are numeric - coercing to characters. For proper function, ensure they are characters in other files and they match this output.")
  }

  #check chromosomes are numeric
  if(!is.numeric(map[,2])){
    map[,2] = as.numeric(as.factor(map[,2]))
    warning("Chromosomes were not numeric - coercing to numeric. Check the output is correct. For proper function, ensure they are numeric in other files and they match this output.")
  }

  #check that positions are numeric
  if(!is.numeric(map[,3])){
    map[,3] = as.numeric(map[,3])
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


