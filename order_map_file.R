#####load dependencies####
library(purrr)
library(progressr)

load(file = "Example_Files/gapit_map.R")


####order_chromo() - function to order the map file for an individual chromosome (called inside of overall map order function)####
#Given the map of a chromosome (take whole map file and split by chromo), order it based on marker position
order_chromo = function(chromo){
  chromo = chromo[order(chromo$pos),]
  return(chromo)
}


####order_map() - order the entire map file, which calls the order_chromo() function####
#provide the whole map file with columns "SNP" (name of the snp), "chrom", and "pos"
order_map = function(map){
  #create a progress bar - might not be needed as it's so fast
  handlers("txtprogressbar")
  
  #split the map file up by chromosome and order within chromosome
  map_split = split(map, map$chromo)
  
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


