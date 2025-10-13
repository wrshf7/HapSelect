####load dependencies####

library(purrr)
#library(furrr)
#library(future)
#library(parallel)
library(dplyr)
library(progressr)


#####function to extend the block left####

extend_left = function(first.mkr.blk, snps.chr, assigned, ld.chr, block, method, threshold, tolerance, tol_reset){
  
  #start counting LD pairs that didn't meet the threshold and track the marker for the next iteration to see if it should be included
  tolerance_counter = 0
  failed_markers = c()
  
  #repeat left extension until break point is reached (tolerance and threshold)
  repeat {
    
    #extract index of the first marker in the block from the SNP pair
    first.idx = match(first.mkr.blk, snps.chr)
    
    #don't extend left if it was the first marker
    if (first.idx == 1) break
    
    #get index of marker to the left of the first SNP in the block - modify for skipped markers from tolerance fail
    left.mkr = snps.chr[first.idx - (length(failed_markers) + 1) ]
    
    #added check if tolerance checks for failed_markers was out of bounds
    if(length(left.mkr) == 0) break
    if(is.na(left.mkr)) break
    
    #if marker to the left is already in a block, break
    if (assigned[left.mkr]) break
    
    #if utilizing the flanking method, pull the LD of the SNP to the left idetnified from left.mkr
    if (method=="flanking") {
      #singular value?
      ld.vals = ld.chr[ld.chr$Name1 == left.mkr & ld.chr$Name2 == first.mkr.blk, "LD"]
    }
    
    #if utilizing the average method, then compute the average LD of the new SNP considered with all markers in the block currently
    if (method=="average") {
      #could be multiple values if the block is greater than
      ld.vals = ld.chr[ld.chr$Name1 == left.mkr & ld.chr$Name2 %in% block, "LD"]
    }
    
    #check that the LD value (flanking) or average LD value (average) meets the threshold - mean of a single marker is itself, so it works for both
    if (mean(ld.vals) <= threshold){
      
      #check to see if the tolerance threshold has been met, if not add one to the counter
      if(tolerance_counter >= tolerance){
        break
      }
      else{
        tolerance_counter = tolerance_counter + 1
        failed_markers = c(left.mkr, failed_markers)
        next
      }
    } 
    
    #reset tolerance counter if option is provided
    if(tol_reset){
      tolerance_counter = 0
    }
    
    #if the threshold is met, add it to the block - check to see if a marker had been skipped as part of the tolerance and threshold check
    if(length(failed_markers) > 0){
      block = c(left.mkr, failed_markers, block)
      
      #reset to track for further extension
      failed_markers = c()
    } else{
      block = c(left.mkr, block)
    }
    
    
    #update the first marker and repeat the loop
    first.mkr.blk = left.mkr
  }
  
  
  #return relevant objects
  return(block)
}

#####function to extend the block right#####

extend_right = function(last.mkr.blk, snps.chr, assigned, ld.chr, block, method, threshold, tolerance, tol_reset){
  
  tolerance_counter = 0
  failed_markers = c()
  
  #keep extending until the break signal is initiated
  repeat{
    #starting from the second marker of the original pair (last marker)
    last.idx = match(last.mkr.blk, snps.chr)
    
    #if it's the last marker in the chromo, don't extend
    if (last.idx == length(snps.chr)) break
    
    #if not, find the next marker - modify by failed markers when considering skipped markers from tolerance
    right.mkr = snps.chr[last.idx + (length(failed_markers) + 1)]

    #added check if tolerance checks for failed_markers was out of bounds
    if(length(right.mkr) == 0) break
    if(is.na(right.mkr)) break
    
    #check if the next marker is already in a block
    if (assigned[right.mkr]) break
    
    #same as left extension - pull the LD of the marker to the right or the average of the new marker with the rest of the markers in the block originating from the original pair
    if (method=="flanking") {
      ld.vals = ld.chr[ld.chr$Name1 == last.mkr.blk & ld.chr$Name2 == right.mkr, "LD"]
    }
    
    if (method=="average") {
      ld.vals = ld.chr[ld.chr$Name1 %in% block & ld.chr$Name2 == right.mkr, "LD"]
    }
    
    #check to see if it meets threshold criteria and check tolerance counter and threshold like left block extension
    if (mean(ld.vals) <= threshold){
      
      #check to see if the tolerance threshold has been met, if not add one to the counter
      if(tolerance_counter >= tolerance){
        break
      }
      else{
        tolerance_counter = tolerance_counter + 1
        failed_markers = c(failed_markers, right.mkr)
        next
      }
    } 
    
    #reset tolerance counter if option is provided
    if(tol_reset){
      tolerance_counter = 0
    }

    #if the threshold is met, add it to the block - check to see if a marker had been skipped as part of the tolerance and threshold check
    if(length(failed_markers) > 0){
      block = c(block, failed_markers, right.mkr)
      
      #reset to track for further extension
      failed_markers = c()
      
    } else{
      block = c(block, right.mkr)
    }
    
    #update the last marker
    last.mkr.blk = right.mkr
    
  }
  return(block)
}


#####function to find a SNP pair to start block extension#####

make_blocks = function(ld.chr, ld.adj.chr, snps.chr, snps.pos.chr, first.mkr.chr, last.mkr.chr, 
                       assigned, method, threshold, tolerance, tol_reset, start){
  
  chromo_blocks = list()
  
  if(start == "LD"){
    #keep making blocks until there is no pair with high enough LD
    repeat {
      
      #identify max LD pair and pull its id
      max.ld.idx = which.max(ld.adj.chr$LD)
      
      #extract that LD value and check it
      max.ld = ld.adj.chr$LD[max.ld.idx]
      if (max.ld <= threshold) break
      
      #Pull the two marker names in the LD pair
      first.mkr.blk = as.character(ld.adj.chr$Name1[max.ld.idx])
      last.mkr.blk  = as.character(ld.adj.chr$Name2[max.ld.idx])
      
      #if assigned is true for either marker in the LD pair (assigned to a block already), remove it from consideration - updated at the end of the block formation
      if (assigned[first.mkr.blk] || assigned[last.mkr.blk]) {
        ld.adj.chr = ld.adj.chr[-max.ld.idx, ]
        next
      }
      
      #extend the block
      block = c(first.mkr.blk, last.mkr.blk)
      
      #extend the block
      block = extend_left(first.mkr.blk = first.mkr.blk, snps.chr = snps.chr, assigned = assigned, 
                          ld.chr = ld.chr, block = block, method = method, threshold = threshold, tolerance = tolerance, tol_reset = tol_reset)
      
      block = extend_right(last.mkr.blk = last.mkr.blk, snps.chr = snps.chr, assigned = assigned, 
                           ld.chr = ld.chr, block = block, method = method, threshold = threshold, tolerance = tolerance, tol_reset = tol_reset)
      
      #update the assigned variable to indicate that SNP have been assigned to a block
      assigned[block] = TRUE
      
      
      chromo_blocks[[length(chromo_blocks) + 1]] = block
      
      
      #ask victor about this code - we need to figure out whether we can use this or assigned
      #probable answer - the last snp doesn't have an "adjacent" snp to the right, so ld.adj.chr has 1 less
      #snp than snp.chr. The used variable is therefore 1 shorter
      used = with(ld.adj.chr, Name1 %in% block | Name2 %in% block)
      ld.adj.chr = ld.adj.chr[!used, ]
      
      
      if (nrow(ld.adj.chr) == 0) break
      
      
    }
  } else if(start == "beginning"){
    
    position = 1
    total_SNP = length(snps.chr)

    repeat {
      
      #Pull the two marker names in the LD pair
      last.mkr.blk = as.character(snps.chr[position])
      
      block = last.mkr.blk
      
      #extend the block
      block = extend_right(last.mkr.blk = last.mkr.blk, snps.chr = snps.chr, assigned = assigned, 
                           ld.chr = ld.chr, block = block, method = method, threshold = threshold, tolerance = tolerance, tol_reset = tol_reset)
      
      #update the assigned variable to indicate that SNP have been assigned to a block
      assigned[block] = TRUE
      
      
      chromo_blocks[[length(chromo_blocks) + 1]] = block
      
      #update position
      position = match(block[length(block)], snps.chr) + 1
      
      if (position > length(snps.chr)) break
    }
  }
  chromo_blocks = list(chromo_blocks, assigned)
  return(chromo_blocks)
}


#####master blocking function at the chromosome level#####
chromo_blocking = function(chr, ld, map, method, tolerance, tol_reset, threshold, start){
  
  #pull the LD information for the chromosome
  ld.chr = ld[ld$Chrom == chr, ]
  
  #pulling the immediately adjacent SNP to find all SNP pairs
  ld.adj.chr = ld.chr[ld.chr$Locus2 == ld.chr$Locus1 + 1, ]
  
  #pull all snp name on the chromosome
  snps.chr = unique(c(ld.chr$Name1, ld.chr$Name2))
  
  #match snp name to the position in the map file
  snps.pos.chr = map$Position[match(snps.chr, map$SNP)]
  
  #order the snp based on position - should already be ordered
  snps.chr = snps.chr[order(snps.pos.chr)]
  
  #find the first and last SNP on the chromo - utilized to check whether blocks have reached the end of the chromosome 
  first.mkr.chr = snps.chr[1]
  last.mkr.chr = snps.chr[length(snps.chr)]
  
  #tracker utilized later to determine whether a SNP has been added to a block
  assigned = rep(FALSE, length(snps.chr))
  names(assigned) = snps.chr
  
  #return multiple objects - destructure after
  chromo_blocks = make_blocks(ld.chr = ld.chr, ld.adj.chr = ld.adj.chr, snps.chr = snps.chr, snps.pos.chr = snps.pos.chr,
             first.mkr.chr = first.mkr.chr, last.mkr.chr = last.mkr.chr, assigned = assigned, method = method, 
             threshold = threshold, tolerance = tolerance, tol_reset = tol_reset, start = start)

  assigned = chromo_blocks[[2]]
  chromo_blocks = chromo_blocks[[1]]
  
  unassigned = names(assigned)[!assigned]
  
  #add these as blocks
  for (snp in unassigned) {
    chromo_blocks[[length(chromo_blocks) + 1]] = snp
  }
  
  return(chromo_blocks)
}
  

#####overall function to track progress, take input, and call blocking functions#####
def_blocks = function(ld, map, method = "flanking", tolerance = 1, tol_reset = TRUE, threshold = 0.7, start = c("LD", "beginning")){
  
  start = match.arg(start)
  
  #pull the list of chromosomes to parallelize
  chromosomes = sort(unique(ld$Chrom))
  
  #setup parallelization using futures and parallel package and utilize all but 1 core
  #future::plan(multisession, workers = parallel::detectCores() - 1)
  
  #track progress across chromosomes and set up progress bar
  handlers("txtprogressbar")
  with_progress({
    
    #while tracking, set up the progressor for when the progress bar advances
    p = progressor(steps = length(chromosomes))
    
    #form blocks across chromosomes - parallelize each chromosome
    blocks_list = map(chromosomes, function(chr){
      
      #call chromosome blocking function
      chromo_blocks = chromo_blocking(chr = chr, ld = ld, map = map, method = method, tolerance = tolerance, 
                                      threshold = threshold, start = start, tol_reset = tol_reset)
      
      #after blocking on the chromosome is finished iterate the progress bar
      p()
      
      return(chromo_blocks)
    })
  })
  #plan(sequential)
  names(blocks_list) = as.character(chromosomes)
  
  return(blocks_list)
}






#####Turn blocks into a dataframe####
#function to pull SNP locations and chromosome of the block and then order the blocks


#function to condense chromosome blocks into a df
chromo_blocks_to_df = function(chromo, map){
  
  #extract first and last snp from a block, concatenate markers from a vector to a string
  #and return as a df
  chromo_df = map_dfr(chromo, function(block){
    first_SNP = block[1]
    last_SNP = block[length(block)]
    block_length = length(block)
    SNPs = paste(block, collapse = ";")
    chromo_df = data.frame(
      Block = SNPs, 
      Num_SNP = block_length,
      First_SNP = first_SNP,
      Last_SNP = last_SNP
    )
    return(chromo_df)
  })
  
  #use the map file to define start and end locations of the block
  chromo_df = left_join(chromo_df, map, c("First_SNP" = "SNP"))
  
  map = map[,c("SNP", "Position")]
  
  chromo_df = left_join(chromo_df, map, c("Last_SNP" = "SNP"))
  
  colnames(chromo_df) = c(colnames(chromo_df)[1:4], "Chrom", "Start_Pos", "End_Pos")
  
  #order the blocks based on the first SNP location
  chromo_df = chromo_df[order(chromo_df$Start_Pos), ]
  
  #Give the blocks their ID and show which chromosome they come from
  chromo_df$Block_ID = 1:nrow(chromo_df)
  chromo_df$Block_ID = paste(chromo_df$Chrom, chromo_df$Block_ID, sep = ":")
  
  return(chromo_df)
}


#overall function to return block object and apply chromosome-level condensing
block_obj_to_df = function(block_obj, map){
  #Extract first and last SNP and block markers
  block_df = map_dfr(block_obj, function(chromo){
    chromo_df = chromo_blocks_to_df(chromo, map)
  })
  
  #give each block an ID
  #block_df$Block_ID = 1:nrow(block_df)
  
  #rearrange columns
  block_df = block_df[,c(1,8,2:7)]
  
  block_df$Physical_Distance_kb = (block_df$End_Pos - block_df$Start_Pos) / 1000
  
  return(block_df)
}

