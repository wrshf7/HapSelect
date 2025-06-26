####load dependencies####

#notes - check returns, assign values, etc.

library(purrr)
library(furrr)
library(future)
library(parallel)
library(dplyr)
library(progressr)

load(file = "Example_Files/ld.R")
load(file = "Example_Files/map.R")

ld[,1:3] = lapply(ld[,1:3], as.numeric)
ld[,4:5] = lapply(ld[,4:5], as.character)



#####function to extend the block left####

extend_left = function(first.mkr.blk, snps.chr, assigned, ld.chr, block, method){
  
  #repeat left extension until break point is reached (tolerance and threshold)
  repeat {
    
    #extract index of the first marker in the block from the SNP pair
    first.idx <- match(first.mkr.blk, snps.chr)
    
    #don't extend left if it was the first marker
    if (first.idx == 1) break
    
    #get index of marker to the left of the first SNP in the block
    left.mkr <- snps.chr[first.idx - 1]
    
    #if marker to the left is already in a block, break
    if (assigned[left.mkr]) break
    
    #if utilizing the flanking method, pull the LD of the SNP to the left idetnified from left.mkr
    if (method=="flanking") {
      #singular value?
      ld.vals <- ld.chr[ld.chr$Name1 == left.mkr & ld.chr$Name2 == first.mkr.blk, "LD"]
    }
    
    #if utilizing the average method, then compute the average LD of the new SNP considered with all markers in the block currently
    if (method=="average") {
      #could be multiple values if the block is greater than
      ld.vals <- ld.chr[ld.chr$Name1 == left.mkr & ld.chr$Name2 %in% block, "LD"]
    }
    
    #check that the LD value (flanking) or average LD value (average) meets the threshold - mean of a single marker is itself, so it works for both
    if (mean(ld.vals$LD) <= threshold) break
    
    #if the threshold is met, add it to the block
    block <- c(left.mkr, block)
    
    #update the first marker and repeat the loop
    #first.mkr.blk <- left.mkr
  }
  
  
  #return relevant objects
  return(block)
}

#####function to extend the block right#####

extend_right = function(last.mkr.blk, snps.chr, assigned, ld.chr, block, method){
  
  #keep extending until the break signal is initiated
  repeat{
    #starting from the second marker of the original pair (last marker)
    last.idx <- match(last.mkr.blk, snps.chr)
    
    #if it's the last marker in the chromo, don't extend
    if (last.idx == length(snps.chr)) break
    
    #if not, find the next marker
    right.mkr <- snps.chr[last.idx + 1]
    
    #check if the next marker is already in a block
    if (assigned[right.mkr]) break
    
    #same as left extension - pull the LD of the marker to the right or the average of the new marker with the rest of the markers in the block originating from the original pair
    if (method=="flanking") {
      ld.vals <- ld.chr[ld.chr$Name1 == last.mkr.blk & ld.chr$Name2 == right.mkr, "LD"]
    }
    
    if (method=="average") {
      ld.vals <- ld.chr[ld.chr$Name1 %in% block & ld.chr$Name2 == right.mkr, "LD"]
    }
    
    #check to see if it meets threshold criteria
    if (mean(ld.vals$LD) <= threshold) break
    
    #add it to the block
    block <- c(block, right.mkr)
    
    #update the last marker
    #last.mkr.blk <- right.mkr
  }
  return(block)
}


#####function to find a SNP pair to start block extension#####

make_block = function(ld.chr, ld.adj.char, snps.chr, snps.pos.chr, snps.chr, first.mkr.chr, last.mkr.chr, assigned, threshold, tolerance){
  #keep making blocks until there is no pair with high enough LD
  repeat {
    
    #identify max LD pair and pull its id
    max.ld.idx <- which.max(ld.adj.chr$LD)
    
    #extract that LD value and check it
    max.ld     <- ld.adj.chr$LD[max.ld.idx]
    if (max.ld <= threshold) break
    
    #Pull the two marker names in the LD pair
    first.mkr.blk <- as.character(ld.adj.chr$Name1[max.ld.idx])
    last.mkr.blk  <- as.character(ld.adj.chr$Name2[max.ld.idx])
    
    #if assigned is true for either marker in the LD pair (assigned to a block already), remove it from consideration - updated at the end of the block formation
    if (assigned[first.mkr.blk] || assigned[last.mkr.blk]) {
      ld.adj.chr <- ld.adj.chr[-max.ld.idx, ]
      next
    }
    
    #extend the block
    block <- c(first.mkr.blk, last.mkr.blk)
    
    #extend the block
    extend_left()
    extend_right()
    
    #update the assigned variable to indicate that SNP have been assigned to a block
    assigned[block] = TRUE
    
    
    #blocks_list[[length(blocks_list) + 1]] <- block
    
    
    #ask victor about this code
    used = with(ld.adj.chr, Name1 %in% block | Name2 %in% block)
    ld.adj.chr = ld.adj.chr[!used, ]
    
    
    if (nrow(ld.adj.chr) == 0) break
    
    
  }
}


#####master blocking function at the chromosome level#####
chromo_blocking = function(chr, ld, map, method, tolerance){
  #pull the LD information for the chromosome
  ld.chr = ld[ld$Chrom == chr, ]
  
  #pulling the immediately adjacent SNP to find all SNP pairs
  ld.adj.chr = ld.chr[ld.chr$Locus2 == ld.chr$Locus1 + 1, ]
  
  #pull all snp name on the chromosome
  snps.chr = unique(c(ld.chr$Name1, ld.chr$Name2))
  
  #match snp name to the position in the map file
  snps.pos.chr = map$pos[match(snps.chr, map$SNP)]
  
  #order the snp based on position - should already be ordered
  snps.chr = snps.chr[order(snps.pos.chr)]
  
  #find the first and last SNP on the chromo - utilized to check whether blocks have reached the end of the chromosome 
  first.mkr.chr = snps.chr[1]
  last.mkr.chr = snps.chr[length(snps.chr)]
  
  #tracker utilized later to determine whether a SNP has been added to a block
  assigned = rep(FALSE, length(snps.chr))
  names(assigned) <- snps.chr
  
  make_block()

  unassigned = names(assigned)[!assigned]
  for (snp in unassigned) {
    blocks_list[[length(blocks_list) + 1]] = snp
  }
  
  return(blocks_list)
}
  





#####overall function to track progress, take input, and call blocking functions#####
def_blocks = function(ld, map, method = "flanking", tolerance = 1){
  
  blocks_list = list()
  #pull the list of chromosomes to parallelize
  chromosomes = sort(unique(ld$Chrom))
  
  #setup parallelization using futures and parallel package and utilize all but 1 core
  futures::plan(multisession, workers = parallel::detectCores() - 1)
  
  #track progress across chromosomes and set up progress bar
  handlers("txtprogressbar")
  with_progress({
    
    #while tracking, set up the progressor for when the progress bar advances
    p = progressor(steps = length(chromosomes))
    
    #form blocks across chromosomes - parallelize each chromosome
    furrr::map(chromosomes, function(chr){
      
      #call chromosome blocking function
      chromo_blocking(chr, ld, map, method, tolerance)
      
      #after blocking on the chromosome is finished iterate the progress bar
      p()
    })
  })
  
  max_len = max(lengths(blocks_list))
  blocks_df = as.data.frame(do.call(rbind, lapply(blocks_list, function(x) {
    length(x) = max_len  # remplit avec NA
    x
  })), stringsAsFactors = FALSE)
  
  blocks_df$chr = map$chrom[match(blocks_df$V1, map$SNP)]
  blocks_df$pos = map$pos[match(blocks_df$V1, map$SNP)]
  blocks_df = blocks_df[order(blocks_df$pos),]
  
  return(blocks_df)
}



def_blocks = function(ld, map, method="none"){
  
  blocks_list = list()
  
  chromosomes = sort(unique(ld$Chrom))
  
  
  
  
  for(chr in chromosomes){print(chr)
    
    ld.chr        <- ld[ld$Chrom == chr, ]
    
    #doesn't this only consider the immediately adjacent snp? Are blocks being extended based on LD between the new SNP and the original SNP pair (depending on if you're moving left or right) 
    #or the SNP last incorporated?
    #I think it needs to be based on the original SNP pair
    ld.adj.chr    <- ld.chr[ld.chr$Locus2 == ld.chr$Locus1 + 1, ]
    snps.chr      <- unique(c(ld.chr$Name1, ld.chr$Name2))
    snps.pos.chr  <- map$pos[match(snps.chr, map$SNP)]
    snps.chr      <- snps.chr[order(snps.pos.chr)]
    
    first.mkr.chr <- snps.chr[1]
    last.mkr.chr  <- snps.chr[length(snps.chr)]
    
    assigned <- rep(FALSE, length(snps.chr))
    names(assigned) <- snps.chr
    
    
    ## BUILD THE HAPLOBLOCKS
    ########################
    repeat {
     
      #is this necessary to identify if there is even a pair to advance from first?
      
      #identify max LD between a SNP and the one immediately following?
      #need to consider moving backwards as well?
      max.ld.idx <- which.max(ld.adj.chr$LD)
      
      #extract that LD value and check it
      max.ld     <- ld.adj.chr$LD[max.ld.idx]
      if (max.ld <= threshold) break
      
      #Pull the two marker names in the LD pair
      first.mkr.blk <- as.character(ld.adj.chr$Name1[max.ld.idx])
      last.mkr.blk  <- as.character(ld.adj.chr$Name2[max.ld.idx])
      
      #if assigned is true for either marker in the LD pair, remove it from consideration, but when is assigned updated?
      if (assigned[first.mkr.blk] || assigned[last.mkr.blk]) {
        ld.adj.chr <- ld.adj.chr[-max.ld.idx, ]
        next
      }
      
      #extend the block
      block <- c(first.mkr.blk, last.mkr.blk)
      
      # Extension to the left
      #I think this answers my question about extending to the left
      repeat {
        first.idx <- match(first.mkr.blk, snps.chr)
        if (first.idx == 1) break
        
        left.mkr <- snps.chr[first.idx - 1]
        
        if (assigned[left.mkr]) break
        
        if (method=="flanking") {
          #singular value?
          ld.vals <- ld.chr[ld.chr$Name1 == left.mkr & ld.chr$Name2 == first.mkr.blk, "LD"]
        }
        if (method=="average") {
          #could be multiple values?
          ld.vals <- ld.chr[ld.chr$Name1 == left.mkr & ld.chr$Name2 %in% block, "LD"]
        }
        if (mean(ld.vals$LD) <= threshold) break
        
        block <- c(left.mkr, block)
        first.mkr.blk <- left.mkr
      }
      
      # Extension to the right
      
      first.mkr.blk = block[1]
      repeat {
        
        #starting from the second marker of the original pair
        last.idx <- match(last.mkr.blk, snps.chr)
        if (last.idx == length(snps.chr)) break
        
        right.mkr <- snps.chr[last.idx + 1]
        
        if (assigned[right.mkr]) break
        
        if (method=="flanking") {
          ld.vals <- ld.chr[ld.chr$Name1 == last.mkr.blk & ld.chr$Name2 == right.mkr, "LD"]
        }
        if (method=="average") {
          ld.vals <- ld.chr[ld.chr$Name1 %in% block & ld.chr$Name2 == right.mkr, "LD"]
        }
        
        
        if (mean(ld.vals$LD) <= threshold) break
        
        block <- c(block, right.mkr)
        last.mkr.blk <- right.mkr
      }
      
      #start here
      assigned[block] <- TRUE
      blocks_list[[length(blocks_list) + 1]] <- block
      
      used <- with(ld.adj.chr, Name1 %in% block | Name2 %in% block)
      ld.adj.chr <- ld.adj.chr[!used, ]
      if (nrow(ld.adj.chr) == 0) break
    }
    
    
    ## COMBINE RESULTS 
    ##################
    
    unassigned <- names(assigned)[!assigned]
    for (snp in unassigned) {
      blocks_list[[length(blocks_list) + 1]] <- snp
    }
    
    
  }
  
  max_len <- max(lengths(blocks_list))
  blocks_df <- as.data.frame(do.call(rbind, lapply(blocks_list, function(x) {
    length(x) <- max_len  # remplit avec NA
    x
  })), stringsAsFactors = FALSE)
  
  blocks_df$chr = map$chrom[match(blocks_df$V1, map$SNP)]
  blocks_df$pos = map$pos[match(blocks_df$V1, map$SNP)]
  blocks_df = blocks_df[order(blocks_df$pos),]
  
  
  return(blocks_df)
}
