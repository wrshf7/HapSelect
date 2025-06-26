####load dependencies####

library(purrr)
library(furrr)
library(future)
library(parallel)
library(dplyr)
library(progressr)

load(file = "../../Trainings/Hapotype_Stacking_Example/output/ld.R")
load(file = "../../Trainings/Hapotype_Stacking_Example/output/map.R")

ld[,1:3] = lapply(ld[,1:3], as.numeric)
ld[,4:5] = lapply(ld[,4:5], as.character)

# Hello Will


#####function to find a SNP pair#####



#####function to extend the block left####



#####function to extend the block right#####


#####master blocking function at the chromosome level#####
chromo_block = function(chr, ld, map, method, tolerance){
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
      chromo_block(chr, ld, map, method, tolerance)
      
      #after blocking on the chromosome is finished iterate the progress bar
      p()
    })
    
    
  })
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
