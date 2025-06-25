def_blocks = function(ld, map, method="none"){
  
  blocks_list = list()
  
  chromosomes = sort(unique(ld$Chrom))
  
  for(chr in chromosomes){print(chr)
    
    ld.chr        <- ld[ld$Chrom == chr, ]
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
      max.ld.idx <- which.max(ld.adj.chr$LD)
      max.ld     <- ld.adj.chr$LD[max.ld.idx]
      if (max.ld <= threshold) break
      
      first.mkr.blk <- as.character(ld.adj.chr$Name1[max.ld.idx])
      last.mkr.blk  <- as.character(ld.adj.chr$Name2[max.ld.idx])
      
      if (assigned[first.mkr.blk] || assigned[last.mkr.blk]) {
        ld.adj.chr <- ld.adj.chr[-max.ld.idx, ]
        next
      }
      
      block <- c(first.mkr.blk, last.mkr.blk)
      
      # Extension to the left
      repeat {
        first.idx <- match(first.mkr.blk, snps.chr)
        if (first.idx == 1) break
        
        left.mkr <- snps.chr[first.idx - 1]
        
        if (assigned[left.mkr]) break
        
        if (method=="flanking") {
          ld.vals <- ld.chr[ld.chr$Name1 == left.mkr & ld.chr$Name2 == first.mkr.blk, "LD"]
        }
        if (method=="average") {
          ld.vals <- ld.chr[ld.chr$Name1 == left.mkr & ld.chr$Name2 %in% block, "LD"]
        }
        if (mean(ld.vals$LD) <= threshold) break
        
        block <- c(left.mkr, block)
        first.mkr.blk <- left.mkr
      }
      
      # Extension to the right
      first.mkr.blk = block[1]
      repeat {
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
