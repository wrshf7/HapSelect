####load dependencies####

library(purrr)
#library(furrr)
#library(future)
#library(parallel)
library(dplyr)
library(progressr)




####Functions to compute the haplotype effect, PEV, and P-Value####
haplotype_effect_calc = function(block_marker_effects, haplotype){
  #haplotype estimated effect (linear contract of marker effects)
  haplotype_effect = (t(block_marker_effects) %*% haplotype)[1,1]
  return(haplotype_effect)
}

haplotype_PEV_calc = function(haplotype, haplotype_pecov){
  PEV = (t(haplotype) %*% haplotype_pecov %*% haplotype)[1,1]
  return(PEV)
}

haplotype_p_value_calc = function(haplotype_effect, haplotype_PEV){
  p_value = 2 * (1 - pnorm(abs(haplotype_effect / sqrt(haplotype_PEV))))
}

####Function to compute local GEBV for all individuals at a given haploblock, define unique haplotypes, and compute other metrics######
#unique haplotypes
local_GEBV_haploblock = function(haploblock_ID, markers, geno_markers, marker_pecov, haplo_test){

  #extract all haplotypes and unique haplotypes
  haplotypes = purrr::map(geno_markers, function(haplotype){
    haplotype = paste(haplotype, collapse = "")
  }) %>% do.call(c,.)
  
  unique_haplotypes = unique(haplotypes)
  unique_haplotypes = unique_haplotypes[!grepl("NA", unique_haplotypes)]

  
  
  #haploblock markers
  block_marker_effects = as.matrix(markers[,2])
  
  
  #vector of unique haplotype effects, haplotype variances, and p-values
  if(!haplo_test){
    
    haplotype_effects = purrr::map(unique_haplotypes, function(haplotype){
      haplotype = as.numeric(strsplit(haplotype, "")[[1]])
      
      #haplotype estimated effect (linear contract of marker effects)
      haplotype_effect = haplotype_effect_calc(block_marker_effects = block_marker_effects, haplotype = haplotype)
      
      return(haplotype_effect)
      
    }) %>% do.call(c,.)
    
    haplotype_effects = data.frame(
      Haplotype_Effect = haplotype_effects
    )
    
  } else if(haplo_test){
    
    haplotype_effects = purrr::map_dfr(unique_haplotypes, function(haplotype){
      haplotype = as.numeric(strsplit(haplotype, "")[[1]])
      
      #haplotype estimated effect (linear contract of marker effects)
      haplotype_effect = haplotype_effect_calc(block_marker_effects = block_marker_effects, haplotype = haplotype)
      
      haplotype_pecov = marker_pecov[which(colnames(marker_pecov) %in% markers$SNP), which(colnames(marker_pecov) %in% markers$SNP)]
    
      #compute the PEV of the haplotype
      haplotype_variance = haplotype_PEV_calc(haplotype = haplotype, 
                                              haplotype_pecov = haplotype_pecov)
    
      #compute the p-value
      haplotype_p_value = haplotype_p_value_calc(haplotype_effect = haplotype_effect, 
                                                 haplotype_PEV = haplotype_variance)
    
      return_df = data.frame(
        Haplotype_Effect = haplotype_effect,
        Haplotype_PEV = haplotype_variance,
        Haplotype_Effect_p_value = haplotype_p_value
      )
    
      return(return_df)
    })
  }
  
  

  #haploblock effect variance (population variance)
  unique_var_haplo_effects = mean((haplotype_effects$Haplotype_Effect - mean(haplotype_effects$Haplotype_Effect))^2)
  
  #number of haplotypes within a block
  num_haplo = length(unique_haplotypes)
  
  #return all of the values
  HaploID = paste(haploblock_ID, 1:num_haplo, sep = ":")
  haplotype_df = data.frame(
    Block_ID = rep(haploblock_ID, num_haplo),
    Haplo_ID = HaploID,
    Haplotype = unique_haplotypes
  )
  
  haplotype_df = cbind(haplotype_df, haplotype_effects)
  
  #make a vector of all individuals' haplotype IDs at this block
  haplotype_IDs = HaploID[match(haplotypes, unique_haplotypes)]
  
  #vector of all haplotype effects
  all_haplotype_effects = haplotype_df[match(haplotype_IDs, haplotype_df$Haplo_ID),"Haplotype_Effect"]
  
  #var of all haplotyple effects for block var
  var_haplo_effects = mean((all_haplotype_effects - mean(all_haplotype_effects, na.rm = TRUE))^2, na.rm = TRUE)
  
  #put all values in a list and return them
  return_list = list(haplotype_df, haplotype_IDs, num_haplo, var_haplo_effects, unique_var_haplo_effects, all_haplotype_effects)
  
  return(return_list)
}

####check marker effects file####
check_effects = function(marker_effects){
  if(!is.character(marker_effects[,1])){
    marker_effects[,1] = as.character(marker_effects[,1])
    warning("Coerced SNP names to characters.")
  }
  
  marker_effects = purrr::imap_dfc(marker_effects, function(trait, index){
    if(colnames(marker_effects)[1] != index){
      if(!is.numeric(marker_effects[,index])){
        warning(paste0("Coercing marker effects for trait ", index, " to numeric."))
      }
      trait = as.numeric(trait)
    }
    return(trait)
  }) %>% as.data.frame()
  
  
  return(marker_effects)
}

####Head function to split by block and compute local GEBV per individual####
compute_local_GEBV = function(geno, marker_effects, haploblocks_df, marker_pecov){
  
  #check marker effects
  marker_effects = check_effects(marker_effects)
  
  #ensure markers are in the same order
  geno = geno[match(marker_effects[,1], geno[,1]), ]
  
  #set up
  #future::plan(multisession, workers = parallel::detectCores() - 1)
  
  #check if marker_pecov is provided
  if(!missing(marker_pecov) & ncol(marker_effects) == 2){
    haplo_test = TRUE
  } else{
    haplo_test = FALSE
  }
  
  #lookup index
  snp_index = setNames(seq_len(nrow(geno)), geno[,1])
  marker_index = setNames(seq_len(nrow(marker_effects)), marker_effects[,1])
  
  
  #set up progress bar
  progressr::handlers("txtprogressbar")
  progressr::with_progress({
    #define number of steps in progress bar
    p = progressr::progressor(steps = nrow(haploblocks_df))
    
    #spread the CPU love across the different haploblocks
    local_GEBV = purrr::map(haploblocks_df$Block_ID, function(haploblock){
      
      #extract the markers in the specific block ID
      #markers = strsplit(haploblocks_df[haploblocks_df$Block_ID == haploblock, "Block"], ";")[[1]]
      #markers = marker_effects[marker_effects[,1] %in% markers, ]
      
      marker_ids = strsplit(haploblocks_df[haploblocks_df$Block_ID == haploblock, "Block"], ";")[[1]]
      markers = marker_effects[marker_index[marker_ids], ]
      
      #extract genotypes of the markers
      geno_markers = geno[geno[,1] %in% markers[,1], 4:ncol(geno)]
      
      #compute all of the necessary stats and values
      if(haplo_test){
        local_GEBV_obj = local_GEBV_haploblock(haploblock_ID = haploblock, 
                                               markers = markers, 
                                               geno_markers = geno_markers, 
                                               marker_pecov = marker_pecov,
                                               haplo_test = TRUE)
      }else{
        local_GEBV_obj = local_GEBV_haploblock(haploblock_ID = haploblock, 
                                               markers = markers, 
                                               geno_markers = geno_markers,
                                               haplo_test = FALSE)
      }
      
      #progress the progress bar
      p()
      
      return(local_GEBV_obj)
      
    })
  })
  
  cat("Formatting Data\n")
  
  #extract number of haplotypes per block
  haploblocks_df$Num_Uniq_Hap = purrr::map_int(local_GEBV, function(x){
    x = x[[3]]
    return(x)
  })
  
  
  #extract variance of haplotype effects within a block
  haploblocks_df$Block_Var = purrr::map_dbl(local_GEBV, function(x){
    x = x[[4]]
    return(x)
  })
  
  #extract variance of unique haplotype effects within a block
  haploblocks_df$Unique_Haplo_Block_Var = purrr::map_dbl(local_GEBV, function(x){
    x = x[[5]]
    return(x)
  })
  
  #extract haplotype IDs of all individuals for all blocks
  haploblocks_ind_haplotypes = purrr::map(local_GEBV, function(x){
    x = x[[2]]
    return(x)
  }) %>% do.call(rbind, .) %>% as.data.frame(.)
  
  colnames(haploblocks_ind_haplotypes) = colnames(geno[,4:ncol(geno)])
  row.names(haploblocks_ind_haplotypes) = haploblocks_df$Block_ID
  
  #extract haplotype IDs, their effects, PEV, and p-values
  haplotype_effects = purrr::map(local_GEBV, function(x){
    x = x[[1]]
    return(x)
  }) %>% do.call(rbind,.)
  
  #extract haplotype effect matrix
  haplotype_effect_mat = purrr::map(local_GEBV, function(x){
    x = x[[6]]
    return(x)
  }) %>% do.call(rbind, .) %>% as.data.frame(.)
  

  row.names(haplotype_effect_mat) = row.names(haploblocks_ind_haplotypes)
  colnames(haplotype_effect_mat) = colnames(haploblocks_ind_haplotypes)
  
  return_obj = list(haploblocks_df, haploblocks_ind_haplotypes, haplotype_effect_mat, haplotype_effects)
  
  names(return_obj) = c("Haploblocks", "Haplotype_ID_Matrix", "Haplotype_Effect_Matrix", "Haplotypes")
  
  return(return_obj)
}
  
