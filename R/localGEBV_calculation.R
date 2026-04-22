####load dependencies####

####Functions to compute the haplotype effect, PEV, and P-Value####
haplotype_effect_calc = function(block_marker_effects, haplotype, marker_means_block, center){
  #haplotype estimated effect (linear contract of marker effects)
  
  if(center){
    haplotype = haplotype - marker_means_block
  }
  
  haplotype[is.na(haplotype)] = 0
  haplotype_effect = (t(block_marker_effects) %*% haplotype)[1,1]
  return(haplotype_effect)
}

haplotype_PEV_calc = function(haplotype, haplotype_pecov, marker_means_block, center){
  if(center){
    haplotype = haplotype - marker_means_block
  }
  PEV = (t(haplotype) %*% haplotype_pecov %*% haplotype)[1,1]
  return(PEV)
}

haplotype_p_value_calc = function(haplotype_effect, haplotype_PEV){
  p_value = 2 * (1 - pnorm(abs(haplotype_effect / sqrt(haplotype_PEV))))
}

# Computes local GEBV for all individuals at a single haploblock.
# For each unique haplotype observed in the block, the effect is the dot product
# of the (optionally centred) genotype vector with the marker effects vector.
# Individual effects are then assigned by matching each individual's haplotype
# configuration back to the unique haplotype lookup table.
#
# haploblock_ID      - character ID of the block being processed (e.g. "B1")
# markers            - data frame: rows = markers in this block, col 1 = SNP name, col 2 = estimated effect
# geno_markers       - data frame: rows = markers in this block, cols = individual dosage values (0/1/2/...9)
# marker_pecov       - prediction error covariance matrix for markers; only used when haplo_test = TRUE
# haplo_test         - logical; if TRUE compute PEV and p-values per unique haplotype (requires marker_pecov)
# marker_means_block - named numeric vector of per-marker mean dosage across all individuals, used for centering
# set_missing_NA     - logical; if TRUE any individual with >= 1 missing genotype in the block gets NA effect;
#                      if FALSE only individuals with all genotypes missing get NA (partial NA imputed to 0)
# center             - logical; if TRUE subtract marker_means_block from each marker column before the multiply;
#                      must match how marker effects were originally estimated
local_GEBV_haploblock = function(haploblock_ID, markers, geno_markers, marker_pecov,
                                 haplo_test, marker_means_block, set_missing_NA, center){

  block_marker_effects = as.matrix(markers[, 2])

  geno_matrix    = t(as.matrix(geno_markers))
  haplotype_keys = apply(geno_matrix, 1, paste, collapse = ",")
  unique_keys    = unique(haplotype_keys)

  has_missing = apply(geno_matrix, 1, anyNA)
  all_missing = apply(geno_matrix, 1, function(x) all(is.na(x)))

  filtered_matrix = geno_matrix
  if (center) filtered_matrix = sweep(geno_matrix, 2, marker_means_block, "-")
  filtered_matrix[is.na(filtered_matrix)] = 0

  all_effects = drop(filtered_matrix %*% block_marker_effects)
  if (set_missing_NA) all_effects[has_missing] = NA else all_effects[all_missing] = NA

  num_haplo      = length(unique_keys)
  HaploID        = paste(haploblock_ID, seq_len(num_haplo), sep = ":")
  unique_effects = all_effects[match(unique_keys, haplotype_keys)]
  haplotype_df   = data.frame(
    Block_ID         = rep(haploblock_ID, num_haplo),
    Haplo_ID         = HaploID,
    Haplotype        = unique_keys,
    Haplotype_Effect = unique_effects,
    stringsAsFactors = FALSE
  )

  haplotype_IDs        = HaploID[match(haplotype_keys, unique_keys)]
  var_haplo_effects    = mean((all_effects    - mean(all_effects,    na.rm = TRUE))^2, na.rm = TRUE)
  unique_var_haplo_effects = mean((unique_effects - mean(unique_effects, na.rm = TRUE))^2, na.rm = TRUE)

  return_list = list(
    haplotype_df          = haplotype_df,
    haplotype_IDs         = haplotype_IDs,
    num_haplo             = num_haplo,
    var_haplo_effects     = var_haplo_effects,
    unique_var_haplo      = unique_var_haplo_effects,
    all_haplotype_effects = all_effects
  )

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
compute_local_GEBV = function(geno, marker_effects, haploblocks_df, marker_pecov, set_missing_NA = TRUE, center = TRUE){
  
  #check marker effects
  marker_effects = check_effects(marker_effects)
  
  #ensure markers are in the same order
  geno = geno[match(marker_effects[,1], geno[,1]), ]
  
  marker_means = apply(geno[,4:ncol(geno)], 1, function(x){
    row_mean = sum(x, na.rm = TRUE) / length(x[!is.na(x)])
    return(row_mean)
  })

  names(marker_means) = geno[,1]

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
      marker_means_block = marker_means[marker_ids]

      #compute all of the necessary stats and values
      local_GEBV_obj = local_GEBV_haploblock(
        haploblock_ID = haploblock,
        markers = markers,
        geno_markers = geno_markers,
        marker_pecov = marker_pecov,
        haplo_test = haplo_test,
        marker_means_block = marker_means_block,
        set_missing_NA = set_missing_NA, 
        center = center
      )
      
      #progress the progress bar
      p()
      
      return(local_GEBV_obj)
      
    })
  })
    
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
  
