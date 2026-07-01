##################################
###### Computing localGEBV ######
##################################

# Computes local GEBV for each individual at each haploblock using the genotype/dosage approach.
# Marker effects should come from a marker-level genomic prediction model (e.g. RR-BLUP / GBLUP)
# and genotypes should be dosage values (0/1/2).
# mean_adjust defaults to TRUE because marker effects are typically estimated on mean-centered
# genotypes; subtracting the per-marker mean at prediction time matches that estimation scale
# and ensures the average individual receives a local effect of zero at each block.
# Inputs:
# geno - data frame of genotypes, rows = markers, cols = individuals (starting from 4th col), first col = marker ID
# marker_effects - data frame with at least 2 columns: first col = marker ID (matching geno), second col = marker effect size
# haploblocks_df - data frame defining haploblocks, must have at least Block_ID and Block columns; Block column should list marker IDs in the block separated by ";"
# marker_pecov - optional prediction error covariance matrix for markers; required if haplo_test is TRUE
# set_missing_NA - logical; if TRUE any individual with >= 1 missing genotype in the block gets NA effect; if FALSE only individuals with all genotypes missing get NA (partial NA imputed to 0)
# mean_adjust - logical; if TRUE subtract marker means from genotypes before multiplying by effects; must match how marker effects were originally estimated
# parallel - logical; if TRUE, process haploblocks in parallel using future; if FALSE, process sequentially with progress bar
# chunk_size - integer; maximum number of blocks per chunk in parallel mode; caps the genotype data sent to each worker (default 100)
compute_local_GEBV = function(geno, marker_effects, haploblocks_df, marker_pecov,
                              set_missing_NA = TRUE, mean_adjust = TRUE,
                              parallel = FALSE, chunk_size = 100){
  .compute_local_block_effects(
    geno = geno,
    marker_effects = marker_effects,
    haploblocks_df = haploblocks_df,
    marker_pecov = marker_pecov,
    set_missing_NA = set_missing_NA,
    mean_adjust = mean_adjust,
    parallel = parallel,
    chunk_size = chunk_size
  )
}

# Computes local haplotype effects for each individual at each haploblock using the true haplotype approach.
# Each row in geno should represent a single chromosome (phased, binary 0/1 allele values) rather
# than a diploid dosage. Marker effects should be per-marker contributions derived from a
# haplotype-level prediction model.
# mean_adjust defaults to TRUE because the same logic applies as the dosage approach: effects are
# estimated on mean-centered input (allele frequencies play the role of the per-marker mean), so
# binary allele values must be centered before multiplication to avoid the population mean leaking
# into every individual's local effect.
# Inputs:
# geno - data frame of phased genotypes, rows = markers, cols = chromosomes (starting from 4th col), first col = marker ID; allele values should be binary (0/1)
# marker_effects - data frame with at least 2 columns: first col = marker ID (matching geno), second col = marker effect size derived from a haplotype-level model
# haploblocks_df - data frame defining haploblocks, must have at least Block_ID and Block columns; Block column should list marker IDs in the block separated by ";"
# marker_pecov - optional prediction error covariance matrix for markers; required if haplo_test is TRUE
# set_missing_NA - logical; if TRUE any chromosome with >= 1 missing allele in the block gets NA effect; if FALSE only chromosomes with all alleles missing get NA (partial NA imputed to 0)
# mean_adjust - logical; if TRUE subtract per-marker allele frequencies from allele values before multiplying by effects; must match how marker effects were originally estimated
# parallel - logical; if TRUE, process haploblocks in parallel using future; if FALSE, process sequentially with progress bar
# chunk_size - integer; maximum number of blocks per chunk in parallel mode; caps the genotype data sent to each worker (default 100)
compute_haplotype_effects = function(geno, marker_effects, haploblocks_df, marker_pecov,
                                     set_missing_NA = TRUE, mean_adjust = TRUE,
                                     parallel = FALSE, chunk_size = 100){
  .compute_local_block_effects(
    geno = geno,
    marker_effects = marker_effects,
    haploblocks_df = haploblocks_df,
    marker_pecov = marker_pecov,
    set_missing_NA = set_missing_NA,
    mean_adjust = mean_adjust,
    parallel = parallel,
    chunk_size = chunk_size
  )
}

# Internal computation engine. See compute_local_GEBV() or compute_haplotype_effects() for
# full parameter documentation.
.compute_local_block_effects = function(geno, marker_effects, haploblocks_df, marker_pecov, set_missing_NA = TRUE, mean_adjust = TRUE, parallel = FALSE, chunk_size = 100){

  # If users submit a custom haploblocks_df file, ensure it was properly formatted
  check_haploblocks_df(haploblocks_df)

  # Prepare inputs and pre-calculate any necessary intermediate objects for the local GEBV calculation.
  # This includes aligning marker effects with genotype markers, creating lookup tables for block markers, and determining if haplotype testing is needed based on the presence of marker_pecov.
  prep = prepare_local_gebv_inputs(
    geno = geno,
    marker_effects = marker_effects,
    haploblocks_df = haploblocks_df,
    marker_pecov = marker_pecov
  )

  # Define the block mapping function that will be applied to each haploblock.
  # This function calls the core calculation function for a single block, passing the relevant subset of markers and genotypes, as well as the pre-calculated inputs.
  block_mapper = function(haploblock, p = NULL){
    calculate_gebv_block(
      haploblock = haploblock,
      marker_idx = prep$block_marker_idx[[haploblock]],
      marker_ids = prep$marker_ids,
      effect_vec = prep$effect_vec,
      geno_matrix = prep$geno_matrix,
      marker_means = prep$marker_means,
      marker_pecov = prep$marker_pecov,
      haplo_test = prep$haplo_test,
      set_missing_NA = set_missing_NA,
      mean_adjust = mean_adjust,
      p = p
    )
  }

  # Process each haploblock to compute local GEBV
  # If in serial mode, use purrr::map with a progress bar to apply the block_mapper to each block ID in haploblocks_df.
  if(!parallel){
    progressr::handlers("txtprogressbar")
    progressr::with_progress({
      p = progressr::progressor(steps = nrow(haploblocks_df))
      local_GEBV = purrr::map(haploblocks_df$Block_ID, function(haploblock){
        block_mapper(haploblock, p = p)
      })
    })
  # Otherwise, if in parallel mode, build payloads for each chunk of haploblocks based on estimated computational cost, and use furrr::future_map to process each chunk in parallel.
  # The results from all chunks are then combined back into the final local_GEBV list.
  } else {
    # Use all cores available except for 1 to avoid overloading the system
    workers = max(1L, as.integer(future::availableCores()) - 1L)

    progressr::handlers("txtprogressbar")

    # Build payloads for each chunk of haploblocks based on estimated computational cost, and target chunk size. This will determine how the haploblocks are grouped for parallel processing.
    chunk_payloads = build_gebv_chunk_payloads(
      prep = prep,
      haploblocks_df = haploblocks_df,
      workers = workers,
      chunk_size = chunk_size
    )

    # Limit the number of workers to the number of chunks to avoid idle workers for small datasets
    # Unlikely to be an issue, but this is a safeguard against cases where chunk_size is set very high or the dataset is small, resulting in fewer chunks than available cores.
    workers = min(workers, length(chunk_payloads))
    progressr::with_progress({
      p = progressr::progressor(steps = nrow(haploblocks_df))

      # If workers is 1 or less, process the chunks sequentially without parallelization to avoid overhead.
      if(workers <= 1L){
        chunk_results = purrr::map(
          chunk_payloads,
          process_gebv_chunk,
          haplo_test = prep$haplo_test,
          set_missing_NA = set_missing_NA,
          mean_adjust = mean_adjust,
          calculate_block_fn = calculate_gebv_block,
          p = p
        )
      } else{
        # Otherwise use future to process the chunks in parallel, with the specified number of workers.
        # furrr forwards progressr signals from workers back to the main session automatically.
        future::plan(future::multisession, workers = workers)
        on.exit(future::plan(future::sequential), add = TRUE)

        chunk_results = furrr::future_map(
          chunk_payloads,
          process_gebv_chunk,
          haplo_test = prep$haplo_test,
          set_missing_NA = set_missing_NA,
          mean_adjust = mean_adjust,
          calculate_block_fn = calculate_gebv_block,
          p = p
        )
      }
    })

    # Combine the results from all chunks back into the final local_GEBV list.
    local_GEBV = purrr::flatten(chunk_results)
    local_GEBV = local_GEBV[haploblocks_df$Block_ID]
  }


  # Add summary information about the number of unique haplotypes and variance explained by haplotypes for each block back into the haploblocks_df for convenient access.

  # Number of unique haplotypes observed in each block
  haploblocks_df$Num_Uniq_Hap = purrr::map_int(local_GEBV, function(x){
    x = x[[3]]
    return(x)
  })

  # Block variance explained by haplotypes (population variance of haplotype effects)
  haploblocks_df$Block_Var = purrr::map_dbl(local_GEBV, function(x){
    x = x[[4]]
    return(x)
  })

  # Unique haplotype variance
  haploblocks_df$Unique_Haplo_Block_Var = purrr::map_dbl(local_GEBV, function(x){
    x = x[[5]]
    return(x)
  })

  # Independent haplotype effects for each block
  haploblocks_ind_haplotypes = purrr::map(local_GEBV, function(x){
    x = x[[2]]
    return(x)
  }) %>% do.call(rbind, .)

  # Set row and column names for the haplotype ID matrix, which identifies each individual's haplotype configuration at each block.
  # Rows correspond to unique haplotypes (with Block_ID as row names), and columns correspond to individuals (with individual IDs as column names).
  colnames(haploblocks_ind_haplotypes) = prep$individual_ids
  row.names(haploblocks_ind_haplotypes) = haploblocks_df$Block_ID

  # Combine the haplotype effects for each block into a single data frame, and set row names to NULL for easier downstream use.
  haplotype_effects = purrr::map(local_GEBV, function(x){
    x = x[[1]]
    return(x)
  }) %>% do.call(rbind,.)
  row.names(haplotype_effects) = NULL

  # Matrix of haplotype effects for each unique haplotype observed in each block, with row names corresponding to the unique haplotype IDs and column names corresponding to individuals.
  haplotype_effect_mat = purrr::map(local_GEBV, function(x){
    x = x[[6]]
    return(x)
  }) %>% do.call(rbind, .)

  # Set row and column names for the haplotype effect matrix, which contains the effect of each unique haplotype configuration for each individual.
  row.names(haplotype_effect_mat) = row.names(haploblocks_ind_haplotypes)
  colnames(haplotype_effect_mat) = colnames(haploblocks_ind_haplotypes)

  # Construct and return the results
  results = list(haploblocks_df, haploblocks_ind_haplotypes, haplotype_effect_mat, haplotype_effects)
  names(results) = c("Haploblocks", "Haplotype_ID_Matrix", "Haplotype_Effect_Matrix", "Haplotypes")

  return(results)
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
# mean_adjust       - logical; if TRUE subtract marker_means_block from each marker column before the multiply;
#                      must match how marker effects were originally estimated
local_GEBV_haploblock = function(haploblock_ID, markers, geno_markers, marker_pecov,
                                 haplo_test, marker_means_block, set_missing_NA, mean_adjust){
  block_marker_effects = as.matrix(markers[, 2])
  geno_matrix    = t(as.matrix(geno_markers))

  # Extract all haplotype keys for the individuals in this block, and the unique haplotypes observed
  haplotype_keys = apply(geno_matrix, 1, paste, collapse = ",")
  unique_keys    = unique(haplotype_keys)

  # Matrices that identify which haplotypes have missing genotypes, and which have all genotypes missing
  has_missing = apply(geno_matrix, 1, anyNA)
  all_missing = apply(geno_matrix, 1, function(x) all(is.na(x)))

  filtered_matrix = geno_matrix

  # If mean_adjust is TRUE, center the genotype matrix by subtracting the marker means (per marker) before multiplying by effects.
  if (mean_adjust) {
    filtered_matrix = sweep(geno_matrix, 2, marker_means_block, "-")
    # For individuals with any missing genotypes, set those genotypes to 0 to prevent calculation issues.
    filtered_matrix[is.na(filtered_matrix)] = 0
  } else {
    # Replace NAs with the corresponding column means
    col_indices = which(is.na(filtered_matrix), arr.ind = TRUE)
    filtered_matrix[col_indices] = marker_means_block[col_indices[, 2]]
  }


  # Calculate the effect of each unique haplotype
  all_effects = drop(filtered_matrix %*% block_marker_effects)

  # Add NA back to individuals with missing genotypes according to earlier has_missing and all_missing matrices.
  # If set_missing_NA is TRUE, any individual with >= 1 missing genotype gets NA
  if (set_missing_NA) {
    all_effects[has_missing] = NA
  # Otherwise, only individuals with all genotypes missing get NA (partial NA imputed to 0)
  } else {
    all_effects[all_missing] = NA
  }

  # The number of unique haplotypes observed in this block
  num_haplo = length(unique_keys)
  # Create a lookup table to match each unique haplotype back to its effect
  HaploID = paste(haploblock_ID, seq_len(num_haplo), sep = ":")
  unique_effects = all_effects[match(unique_keys, haplotype_keys)]

  # Only perform PECOV and p-value calculations if haplo_test is TRUE AND marker_pecov is provided
  if (haplo_test && !missing(marker_pecov)) {
    haplotype_pecov = marker_pecov[markers[,1], markers[,1]]

    # Create a raw matrix of haplotype configurations to be used in the PEV calculation.
    raw_matrix = geno_matrix
    raw_matrix[is.na(raw_matrix)] = 0
    unique_raw = raw_matrix[match(unique_keys, haplotype_keys), , drop = FALSE]

    # Calculate PEV for each unique haplotype
    unique_PEV  = apply(unique_raw, 1, function(h) {
      haplotype_PEV_calc(h, haplotype_pecov, marker_means_block, mean_adjust)
    })
    # Calculate p-significance for each unique haplotype
    unique_pval = haplotype_p_value_calc(unique_effects, unique_PEV)
  }

  # Create data frame with common haplotype information to be returned for this block
  haplotype_df = data.frame(
    Block_ID         = rep(haploblock_ID, num_haplo),
    Haplo_ID         = HaploID,
    Haplotype        = unique_keys,
    Haplotype_Effect = unique_effects,
    stringsAsFactors = FALSE
  )

  # If PEV and p-value were calculated, add those to the haplotype_df before returning
  if (haplo_test && !missing(marker_pecov)) {
    haplotype_df$Haplotype_PEV     = unique_PEV
    haplotype_df$Haplotype_P_Value = unique_pval
  }

  #make a vector of all individuals' haplotype IDs at this block
  haplotype_IDs = HaploID[match(haplotype_keys, unique_keys)]
  #var of all haplotyple effects for block var
  var_haplo_effects = mean((all_effects - mean(all_effects,    na.rm = TRUE))^2, na.rm = TRUE)
  #haploblock effect variance (population variance)
  unique_var_haplo_effects = mean((unique_effects - mean(unique_effects, na.rm = TRUE))^2, na.rm = TRUE)

  # Return the final constructed list of outputs for this block
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

# Generates intermediate objects for local GEBV calculation, such as aligning marker effects with genotype markers
# eg: lookup tables for block markers, and determining if haplotype testing is needed based on the presence of marker_pecov.
# Inputs
# geno - data frame of genotypes, rows = markers, cols = individuals (starting from 4th col), first col = marker ID
# marker_effects - data frame with at least 2 columns: first col = marker ID
#                  (matching geno), second col = marker effect size
# haploblocks_df - data frame defining haploblocks, must have at least Block_ID and Block columns; Block column should list marker IDs in the block separated by ";"
# marker_pecov - optional prediction error covariance matrix for markers; required if haplo_test is TRUE
prepare_local_gebv_inputs = function(geno, marker_effects, haploblocks_df, marker_pecov){
  # Check marker effects input
  marker_effects = check_effects(marker_effects)
  # Align marker_effects with geno markers, ensuring that the marker IDs match and are in the same order.
  geno = geno[match(marker_effects[,1], geno[,1]), ]

  # Create a lookup table for block markers, where each block ID maps to the indices of the markers in that block.
  geno_matrix = as.matrix(geno[,4:ncol(geno), drop = FALSE])
  storage.mode(geno_matrix) = "double"

  # Marker IDs from the genotype data
  marker_ids = geno[,1]
  rownames(geno_matrix) = marker_ids

  # Marker means across all individuals for each marker, used for centering if mean_adjust is TRUE. Named vector with marker IDs as names.
  marker_means = rowMeans(geno_matrix, na.rm = TRUE)
  names(marker_means) = marker_ids

  # Effect vector for markers, with names corresponding to marker IDs. This allows for easy lookup of effect sizes by marker ID during the block calculations.
  effect_vec = as.numeric(marker_effects[,2])
  names(effect_vec) = marker_effects[,1]

  # SNP index for quick lookup of marker positions in the genotype matrix
  snp_index = setNames(seq_len(length(marker_ids)), marker_ids)

  # Block marker IDs as a list
  block_marker_ids = setNames(
    strsplit(haploblocks_df$Block, ";"),
    haploblocks_df$Block_ID
  )
  # Block marked index where each block ID maps to the indices of the markers in that block, based on the SNP index lookup.
  block_marker_idx = lapply(block_marker_ids, function(ids){
    unname(snp_index[ids])
  })

  # Perform input checks
  if(any(vapply(block_marker_idx, anyNA, logical(1)))){
    stop("Some haploblock markers in the haploblock data frame are missing from the genotype matrix.")
  }

  effect_vec_test = lapply(block_marker_ids, function(ids){
    unname(snp_index[names(effect_vec)])
  })

  if(any(vapply(effect_vec_test, anyNA, logical(1)))){
    stop("Some haploblock markers in the haploblock data frame are missing from the marker effects input.")
  }

  if(!missing(marker_pecov) && !all(dim(marker_pecov) == length(effect_vec))) {
    stop("Marker PECOV matrix must have number of rows and columns equal to the number of markers in the effects file.")
  }

  # Determine if haplotype testing is needed based on the presence of marker_pecov and the number of columns in marker_effects.
  haplo_test = !missing(marker_pecov) && ncol(marker_effects) == 2

  # Return to constructed list of prepared inputs for local GEBV calculation
  return(list(
    geno_matrix = geno_matrix,
    marker_ids = marker_ids,
    individual_ids = colnames(geno_matrix),
    effect_vec = effect_vec,
    marker_means = marker_means,
    block_marker_idx = block_marker_idx,
    haplo_test = haplo_test,
    marker_pecov = if(haplo_test) marker_pecov else NULL
  ))
}

# Estimate the computational cost of processing a haploblock based on the number of markers in the block, the number of individuals, and whether haplotype testing is being performed.
# Used to inform the chunking strategy for parallel processing, with the goal of creating balanced chunks of haploblocks that have similar total computational cost.
estimate_block_costs = function(block_marker_idx, n_individuals, haplo_test) {
  n_markers = lengths(block_marker_idx)
  block_costs = n_markers * n_individuals

  if(haplo_test){
    block_costs = block_costs + n_markers^2
  }

  return(block_costs)
}

# Make a set of balanced chunks of haploblocks based on their estimated computational costs, using a greedy algorithm to assign haploblocks to chunks while trying to balance the total cost across chunks.
make_balanced_chunks = function(block_ids, block_costs, n_chunks) {
  n_chunks = max(1L, min(as.integer(n_chunks), length(block_ids)))

  ordered_blocks = block_ids[order(block_costs[block_ids], decreasing = TRUE)]
  chunks = vector("list", n_chunks)
  chunk_loads = numeric(n_chunks)

  for(block_id in ordered_blocks){
    chunk_id = which.min(chunk_loads)
    chunks[[chunk_id]] = c(chunks[[chunk_id]], block_id)
    chunk_loads[chunk_id] = chunk_loads[chunk_id] + block_costs[block_id]
  }

  chunks[lengths(chunks) > 0]
}

# Allocates how many chunks each block group should be split into, proportional to its total
# computational cost. Every group receives at least 1 chunk. Remaining chunks are distributed
# by cost weighting using floor + largest-remainder rounding, then each group is capped at its own
# block count. Any freed slots from capping are redistributed to the highest-cost underfilled groups.
#
# group_costs     - named numeric vector: total estimated cost per group
# group_sizes     - named integer vector: number of blocks per group (upper bound on chunks per group)
# target_n_chunks - integer: total chunks to allocate across all groups
allocate_chunk_counts = function(group_costs, group_sizes, target_n_chunks){
  n_groups = length(group_costs)
  target_n_chunks = max(n_groups, as.integer(target_n_chunks))

  chunk_counts = rep(1L, n_groups)
  remaining = target_n_chunks - n_groups

  if(remaining > 0){
    weights = group_costs / sum(group_costs)
    weights[!is.finite(weights)] = 1 / n_groups

    raw_extra = remaining * weights
    extra = floor(raw_extra)
    leftover = remaining - sum(extra)

    if(leftover > 0){
      order_idx = order(raw_extra - extra, decreasing = TRUE)
      extra[order_idx[seq_len(leftover)]] = extra[order_idx[seq_len(leftover)]] + 1L
    }

    chunk_counts = chunk_counts + extra
  }

  spare = sum(pmax(chunk_counts - group_sizes, 0L))
  chunk_counts = pmin(chunk_counts, group_sizes)

  while(spare > 0){
    underfilled = which(chunk_counts < group_sizes)
    if(length(underfilled) == 0){
      break
    }

    target_groups = underfilled[order(group_costs[underfilled], decreasing = TRUE)]
    for(group_id in target_groups){
      if(spare == 0){
        break
      }
      chunk_counts[group_id] = chunk_counts[group_id] + 1L
      spare = spare - 1L
    }
  }

  return(chunk_counts)
}

# Builds self-contained data payloads for parallel dispatch. Groups blocks by chromosome (or treats all
# blocks as one group if no Chrom column), allocates chunk counts proportionally by cost via
# allocate_chunk_counts, then packs blocks within each group into balanced chunks. For each chunk,
# the genotype matrix, effects, marker means, and (optionally) marker_pecov are subsetted to only
# the rows used by that chunk's blocks. Block marker positions are re-indexed to be local to the
# chunk's compact matrix slice, so each worker receives only the data it needs rather than a sparse
# view of the full genome-wide matrix — critical on Windows where multisession copies all data to each worker.
#
# prep           - list from prepare_local_gebv_inputs
# haploblocks_df - data frame with Block_ID and optionally Chrom
# workers        - integer: number of parallel workers; sets a floor on chunk count (workers * 2)
# chunk_size     - integer: maximum blocks per chunk; sets a floor on chunk count (ceiling(n_blocks / chunk_size))
build_gebv_chunk_payloads = function(prep, haploblocks_df, workers, chunk_size){
  block_ids = haploblocks_df$Block_ID
  block_costs = estimate_block_costs(
    block_marker_idx = prep$block_marker_idx,
    n_individuals = ncol(prep$geno_matrix),
    haplo_test = prep$haplo_test
  )
  names(block_costs) = block_ids

  # target_n_chunks is the larger of two floors:
  # - workers * 2: enough chunks for good load balancing across workers
  # - ceiling(n_blocks / chunk_size): enough chunks to keep each chunk within the memory cap
  target_n_chunks = max(
    workers * 2L,
    ceiling(length(block_ids) / max(1L, as.integer(chunk_size)))
  )
  target_n_chunks = min(length(block_ids), target_n_chunks)

  # Group blocks by chromosome so related blocks stay together; fall back to a single group if no Chrom column
  if("Chromosome" %in% names(haploblocks_df)){
    block_groups = split(haploblocks_df$Block_ID, haploblocks_df$Chromosome)
  } else{
    block_groups = list(all = haploblocks_df$Block_ID)
  }

  # Allocate chunk counts across chromosome groups proportional to cost, then pack blocks within each group
  group_costs = vapply(block_groups, function(ids) sum(block_costs[ids]), numeric(1))
  group_sizes = lengths(block_groups)
  chunk_counts = allocate_chunk_counts(group_costs, group_sizes, target_n_chunks)

  block_chunks = purrr::imap(block_groups, function(group_block_ids, group_name){
    group_index = match(group_name, names(block_groups))
    make_balanced_chunks(
      block_ids = group_block_ids,
      block_costs = block_costs,
      n_chunks = chunk_counts[group_index]
    )
  }) %>% unlist(recursive = FALSE)

  progressr::with_progress({
    p = progressr::progressor(steps = length(block_chunks), label = "Preparing chunks")
    purrr::map(block_chunks, function(chunk_block_ids){
      # Identify the union of all marker rows needed by this chunk's blocks
      chunk_marker_idx = sort(unique(unlist(prep$block_marker_idx[chunk_block_ids], use.names = FALSE)))

      # Build a local position lookup: maps global row index → position within this chunk's compact matrix
      local_index = seq_along(chunk_marker_idx)
      names(local_index) = as.character(chunk_marker_idx)

      # Re-express each block's marker positions as local indices into the chunk matrix
      block_local_idx = lapply(prep$block_marker_idx[chunk_block_ids], function(idx){
        unname(local_index[as.character(idx)])
      })

      payload = list(
        block_ids = chunk_block_ids,
        block_local_idx = block_local_idx,
        chunk_marker_ids = prep$marker_ids[chunk_marker_idx],
        geno_chunk = prep$geno_matrix[chunk_marker_idx, , drop = FALSE],
        effect_chunk = prep$effect_vec[chunk_marker_idx],
        mean_chunk = prep$marker_means[chunk_marker_idx]
      )

      if(prep$haplo_test){
        payload$pecov_chunk = prep$marker_pecov[chunk_marker_idx, chunk_marker_idx, drop = FALSE]
      }

      p()
      payload
    })
  })
}

# Processes one chunk payload on a single worker. Iterates over the blocks assigned to the chunk,
# calling calculate_block_fn for each using the chunk-local genotype data. Returns a named list
# of per-block results keyed by block ID.
#
# payload            - list from build_gebv_chunk_payloads: contains block_ids, block_local_idx,
#                      chunk_marker_ids, geno_chunk, effect_chunk, mean_chunk, and optionally pecov_chunk
# haplo_test         - logical; passed through to calculate_block_fn
# set_missing_NA     - logical; passed through to calculate_block_fn
# mean_adjust        - logical; passed through to calculate_block_fn
# calculate_block_fn - function used to compute a single block (calculate_gebv_block)
# p                  - optional progressr progressor; called once per block to report progress
process_gebv_chunk = function(payload, haplo_test, set_missing_NA, mean_adjust, calculate_block_fn, p = NULL){
  chunk_results = purrr::map2(payload$block_ids, payload$block_local_idx, function(block_id, local_idx){
    calculate_block_fn(
      haploblock = block_id,
      marker_idx = local_idx,
      marker_ids = payload$chunk_marker_ids,
      effect_vec = payload$effect_chunk,
      geno_matrix = payload$geno_chunk,
      marker_means = payload$mean_chunk,
      marker_pecov = if(haplo_test) payload$pecov_chunk else NULL,
      haplo_test = haplo_test,
      set_missing_NA = set_missing_NA,
      mean_adjust = mean_adjust,
      p = p
    )
  })

  names(chunk_results) = payload$block_ids
  return(chunk_results)
}

# Calculate an individual's local GEBV at a single haploblock by matching their haplotype configuration to the unique haplotype effects calculated in local_GEBV_haploblock.
# haploblock - character ID of the block being processed (e.g. "B1")
# marker_idx - integer vector of marker indices for this block (relative to the full genotype matrix)
# marker_ids - character vector of all marker IDs (same order as rows of geno_matrix)
# effect_vec - named numeric vector of marker effects, names are marker IDs
# geno_matrix - numeric matrix of genotypes, rows are markers (same order as marker_ids), cols are individuals
# marker_means - named numeric vector of mean dosage per marker, names are marker IDs
# marker_pecov - numeric matrix of marker prediction error covariances, only needed if haplo_test = TRUE
# haplo_test - logical; if TRUE compute PEV and p-values per unique haplotype (requires marker_pecov)
# set_missing_NA - logical; if TRUE any individual with >= 1 missing genotype in the block gets NA effect;
#                  if FALSE only individuals with all genotypes missing get NA (partial NA imputed to 0)
# mean_adjust - logical; if TRUE subtract marker_means_block from each marker column before the multiply
# p - optional progress function to call after processing this block
calculate_gebv_block = function(haploblock, marker_idx, marker_ids, effect_vec,
                                geno_matrix, marker_means, marker_pecov, haplo_test,
                                set_missing_NA, mean_adjust, p = NULL){

  # Get the marker IDs for this block based on the provided marker indices, and construct a data frame of markers and their effects for this block.
  block_marker_ids = marker_ids[marker_idx]
  markers = data.frame(
    SNP = block_marker_ids,
    Effect = unname(effect_vec[marker_idx]),
    stringsAsFactors = FALSE
  )

  geno_markers = geno_matrix[marker_idx, , drop = FALSE]
  marker_means_block = marker_means[marker_idx]

  # If haplo_test is TRUE, subset the marker_pecov matrix to only the markers in this block to save memory and computation in the local_GEBV_haploblock function.
  if(haplo_test){
    marker_pecov = marker_pecov[marker_idx, marker_idx, drop = FALSE]
    rownames(marker_pecov) = block_marker_ids
    colnames(marker_pecov) = block_marker_ids
  }

  # Calculate local GEBV for this block by calling the local_GEBV_haploblock function
  local_GEBV_obj = local_GEBV_haploblock(
    haploblock_ID = haploblock,
    markers = markers,
    geno_markers = geno_markers,
    marker_pecov = marker_pecov,
    haplo_test = haplo_test,
    marker_means_block = marker_means_block,
    set_missing_NA = set_missing_NA,
    mean_adjust = mean_adjust
  )

  # Call the progress function if provided (e.g. to update a progress bar in the caller)
  if(!is.null(p)){
    p()
  }

  return(local_GEBV_obj)
}

##################################
##### Supporting Functions #######
##################################
# Check the marker effects input for proper formatting and types, and coerce as needed. The first column should be SNP IDs (characters) and the second column should be numeric effects. If coercion is needed, a warning is issued.
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


##################################
#### Supporting Calculations #####
##################################

# Center genotype matrix by subtracting the mean dosage per marker from each marker column.
center_genotypes = function(geno){
  #code from the localGEBV functions
  marker_means = rowMeans(geno[,4:ncol(geno)], na.rm = TRUE)

  geno[,4:ncol(geno)] = sweep(geno[,4:ncol(geno)], 1, marker_means, "-")

  return(geno)
}

# Calculate prediction error variance for a haplotype
haplotype_PEV_calc = function(haplotype, haplotype_pecov, marker_means_block, mean_adjust){
  if(mean_adjust){
    haplotype = haplotype - marker_means_block
  }
  PEV = drop(t(haplotype) %*% haplotype_pecov %*% haplotype)
  return(PEV)
}

# Calculate p-value for a haplotype effect based on its effect size and PEV.
haplotype_p_value_calc = function(haplotype_effect, haplotype_PEV){
  p_value = 2 * (1 - pnorm(abs(haplotype_effect / sqrt(haplotype_PEV))))
}

# Check that the haploblock data frame contains all necessary information. Sometimes users may compute their own.
check_haploblocks_df = function(haploblocks_df){
  if(!all(c("Block", "Block_ID", "Chrom") %in% colnames(haploblocks_df))){
    stop("Please ensure the haploblock data frame contains the 'Block', 'Block_ID', and 'Chrom' Columns")
  }

  if(!is.character(haploblocks_df$Block) || anyNA(haploblocks_df$Block) ){
    stop("Please ensure each row of the \"Block\" column contains the markers corresponding to a given block ID and are separated by a ';'. Please also ensure that the column is not a factor.")
  }

  if(!any(grepl(";", haploblocks_df$Block))){
    warning("Did not detect the ';' separator in any blocks. Are all blocks single marker blocks? If not, please ensure you use the exact separator, otherwise extracting marker IDs for blocks will fail!")
  }

  split_blocks = strsplit(haploblocks_df$Block, ";")

  if(any(vapply(split_blocks, function(x) any(x == ""), logical(1)))){
    stop("Block strings contain empty marker IDs (e.g., trailing ';', leading ';', or double ';;')")
  }

}
