##################################
###### Building Haploblocks ######
##################################


# extend_block -----------------------------------------------------------------
# NOTE: This is an R implementation and it is much slower than the C++ version in extend_block.cpp. (3.5× to 17.5× depending on method and settings)
#       It is retained here for clarity and testing purposes.
# markers one at a time until the LD threshold is no longer met.
#
# direction    : -1 to extend left, +1 to extend right
# edge_marker  : the marker currently at the boundary of the block in the
#                direction of extension
# marker_names : ordered list of marker names for the chromosome,
#                sorted by physical position
# marker_idx   : named integer lookup of marker name -> position index,
#                pre-built in chromo_blocking for O(1) position lookups
# assigned     : lookup table of marker name -> TRUE/FALSE indicating whether
#                each marker has already been placed in a block
# ld_lookup    : named numeric lookup of "Name1,Name2" -> LD value,
#                pre-built in chromo_blocking for O(1) LD lookups
# block        : marker names currently in the block
# method       : "flanking" — LD of candidate vs edge marker only;
#                "average"  — mean LD of candidate vs all markers in the block
# threshold    : minimum LD value required to extend the block
# tolerance    : number of consecutive below-threshold markers allowed before
#                stopping (absorbed into the block if a later marker meets the threshold)
# tol_reset    : if TRUE, reset the tolerance counter after each accepted marker
extend_block = function(direction, edge_marker, marker_names, marker_idx,
                        assigned, ld_lookup, block, method, threshold,
                        tolerance, tol_reset) {

  tolerance_counter = 0
  failed_markers    = c()

  repeat {

    edge_index = marker_idx[[edge_marker]]

    if (direction == -1 && edge_index == 1)                  break
    if (direction ==  1 && edge_index == length(marker_names))  break

    step = length(failed_markers) + 1
    if (direction == -1) {
      candidate_marker = marker_names[edge_index - step]
    } else {
      candidate_marker = marker_names[edge_index + step]
    }

    if (length(candidate_marker) == 0 || is.na(candidate_marker)) break
    if (assigned[candidate_marker])                             break

    if (method == "flanking") {
      if (direction == -1) {
        # candidate is to the left so Name1 = candidate, Name2 = edge
        ld_vals = ld_lookup[paste(candidate_marker, edge_marker, sep = ",")]
      } else {
        # candidate is to the right so Name1 = edge, Name2 = candidate
        ld_vals = ld_lookup[paste(edge_marker, candidate_marker, sep = ",")]
      }
    }

    if (method == "average") {
      if (direction == -1) {
        # candidate is to the left so Name1 = candidate, Name2 = any block SNP
        ld_vals = ld_lookup[paste(candidate_marker, block, sep = ",")]
      } else {
        # candidate is to the right so Name1 = any block SNP, Name2 = candidate
        ld_vals = ld_lookup[paste(block, candidate_marker, sep = ",")]
      }
    }

    ld_mean = mean(ld_vals, na.rm = TRUE)

    if (ld_mean <= threshold || is.nan(ld_mean)) {
      if (tolerance_counter >= tolerance) {
        break
      } else {
        tolerance_counter = tolerance_counter + 1
        failed_markers = c(candidate_marker, failed_markers)
        next
      }
    }

    if (tol_reset) tolerance_counter = 0

    if (length(failed_markers) > 0) {
      if (direction == -1) {
        block = c(candidate_marker, failed_markers, block)
      } else {
        block = c(block, failed_markers, candidate_marker)
      }
      failed_markers = c()
    } else {
      if (direction == -1) {
        block = c(candidate_marker, block)
      } else {
        block = c(block, candidate_marker)
      }
    }

    edge_marker = candidate_marker
  }

  return(block)
}


# make_blocks ------------------------------------------------------------------
# NOTE: This is an R implementation and it is much slower than the C++ version in make_blocks.cpp (3.5× to 17.5× depending on method and settings)
#       It is retained here for clarity and testing purposes.
# Drives block formation for a single chromosome. Finds seed marker pairs and
# calls extend_block to grow each block left and right. Any markers not absorbed
# into a block remain unassigned and are handled by the caller.
#
# ld_lookup        : named numeric lookup of "Name1,Name2" -> LD value,
#                    pre-built in chromo_blocking for O(1) LD lookups
# ld_adj           : adjacent-only subset of the chromosome LD table
#                    (Locus2 == Locus1 + 1); used to find seed pairs
# marker_names     : ordered list of marker names for the chromosome,
#                    sorted by physical position
# marker_idx       : named integer lookup of marker name -> position index,
#                    pre-built in chromo_blocking for O(1) position lookups
# marker_positions : physical positions of each marker, matching the order of marker_names
# first_marker     : name of the first (leftmost) marker on the chromosome
# last_marker      : name of the last (rightmost) marker on the chromosome
# assigned         : lookup table of marker name -> TRUE/FALSE indicating whether
#                    each marker is already in a block
# method           : LD evaluation method passed through to extend_block
# threshold        : minimum LD value required to seed or extend a block
# tolerance        : consecutive below-threshold markers allowed during extension
# tol_reset        : whether to reset the tolerance counter on each accepted marker
# start            : "LD"        — seed from the highest-LD adjacent pair first;
#                    "beginning" — sweep left to right from the first marker
make_blocks = function(ld_lookup, ld_adj, marker_names, marker_idx,
                       marker_positions, first_marker, last_marker, assigned,
                       method, threshold, tolerance, tol_reset, start) {

  chrom_blocks = list()

  if (start == "LD") {

    repeat {

      max_ld_idx = which.max(ld_adj$LD)
      max_ld     = ld_adj$LD[max_ld_idx]
      if (max_ld <= threshold) break

      seed_marker_left  = as.character(ld_adj$Name1[max_ld_idx])
      seed_marker_right = as.character(ld_adj$Name2[max_ld_idx])

      # if either seed marker is already in a block, drop this pair and move on
      if (assigned[seed_marker_left] || assigned[seed_marker_right]) {
        ld_adj = ld_adj[-max_ld_idx, ]
        next
      }

      block = c(seed_marker_left, seed_marker_right)

      block = extend_block(
        direction    = -1,
        edge_marker  = seed_marker_left,
        marker_names = marker_names,
        marker_idx   = marker_idx,
        assigned     = assigned,
        ld_lookup    = ld_lookup,
        block        = block,
        method       = method,
        threshold    = threshold,
        tolerance    = tolerance,
        tol_reset    = tol_reset
      )

      block = extend_block(
        direction    = 1,
        edge_marker  = seed_marker_right,
        marker_names = marker_names,
        marker_idx   = marker_idx,
        assigned     = assigned,
        ld_lookup    = ld_lookup,
        block        = block,
        method       = method,
        threshold    = threshold,
        tolerance    = tolerance,
        tol_reset    = tol_reset
      )

      assigned[block] = TRUE
      chrom_blocks[[length(chrom_blocks) + 1]] = block

      # remove all adjacent pairs that involve markers now in this block
      pairs_in_block = with(ld_adj, Name1 %in% block | Name2 %in% block)
      ld_adj = ld_adj[!pairs_in_block, ]

      if (nrow(ld_adj) == 0) break
    }

  } else if (start == "beginning") {

    position   = 1
    total_snps = length(marker_names)

    repeat {

      seed_marker = as.character(marker_names[position])
      block          = seed_marker

      block = extend_block(
        direction    = 1,
        edge_marker  = seed_marker,
        marker_names = marker_names,
        marker_idx   = marker_idx,
        assigned     = assigned,
        ld_lookup    = ld_lookup,
        block        = block,
        method       = method,
        threshold    = threshold,
        tolerance    = tolerance,
        tol_reset    = tol_reset
      )

      assigned[block] = TRUE
      chrom_blocks[[length(chrom_blocks) + 1]] = block

      position = marker_idx[[block[length(block)]]] + 1
      if (position > total_snps) break
    }
  }

  return(list(chrom_blocks, assigned))
}

# make_blocks_c ---------------------------------------------------------------
# Wrapper around the C++ implementation of make_blocks. Accepts the same
# arguments and returns the same result, but delegates to compiled code.
# This is much faster than the R implementation, especially on larger chromosomes, and is the default in chromo_blocking.
make_blocks_c = function(ld_lookup, ld_adj, marker_names, marker_idx,
                          marker_positions, first_marker, last_marker, assigned,
                          method, threshold, tolerance, tol_reset, start) {
  make_blocks_cpp(
    ld_lookup    = ld_lookup,
    ld_adj       = ld_adj,
    marker_names = marker_names,
    marker_idx   = marker_idx,
    assigned     = assigned,
    method       = method,
    threshold    = threshold,
    tolerance    = as.integer(tolerance),
    tol_reset    = tol_reset,
    start        = start
  )
}

# chromo_blocking --------------------------------------------------------------
# Coordinates haploblock formation for a single chromosome. Prepares the
# chromosome-level data, builds the O(1) lookup structures, calls make_blocks,
# then adds any unassigned markers as singleton blocks.
#
# chr       : chromosome identifier (matches values in ld$Chrom)
# ld        : pairwise LD table, pre-filtered to this chromosome
# map       : marker map table with columns SNP, Chromosome, Position
# method    : LD evaluation method ("flanking" or "average")
# tolerance : consecutive below-threshold markers allowed during block extension
# tol_reset : whether to reset the tolerance counter on each accepted marker
# threshold : minimum LD value to seed or extend a block
# start     : seed strategy ("LD" or "beginning")
chromo_blocking = function(chr, ld, map, method, tolerance, tol_reset,
                           threshold, start) {

  ld_chrom = ld[ld$Chrom == chr, ]
  ld_adj   = ld_chrom[ld_chrom$Locus2 == ld_chrom$Locus1 + 1, ]

  marker_names     = unique(c(ld_chrom$Name1, ld_chrom$Name2))
  marker_positions = map$Position[match(marker_names, map$SNP)]
  marker_names     = marker_names[order(marker_positions)]

  first_marker = marker_names[1]
  last_marker  = marker_names[length(marker_names)]

  assigned        = rep(FALSE, length(marker_names))
  names(assigned) = marker_names

  # build once per chromosome; passed down to avoid rebuilding on every extension step
  ld_lookup  = setNames(ld_chrom$LD, paste(ld_chrom$Name1, ld_chrom$Name2, sep = ","))
  marker_idx = setNames(seq_along(marker_names), marker_names)

  # call the main blocking function implemented in C++
  result = make_blocks_c(
    ld_lookup        = ld_lookup,
    ld_adj           = ld_adj,
    marker_names     = marker_names,
    marker_idx       = marker_idx,
    marker_positions = marker_positions,
    first_marker     = first_marker,
    last_marker      = last_marker,
    assigned         = assigned,
    method           = method,
    threshold        = threshold,
    tolerance        = tolerance,
    tol_reset        = tol_reset,
    start            = start
  )

  # Alternatively, to use the R implementation of make_blocks, comment out the above call to make_blocks_c and uncomment the line below to call make_blocks instead.
  # result = make_blocks(
  #   ld_lookup        = ld_lookup,
  #   ld_adj           = ld_adj,
  #   marker_names     = marker_names,
  #   marker_idx       = marker_idx,
  #   marker_positions = marker_positions,
  #   first_marker     = first_marker,
  #   last_marker      = last_marker,
  #   assigned         = assigned,
  #   method           = method,
  #   threshold        = threshold,
  #   tolerance        = tolerance,
  #   tol_reset        = tol_reset,
  #   start            = start
  # )

  chrom_blocks = result[[1]]
  assigned     = result[[2]]

  unassigned = names(assigned)[!assigned]
  for (snp in unassigned) {
    chrom_blocks[[length(chrom_blocks) + 1]] = snp
  }

  return(chrom_blocks)
}


# def_blocks -------------------------------------------------------------------
# Top-level function. Splits the LD data by chromosome, optionally parallelises
# across chromosomes, and returns a named list of blocks per chromosome.
#
# ld        : pairwise LD table across all chromosomes
#             (columns: Chrom, Locus1, Locus2, Name1, Name2, LD)
# map       : marker map table with columns SNP, Chromosome, Position
# method    : "flanking" — extend using LD between candidate and edge marker only;
#             "average"  — extend using mean LD between candidate and all markers in the block
# tolerance : number of consecutive below-threshold markers tolerated during extension
# tol_reset : if TRUE, reset tolerance counter each time a marker is accepted
# threshold : minimum LD (r²) required to seed or extend a block
# start     : "LD"        — seed blocks from highest-LD adjacent pairs first;
#             "beginning" — sweep chromosome left to right from the first marker
# parallel  : if TRUE, process chromosomes in parallel using all available cores minus one
def_blocks = function(ld, map, method = c("flanking", "average"), tolerance = 1,
                      tol_reset = TRUE, threshold = 0.7,
                      start = c("LD", "beginning"), parallel = FALSE) {

  start  = match.arg(start)
  method = match.arg(method)

  chromosomes = sort(unique(ld$Chrom))
  ld          = split(ld, ld$Chrom)

  if (parallel == TRUE) {
    future::plan(multisession, workers = parallel::detectCores() - 1)
    map_fun = furrr::future_map
  } else {
    map_fun = purrr::map
  }

  handlers("txtprogressbar")
  with_progress({
    p = progressor(steps = length(chromosomes))

    blocks_list = map_fun(ld, function(chromosome) {
      chr          = unique(chromosome$Chrom)
      chrom_blocks = chromo_blocking(
        chr       = chr,
        ld        = chromosome,
        map       = map,
        method    = method,
        tolerance = tolerance,
        threshold = threshold,
        start     = start,
        tol_reset = tol_reset
      )
      p()
      return(chrom_blocks)
    })
  })

  if (parallel == TRUE) plan(sequential)

  names(blocks_list) = as.character(chromosomes)
  return(blocks_list)
}


# chromo_blocks_to_df ----------------------------------------------------------
# Converts the list of blocks for one chromosome into a table with one row per
# block, including first/last marker names and physical positions.
#
# chrom_blocks : list of blocks, where each block is an ordered list of marker names
# map          : marker map table with columns SNP, Chromosome, Position
chromo_blocks_to_df = function(chrom_blocks, map) {

  block_df = map_dfr(chrom_blocks, function(block) {
    first_marker    = block[1]
    last_marker     = block[length(block)]
    block_length = length(block)
    marker_string   = paste(block, collapse = ";")
    data.frame(
      Block     = marker_string,
      Num_SNP   = block_length,
      First_SNP = first_marker,
      Last_SNP  = last_marker
    )
  })

  block_df = left_join(block_df, map, c("First_SNP" = "SNP"))

  map      = map[, c("SNP", "Position")]
  block_df = left_join(block_df, map, c("Last_SNP" = "SNP"))

  colnames(block_df) = c(colnames(block_df)[1:4], "Chrom", "Start_Pos", "End_Pos")

  block_df          = block_df[order(block_df$Start_Pos), ]
  block_df$Block_ID = 1:nrow(block_df)
  block_df$Block_ID = paste(block_df$Chrom, block_df$Block_ID, sep = ":")

  return(block_df)
}


# block_obj_to_df --------------------------------------------------------------
# Converts the full block list across all chromosomes returned by def_blocks
# into a single flat table with one row per block.
#
# block_obj : per-chromosome block lists as returned by def_blocks
# map       : marker map table with columns SNP, Chromosome, Position
block_obj_to_df = function(block_obj, map) {

  map      = map[, 1:3]
  block_df = map_dfr(block_obj, function(chrom_blocks) {
    chromo_blocks_to_df(chrom_blocks, map)
  })

  block_df = block_df[, c(1, 8, 2:7)]
  block_df$Block_Size = (block_df$End_Pos - block_df$Start_Pos)

  return(block_df)
}


# block_summary ----------------------------------------------------------------
# Computes summary statistics across all blocks in a haploblock table.
#
# block_df : haploblock table as returned by block_obj_to_df
block_summary = function(block_df) {

  mean_snp        = mean(block_df$Num_SNP, na.rm = TRUE)
  max_snp         = max(block_df$Num_SNP,  na.rm = TRUE)
  mean_size       = mean(block_df[block_df$Block_Size != 0, "Block_Size"], na.rm = TRUE)
  max_size        = max(block_df[block_df$Block_Size  != 0, "Block_Size"], na.rm = TRUE)
  singletons      = nrow(block_df[block_df$Num_SNP == 1 & !is.na(block_df$Num_SNP), ])
  perc_singletons = singletons / nrow(block_df[!is.na(block_df$Num_SNP), ])

  data.frame(
    Mean_SNP_per_Block       = mean_snp,
    Max_SNP_per_Block        = max_snp,
    Mean_Block_Size          = mean_size,
    Max_Block_Size           = max_size,
    Singleton_Blocks         = singletons,
    Percent_Singleton_Blocks = perc_singletons * 100
  )
}
