##################################
##### Haploblocking Methods ######
##################################

# base_block_strategy ------------------------------------------------------
# Internal base class shared by all strategy objects. Not created directly —
# use ld_strategy(), window_strategy(), or another *_strategy() constructor.
base_block_strategy = function(params, class) {
  structure(params, class = c(class, "block_strategy"))
}


# ld_strategy --------------------------------------------------------------
# Builds an LD-based blocking strategy for use with def_blocks(). Blocks are
# grown outward from seed marker pairs based on pairwise LD.
#
# ld        : pairwise LD table (columns: Chrom, Locus1, Locus2, Name1, Name2, LD)
# method    : "flanking" — extend using LD between candidate and edge marker only;
#             "average"  — extend using mean LD between candidate and all markers in the block
# tolerance : number of consecutive below-threshold markers tolerated during extension
# tol_reset : if TRUE, reset the tolerance counter each time a marker is accepted
# threshold : minimum LD (r^2) required to seed or extend a block
# start     : "LD"        — seed blocks from highest-LD adjacent pairs first;
#             "beginning" — sweep chromosome left to right from the first marker
# parallel  : if TRUE, process chromosomes in parallel using all available cores minus one
ld_strategy = function(ld, method = c("flanking", "average"), tolerance = 1,
                        tol_reset = TRUE, threshold = 0.7,
                        start = c("LD", "beginning"), parallel = FALSE) {

  method = match.arg(method)
  start  = match.arg(start)

  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("threshold must be a number between 0 and 1.")
  }
  if (!is.numeric(tolerance) || tolerance < 0) {
    stop("tolerance must be a non-negative number.")
  }

  base_block_strategy(
    list(ld = ld, method = method, tolerance = tolerance, tol_reset = tol_reset,
         threshold = threshold, start = start, parallel = parallel),
    class = "ld_strategy"
  )
}

# window_strategy ------------------------------------------------------------
# Builds a fixed-size window blocking strategy for use with def_blocks().
# Blocks are cut straight from the map, with no reference to linkage
# disequilibrium.
#
# window : a marker count for "window_snp", or a distance for "window_map"
# method : "window_snp" — fixed number of markers per block
#          "window_map" — fixed distance per block
window_strategy = function(window, method = c("window_snp", "window_map")) {
  method = match.arg(method)

  base_block_strategy(
    list(window = window, method = method),
    class = "window_strategy"
  )
}


##################################
#### Haploblocking Function ######
##################################

# def_blocks -------------------------------------------------------------------
# Top-level function. Performs blocking given a particular strategy.
#
# map      : marker map table with columns SNP, Chromosome, Position
# strategy : a strategy object built by ld_strategy(), window_strategy(), or
#            another *_strategy() constructor. Any data a strategy needs
#            besides the map (e.g. an LD table) is supplied when building the
#            strategy itself.
def_blocks = function(map, strategy) {

  if (!inherits(strategy, "block_strategy")) {
    stop("strategy must be built with a strategy constructor, e.g. ",
         "ld_strategy() or window_strategy().")
  }

  switch(class(strategy)[1],

    ld_strategy = perform_ld_blocking(
      ld        = strategy$ld,
      map       = map,
      method    = strategy$method,
      tolerance = strategy$tolerance,
      tol_reset = strategy$tol_reset,
      threshold = strategy$threshold,
      start     = strategy$start,
      parallel  = strategy$parallel
    ),

    window_strategy = perform_window_blocking(
      map    = map,
      window = strategy$window,
      method = strategy$method
    ),

    stop("No blocking method defined for strategy class '", class(strategy)[1], "'.")
  )
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
