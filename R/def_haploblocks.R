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

# graph_strategy ---------------------------------------------------------------
# Builds a graph-based blocking strategy for use with def_blocks(). Blocks are
# connected components of a local LD graph, projected back onto the map and cut
# into non-overlapping segments.
#
# theta_core       : minimum r^2 for a core graph edge
# theta_core_by_chr: named numeric vector of per-chromosome theta_core overrides, named by
#                    chromosome as it appears in the map, e.g. c("1" = 0.9, "2" = 0.8);
#                    NULL to use theta_core everywhere
# theta_extend     : minimum r^2 for attaching an ungrouped marker to a block
# theta_bridge     : minimum r^2 for a boundary edge merging two adjacent blocks
# theta_refill     : minimum r^2 for refilling a marker left inside a block's span.
#                    Recomputed from the genotypes, so window_ld and ld_min_r2 do not cap it
# ld_min_r2        : floor for the shared LD edge table, below which edges are not stored
# window_ld        : forward marker window for LD calculation
# window_core      : maximum marker distance for a core edge
# window_extend    : maximum distance to the nearest member of a block being extended
# min_links        : minimum number of supporting edges for an extension
# max_gap_snps     : maximum intervening markers a bridge may cross. Bridging merges the two
#                    blocks; it does not absorb what lies between them
# max_gap_markers  : maximum marker gap within one segment, before it is split
# max_gap_position : maximum map distance within one segment, or NULL for no distance limit
# min_block_snps   : minimum markers for a candidate segment to be kept
# parallel         : if TRUE, process chromosomes in parallel using all available cores
#                    minus one, as ld_strategy() does
graph_strategy = function(theta_core        = 0.80,
                          theta_core_by_chr = NULL,
                          theta_extend      = 0.20,
                          theta_bridge      = 0.20,
                          theta_refill      = 0.80,
                          ld_min_r2         = 0.20,
                          window_ld         = 20,
                          window_core       = 20,
                          window_extend     = 50,
                          min_links         = 1,
                          max_gap_snps      = 2,
                          max_gap_markers   = 2,
                          max_gap_position  = NULL,
                          min_block_snps    = 2,
                          parallel          = FALSE) {

  # Check all thresholds are valid
  thresholds = list(theta_core = theta_core, theta_extend = theta_extend,
                    theta_bridge = theta_bridge, theta_refill = theta_refill,
                    ld_min_r2 = ld_min_r2)
  for (name in names(thresholds)) {
    value = thresholds[[name]]
    if (!is.numeric(value) || length(value) != 1 || is.na(value) || value < 0 || value > 1) {
      stop(name, " must be a single number between 0 and 1, since it is an r^2 threshold.")
    }
  }

  # Three separate ways to get the per-chromosome overrides wrong, so three separate checks
  if (!is.null(theta_core_by_chr)) {
    #Check the overrides are numeric
    if (!is.numeric(theta_core_by_chr) || length(theta_core_by_chr) == 0) {
      stop("theta_core_by_chr must be a numeric vector of r^2 thresholds, or NULL.")
    }

    # Check the overrides are named
    if (is.null(names(theta_core_by_chr)) || any(!nzchar(names(theta_core_by_chr)))) {
      stop("theta_core_by_chr must be named by chromosome, e.g. c(\"1\" = 0.9, ",
           "\"2\" = 0.8). Without names there is nothing to match against the map.")
    }

    # Check all overrides are between 0 and 1
    if (any(is.na(theta_core_by_chr) | theta_core_by_chr < 0 | theta_core_by_chr > 1)) {
      stop("Every theta_core_by_chr value must be between 0 and 1, since they are ",
           "r^2 thresholds.")
    }
  }

  # Check all marker indices or distances are valid, and non fractional
  check_property = function(value, name, minimum) {
    if (!is.numeric(value) || length(value) != 1 || is.na(value) ||
        value != round(value) || value < minimum) {
      stop(name, " must be a single whole number of at least ", minimum,
           ", since it counts markers.")
    }
  }
  check_property(window_ld, "window_ld", 1)
  check_property(window_core, "window_core", 1)
  check_property(window_extend, "window_extend", 1)
  check_property(min_links, "min_links", 1)
  check_property(min_block_snps, "min_block_snps", 1)
  check_property(max_gap_snps, "max_gap_snps", 0)
  check_property(max_gap_markers, "max_gap_markers", 0)

  # A map distance of zero or less would split a segment at every marker
  if (!is.null(max_gap_position)) {
    if (!is.numeric(max_gap_position) || length(max_gap_position) != 1 ||
        is.na(max_gap_position) || max_gap_position <= 0) {
      stop("max_gap_position must be a single positive map distance, or NULL.")
    }
  }

  # Core blocks are the connected components of the theta_core edges
  if (ld_min_r2 > theta_core) {
    stop("ld_min_r2 (", ld_min_r2, ") is above theta_core (", theta_core,
         "), so core blocks would be called at ld_min_r2 instead.")
  }

  # An override stands in for theta_core on the chromosomes it names
  if (!is.null(theta_core_by_chr) && any(ld_min_r2 > theta_core_by_chr)) {
    below = names(theta_core_by_chr)[ld_min_r2 > theta_core_by_chr]
    stop("ld_min_r2 (", ld_min_r2, ") is above theta_core_by_chr for chromosome ",
         paste(below, collapse = ", "), ", so core blocks there would be called at ",
         "ld_min_r2 instead.")
  }

  # Ungrouped markers join a block on the strength of the theta_extend edges
  if (ld_min_r2 > theta_extend) {
    stop("ld_min_r2 (", ld_min_r2, ") is above theta_extend (", theta_extend,
         "), so extension would run at ld_min_r2 instead.")
  }

  # Adjacent blocks merge across a gap on the strength of the theta_bridge edges
  if (ld_min_r2 > theta_bridge) {
    stop("ld_min_r2 (", ld_min_r2, ") is above theta_bridge (", theta_bridge,
         "), so bridging would run at ld_min_r2 instead.")
  }

  # Not wrong, just a no-op, so it should not stop a parameter sweep
  if (window_core > window_ld) {
    warning("window_core (", window_core, ") reaches past window_ld (", window_ld,
            "), and pairs beyond window_ld are never calculated. The effective core ",
            "window is window_ld.")
  }

  base_block_strategy(
    list(theta_core = theta_core, theta_core_by_chr = theta_core_by_chr,
         theta_extend = theta_extend, theta_bridge = theta_bridge,
         theta_refill = theta_refill, ld_min_r2 = ld_min_r2,
         window_ld = window_ld, window_core = window_core,
         window_extend = window_extend, min_links = min_links,
         max_gap_snps = max_gap_snps, max_gap_markers = max_gap_markers,
         max_gap_position = max_gap_position, min_block_snps = min_block_snps,
         parallel = parallel),
    class = "graph_strategy"
  )
}


##################################
#### Haploblocking Function ######
##################################

# def_blocks -------------------------------------------------------------------
# Top-level function. Performs blocking given a particular strategy.
#
# strategy : a strategy object built by ld_strategy(), window_strategy(),
#            graph_strategy(), or another *_strategy() constructor. A strategy
#            holds how to block, not what to block.
# map      : marker map table with columns SNP, Chromosome, Position
# geno     : genotype data frame in the HapSelect layout - SNP ID, chromosome and
#            position in columns 1 to 3, one dosage column per individual after
#            that. Required by graph_strategy(), which computes LD from it directly, and unused by
#            the other strategies.
def_blocks = function(strategy, map, geno = NULL) {

  if (!inherits(strategy, "block_strategy")) {
    stop("strategy must be built with a strategy constructor, e.g. ",
         "ld_strategy() or window_strategy().")
  }

  # Quietly accepting genotypes a strategy never reads would hide the mistake until the
  # blocks came back looking nothing like the data that was passed
  if (!is.null(geno) && !inherits(strategy, "graph_strategy")) {
    warning("geno is only used by graph_strategy(), and is ignored by ",
            class(strategy)[1], ".")
  }

  switch(class(strategy)[1],

    "ld_strategy" = perform_ld_blocking(
      ld        = strategy$ld,
      map       = map,
      method    = strategy$method,
      tolerance = strategy$tolerance,
      tol_reset = strategy$tol_reset,
      threshold = strategy$threshold,
      start     = strategy$start,
      parallel  = strategy$parallel
    ),

    "window_strategy" = perform_window_blocking(
      map    = map,
      window = strategy$window,
      method = strategy$method
    ),

    "graph_strategy" = perform_graph_blocking(
      geno     = geno,
      map      = map,
      strategy = strategy
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
