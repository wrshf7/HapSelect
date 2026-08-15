##########################################
###### Window-Based Haploblocks ##########
##########################################
#
# Fixed-size haploblocks cut straight from the map, with no reference to linkage
# disequilibrium. An alternative to the LD-driven blocking in def_haploblocks.R

# chromo_windows ---------------------------------------------------------------
# Partitions the markers of one chromosome into non-overlapping windows. Returns
# a list of blocks, each an ordered list of marker names, with the blocks
# themselves ordered by position.
#
# chrom_map : marker map for one chromosome (columns SNP, Position), sorted by Position
# window    : a marker count for "window_snp", or a distance for "window_map"
# method    : "window_snp" — consecutive runs of `window` markers. The last block
#                            on a chromosome is short when the marker count does
#                            not divide evenly.
#             "window_map" — half-open intervals [0, w), [w, 2w), ... anchored at
#                            position 0. Intervals holding no markers produce no
#                            block, so the number of blocks follows marker
#                            density rather than chromosome length.
chromo_windows = function(chrom_map, window, method) {

  marker_names = as.character(chrom_map$SNP)

  if (method == "window_snp") {
    # 0-based marker index integer-divided by the window size
    bin = (seq_along(marker_names) - 1) %/% window
  } else {
    # grid index of the half-open interval containing each position
    bin = floor(chrom_map$Position / window)
  }

  # Applies the bin index to the marker names, producing a list of blocks.
  # bin only needs to group and order markers correctly here; its actual
  # value is discarded by unname(), so it doesn't matter that it's 0-based.
  unname(split(marker_names, bin))
}


# perform_window_blocking ------------------------------------------------------------
# Splits the map by chromosome, cuts each chromosome into fixed windows, and
# returns a named list of blocks per chromosome. The return structure matches
# def_blocks(), so the result can be passed to block_obj_to_df(),
# compute_local_GEBV() and the plotting functions.
#
# map    : marker map as returned by order_map(), with columns SNP, Chromosome
#          and Position, sorted by Position within Chromosome. An unsorted map is
#          sorted here with a warning.
# window : window size. A whole number of markers for method = "window_snp", or a
#          distance in the map's own Position units (base pairs or cM) for
#          method = "window_map", may be fractional.
# method : "window_snp" — fixed number of markers per block
#          "window_map" — fixed distance per block
perform_window_blocking = function(map, window, method = c("window_snp", "window_map")) {

  method = match.arg(method)

  # Check the map structure and window argument.
  if (!is.data.frame(map) ||
      !all(c("SNP", "Chromosome", "Position") %in% colnames(map)) ||
      !is.numeric(map$Position)) {
    stop("Map must be a data frame with numeric positions and columns: ",
         "SNP, Chromosome, Position. Use order_map() to generate this object.")
  }

  # Check the window argument.
  if (missing(window) || !is.numeric(window) || length(window) != 1 ||
      is.na(window) || !is.finite(window) || window <= 0) {
    stop("Window must be a single positive number.")
  }

  # Check that the window is a whole number of markers when method = "window_snp".
  if (method == "window_snp" && window != as.integer(window)) {
    stop("Window must be a whole number of markers when method = 'window_snp'.")
  }

  # A marker with no position or no chromosome cannot be placed in any window
  incomplete = is.na(map$Position) | is.na(map$Chromosome)
  if (any(incomplete)) {
    warning(sum(incomplete), " marker(s) with a missing Position or Chromosome were ",
            "dropped and will not appear in any block.")
    map = map[!incomplete, ]
  }

  # Check that the map contains at least one marker with both a position and a chromosome.
  if (nrow(map) == 0) {
    stop("Map contains no markers with both a Position and a Chromosome.")
  }

  chromosomes = sort(unique(map$Chromosome))
  map_split   = split(map, map$Chromosome)

  # Blocks are cut from the map's row order, so markers must run in position
  # order within each chromosome.
  # Note: is.unsorted is a c level function so it is highly efficient.
  unsorted = vapply(map_split, function(cm) is.unsorted(cm$Position), logical(1))

  # Sort any unsorted chromosomes and issue a warning naming the first out-of-order marker.
  if (any(unsorted)) {
    cm = map_split[[which(unsorted)[1]]]
    first = which(diff(cm$Position) < 0)[1] + 1
    warning("Map is not sorted by Position within Chromosome (first out-of-order ",
            "marker: '", cm$SNP[first], "' on chromosome ", cm$Chromosome[first],
            "); sorting it. Run order_map() on your map to avoid this.")

    map_split = lapply(map_split, function(cm) cm[order(cm$Position), ])
  }

  # Split each chromosome into windows and return a named list of blocks per chromosome.
  blocks_list = purrr::map(map_split, function(chrom_map) {
    chromo_windows(chrom_map = chrom_map, window = window, method = method)
  })

  names(blocks_list) = as.character(chromosomes)
  return(blocks_list)
}
