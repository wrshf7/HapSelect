##########################################
####### Graph-Based Haploblocks ##########
##########################################
#
# Haploblocks as connected components of a local LD graph. Where the LD blocking
# in blocking_ld.R grows a block outward from a seed pair one marker at a time,
# this method builds a sparse LD graph over each chromosome and reads blocks off
# its structure, so a block is defined by how its markers hang together as a
# group rather than by a single chain of pairwise comparisons.
#
# The method runs in three stages per chromosome:
#
#   Stage I   graph construction. Compute LD within a forward marker window,
#             keep the edges above theta_core, and take connected components of
#             at least two markers as core blocks. Ungrouped markers are then
#             attached to a neighbouring block when they have enough LD support
#             (theta_extend), and adjacent blocks separated by a small marker
#             gap are merged when their boundary markers are in LD
#             (theta_bridge). Whatever is left over stays a one-marker block.
#
#   Stage II  linearisation. A graph component can be non-contiguous in marker
#             order, which a haploblock cannot be, so each block is projected
#             back onto the map, split where its members are too far apart, and
#             a non-overlapping set of the resulting segments is chosen.
#
#   Stage III refill. Linearisation leaves markers sitting inside a selected
#             block's span without belonging to it. Each is re-tested against
#             that block at theta_refill: it joins the block if it has strong LD
#             to any member, and is otherwise dropped. Markers outside every
#             selected span are kept as one-marker blocks.
#
# This is a port of the reference implementation kept in
# inst/examples/blocking_graph_prototype.R, following it closely enough that the
# two can be compared stage by stage, which is what
# tests/testthat/test-graph-blocking.R does. Where the prototype duplicated
# something the package already had - its own GT to dosage conversion, its own
# VCF reader, its own block summaries - the package's version is used instead.
# The one algorithm swapped out is the connected-component search, which used
# igraph and now uses connected_components_cpp().
#
# The strategy object that configures all of this is graph_strategy(), which
# lives in def_haploblocks.R beside the other blocking strategies.


# geno_matrix ------------------------------------------------------------------
# Turns the HapSelect genotype layout, which has markers as rows, into the
# sample-by-marker double matrix the C++ routines read, with marker names on the
# columns. Both of them index markers by column and samples by row.
#
# geno_chr: genotype rows for one chromosome, markers as rows
geno_matrix = function(geno_chr) {
  dosages = t(as.matrix(geno_chr[, -(1:3), drop = FALSE]))
  storage.mode(dosages) = "double"
  colnames(dosages) = as.character(geno_chr[[1]])
  rownames(dosages) = NULL
  dosages
}


# bidirectional_edges ----------------------------------------------------------
# Expands each undirected edge into both directions and keys the result on the
# marker being looked up, so that "every edge touching marker s" is one keyed
# lookup. Extension and bridging both need this; it adds no evidence, it is only
# a change of representation.
#
# edges : edge table from ld_func_c()
# min_r2: keep only edges at or above this first
bidirectional_edges = function(edges, min_r2 = NULL) {
  dt = as.data.table(edges)
  if (!is.null(min_r2)) dt = dt[LD >= min_r2]

  if (nrow(dt) == 0) {
    out = data.table(from = character(), to = character(), LD = numeric())
    setkey(out, from)
    return(out)
  }

  out = rbind(
    dt[, .(from = Name1, to = Name2, LD)],
    dt[, .(from = Name2, to = Name1, LD)]
  )
  setkey(out, from)
  out
}


# marker_index -----------------------------------------------------------------
# A marker's index is its row in the chromosome's map, which is what the edge
# table's Locus columns report and what every distance and gap in the method is
# measured in. This is the lookup from name to index used throughout.
#
# map_chr: map rows for one chromosome, ordered by position
marker_index = function(map_chr) {
  stats::setNames(seq_len(nrow(map_chr)), as.character(map_chr$SNP))
}


##################################
###### Stage I: LD graph #########
##################################

# core_blocks ------------------------------------------------------------------
# Takes the connected components of the LD graph as the starting blocks. Every
# marker on the chromosome is a vertex; an edge joins two markers when their r^2
# is at or above theta_core and they are within window_core of each other.
#
# A component of one marker is not a block, so it is returned as unassigned for
# extend_blocks() to place. Note that a component can be non-contiguous in marker
# order - markers 4, 5 and 9 can share a component without marker 6 joining it.
# Stage II is what turns that back into something linear.
#
# edges      : edge table from ld_func_c()
# map_chr    : map rows for the chromosome, ordered by position
# theta_core : minimum r^2 for an edge to count
# window_core: maximum marker distance for an edge to count
#
# Returns a list of:
#   blocks     - list of character vectors of marker names, each ordered by marker
#                index, the list itself ordered by each block's first marker
#   unassigned - marker names in no block, in map order
core_blocks = function(edges, map_chr, theta_core, window_core) {

  dt = as.data.table(edges)[LD >= theta_core]
  if (!is.null(window_core)) dt = dt[abs(Locus1 - Locus2) <= window_core]

  if (nrow(dt) == 0) {
    return(list(blocks = list(), unassigned = as.character(map_chr$SNP)))
  }

  snp2idx = marker_index(map_chr)

  # Every marker is a vertex, including the ones no surviving edge touches, so
  # that a marker in no component still gets a label and can be reported as
  # unassigned below.
  membership = connected_components_cpp(
    from = as.integer(dt$Locus1),
    to = as.integer(dt$Locus2),
    n_vertex = nrow(map_chr)
  )

  comp_list = split(as.character(map_chr$SNP), membership)
  comp_list = unname(comp_list[lengths(comp_list) >= 2])

  # A component is a set, not a run of markers, so it is put into marker order
  # here and the components themselves ordered by where each one starts
  if (length(comp_list) > 0) {
    comp_list = lapply(comp_list, function(x) x[order(snp2idx[x])])
    comp_list = comp_list[order(vapply(comp_list, function(x) min(snp2idx[x]), numeric(1)))]
  }

  assigned = if (length(comp_list) == 0) character() else unique(unlist(comp_list, use.names = FALSE))

  list(blocks = comp_list, unassigned = setdiff(as.character(map_chr$SNP), assigned))
}


# extend_blocks ----------------------------------------------------------------
# Attaches ungrouped markers to a neighbouring block when they have enough LD
# support for it. Each candidate marker is scored against every block it has a
# qualifying edge to, and joins the best one.
#
# Assignment is sequential in the order of unassigned, and a marker that joins a
# block immediately becomes a member for the markers considered after it, so an
# earlier marker can pull in a later one that had no edge to the original block.
#
# A block is eligible for a marker when:
#   - at least min_links edges of r^2 >= theta_extend run between them, and
#   - the nearest actual member of that block is at most window_extend marker
#     positions away. Distance is to the nearest member, not to the block's span,
#     so a block with a hole in it does not count as nearby across the hole.
#
# Eligible blocks are ranked by most supporting edges, then highest mean r^2 over
# those edges, then nearest member, and the first is taken.
#
# blocks       : blocks so far, as returned by core_blocks()
# unassigned   : marker names to try to place, in the order they are tried
# edges        : edge table from ld_func_c()
# map_chr      : map rows for the chromosome, ordered by position
# theta_extend : minimum r^2 for an edge to support an attachment
# window_extend: maximum distance to the nearest member of a candidate block
# min_links    : minimum number of supporting edges
#
# Returns a list of blocks and unassigned in the same shape as core_blocks(),
# with each block still ordered by marker index.
extend_blocks = function(blocks, unassigned, edges, map_chr,
                         theta_extend, window_extend, min_links) {

  if (length(blocks) == 0 || length(unassigned) == 0) {
    return(list(blocks = blocks, unassigned = unassigned))
  }

  snp2idx = marker_index(map_chr)
  dt2 = bidirectional_edges(edges, min_r2 = theta_extend)
  if (nrow(dt2) == 0) return(list(blocks = blocks, unassigned = unassigned))

  # Which block each marker currently belongs to, so that an edge can be turned
  # into a vote for a block. Markers that join a block below are recorded here as
  # they go, which is what makes the pass sequential.
  to_block = rep(NA_integer_, nrow(map_chr))
  names(to_block) = as.character(map_chr$SNP)
  for (bid in seq_along(blocks)) to_block[blocks[[bid]]] = bid

  still_unassigned = character()

  for (s in unassigned) {
    s_idx = snp2idx[[s]]
    edges_s = dt2[list(s), nomatch = 0]
    if (nrow(edges_s) == 0) {
      still_unassigned = c(still_unassigned, s)
      next
    }

    bid_vec = to_block[edges_s$to]
    keep = !is.na(bid_vec)
    if (!any(keep)) {
      still_unassigned = c(still_unassigned, s)
      next
    }

    score_dt = data.table(block_id = as.integer(bid_vec[keep]), LD = edges_s$LD[keep])
    score_dt = score_dt[, .(
      n_links = .N,
      max_r2 = max(LD),
      mean_r2 = mean(LD)
    ), by = block_id]

    # Distance is to the nearest actual member, not to the block's span, so a
    # block with a hole in it does not count as nearby across the hole
    score_dt[, min_dist := vapply(block_id, function(bid) {
      min(abs(snp2idx[blocks[[bid]]] - s_idx))
    }, numeric(1))]
    score_dt = score_dt[min_dist <= window_extend & n_links >= min_links]

    if (nrow(score_dt) == 0) {
      still_unassigned = c(still_unassigned, s)
      next
    }

    # max_r2 is carried for inspection but takes no part in the ranking
    setorder(score_dt, -n_links, -mean_r2, min_dist)
    best_block = score_dt$block_id[1]

    blocks[[best_block]] = unique(c(blocks[[best_block]], s))
    blocks[[best_block]] = blocks[[best_block]][order(snp2idx[blocks[[best_block]]])]
    to_block[s] = best_block
  }

  list(blocks = blocks, unassigned = unique(still_unassigned))
}


# bridge_blocks ----------------------------------------------------------------
# Merges adjacent blocks whose facing ends are in LD across a small gap. Blocks
# are walked in marker order, and a block is merged with the one after it when:
#   - at most max_gap_snps marker positions lie strictly between them, and
#   - at least one edge of r^2 >= theta_bridge joins the last three markers of
#     the left block to the first three of the right.
#
# Only the two blocks merge. Markers sitting in the gap are not absorbed, so a
# bridged block can span markers it does not contain - which is exactly what
# Stage II then has to resolve. Merging chains left to right, so a merged block
# is itself a candidate for merging with the block after it.
#
# blocks      : blocks so far
# edges       : edge table from ld_func_c()
# map_chr     : map rows for the chromosome, ordered by position
# theta_bridge: minimum r^2 for a boundary edge
# max_gap_snps: largest gap, in intervening markers, that can still be bridged
#
# Returns the blocks after merging, ordered by first marker, each ordered by
# marker index.
bridge_blocks = function(blocks, edges, map_chr, theta_bridge, max_gap_snps) {

  if (length(blocks) <= 1) return(blocks)

  snp2idx = marker_index(map_chr)
  dt2 = bidirectional_edges(edges, min_r2 = theta_bridge)
  if (nrow(dt2) == 0) return(blocks)

  blocks = blocks[order(vapply(blocks, function(b) min(snp2idx[b]), numeric(1)))]
  merged = list()
  i = 1L

  while (i <= length(blocks)) {
    cur = blocks[[i]]

    if (i == length(blocks)) {
      merged[[length(merged) + 1L]] = cur
      break
    }

    nxt = blocks[[i + 1L]]
    gap_n = min(snp2idx[nxt]) - max(snp2idx[cur]) - 1L
    if (gap_n > max_gap_snps) {
      merged[[length(merged) + 1L]] = cur
      i = i + 1L
      next
    }

    # Only the markers facing the gap can carry a bridge
    cur_tail = utils::tail(cur, min(3L, length(cur)))
    nxt_head = utils::head(nxt, min(3L, length(nxt)))
    support = rbindlist(lapply(cur_tail, function(s) {
      dt2[list(s), nomatch = 0][to %in% nxt_head]
    }), fill = TRUE)

    if (nrow(support) > 0) {
      # The merged block carries forward as the next candidate, so a run of
      # blocks can chain together left to right. Markers lying in the gap are
      # not absorbed, which is what leaves the block with a hole for Stage II.
      blocks[[i + 1L]] = unique(c(cur, nxt))
      blocks[[i + 1L]] = blocks[[i + 1L]][order(snp2idx[blocks[[i + 1L]]])]
      i = i + 1L
    } else {
      merged[[length(merged) + 1L]] = cur
      i = i + 1L
    }
  }

  merged
}


# graph_chromosome_blocks ------------------------------------------------------
# Runs all of Stage I for one chromosome: build the edge table, call core blocks
# from it, extend ungrouped markers into them, bridge across small gaps, and keep
# everything still ungrouped as a one-marker block.
#
# geno_chr: genotype rows for one chromosome, ordered to match map_chr
# map_chr : map rows for that chromosome, ordered by position
# strategy: graph_strategy() object. theta_core is taken from it directly - the
#           per-chromosome override in theta_core_by_chr is resolved by
#           perform_graph_blocking() before this is called.
#
# Returns a list of character vectors of marker names, ordered by first marker.
# Every marker on the chromosome appears in exactly one of them, since the
# leftovers are carried as one-marker blocks rather than discarded.
graph_chromosome_blocks = function(geno_chr, map_chr, strategy) {

  # min_obs = 3 refuses an r^2 of 1 drawn from two shared individuals, which would
  # otherwise arrive as an edge the graph has no business trusting
  edges = ld_func_c(geno_chr, window = strategy$window_ld,
                    min_r2 = strategy$ld_min_r2, min_obs = 3L)

  core = core_blocks(edges, map_chr, strategy$theta_core, strategy$window_core)

  ext = extend_blocks(core$blocks, core$unassigned, edges, map_chr,
                      strategy$theta_extend, strategy$window_extend,
                      strategy$min_links)

  bridged = bridge_blocks(ext$blocks, edges, map_chr, strategy$theta_bridge,
                          strategy$max_gap_snps)

  # Whatever is still ungrouped is carried as a one-marker block, so that Stage I
  # accounts for every marker on the chromosome
  blocks = c(bridged, as.list(ext$unassigned))

  snp2idx = marker_index(map_chr)
  blocks = lapply(blocks, function(b) b[order(snp2idx[b])])
  blocks[order(vapply(blocks, function(b) min(snp2idx[b]), numeric(1)))]
}


##################################
#### Stage II: linearisation #####
##################################

# linearise_blocks -------------------------------------------------------------
# Projects graph blocks back onto the map and cuts them into candidate segments
# that are compact in marker order. A component that jumps across the chromosome
# is not a haploblock; this is what decides where it breaks.
#
# Members are put in marker order and a new segment starts wherever successive
# members are more than max_gap_markers positions apart, or, when
# max_gap_position is given, more than that far apart in map units. Segments with
# fewer than min_block_snps members are dropped. Their markers become uncovered,
# which Stage III then handles.
#
# A segment still describes the interval between its first and last marker, and
# may not contain every marker in it - that is what select_nonoverlapping() and
# the refill step deal with.
#
# blocks          : Stage I blocks for one chromosome
# map_chr         : map rows for that chromosome, ordered by position
# max_gap_markers : largest gap in marker positions within a segment
# max_gap_position: largest gap in map units within a segment, or NULL for no limit
# min_block_snps  : segments with fewer members than this are dropped
#
# Returns a data frame with one row per surviving segment, ordered by
# first_index then last_index, with columns:
#   chromosome, source_block (position of the originating block in blocks),
#   first_index, last_index    - marker indices bounding the segment
#   first_position, last_position, length - map positions and their difference
#   n_markers                  - members in the segment
#   span                       - last_index - first_index + 1
#   density                    - n_markers / span, so 1 when the segment is solid
#   score                      - selection score, see select_nonoverlapping()
#   markers                    - member names, ";"-separated, in marker order
#
# score is the heuristic the reference implementation uses to rank candidates:
#   1000 * n_markers + 100 * density - log1p(length)
# so marker count dominates, density breaks ties between equal counts, and a
# shorter segment is preferred when both match.
linearise_blocks = function(blocks, map_chr, max_gap_markers, max_gap_position,
                            min_block_snps) {

  empty = data.table(
    chromosome = character(), source_block = integer(),
    first_index = integer(), last_index = integer(),
    first_position = numeric(), last_position = numeric(), length = numeric(),
    n_markers = integer(), span = integer(), density = numeric(),
    score = numeric(), markers = character()
  )
  if (length(blocks) == 0) return(empty)

  chromosome = map_chr$Chromosome[1]
  map_keyed = as.data.table(map_chr)
  map_keyed[, idx := .I]
  setkey(map_keyed, SNP)

  segments = rbindlist(lapply(seq_along(blocks), function(b) {
    snps = blocks[[b]]
    snps = snps[!is.na(snps) & nzchar(snps)]
    if (length(snps) == 0) return(NULL)

    tmp = map_keyed[list(snps), .(SNP, Position, idx)]
    tmp = tmp[!is.na(idx)]
    tmp = unique(tmp, by = "SNP")
    if (nrow(tmp) == 0) return(NULL)
    setorder(tmp, idx)

    # A new segment starts at the first member, and wherever the step from the
    # previous member is too long by either measure
    gap_idx = c(0L, diff(tmp$idx))
    gap_pos = c(0, diff(tmp$Position))
    breaks = rep(FALSE, nrow(tmp))
    breaks[1] = TRUE
    if (!is.null(max_gap_markers))  breaks = breaks | (gap_idx > max_gap_markers)
    if (!is.null(max_gap_position)) breaks = breaks | (gap_pos > max_gap_position)

    segs = tmp[, .(
      chromosome = chromosome,
      source_block = b,
      first_index = min(idx),
      last_index = max(idx),
      first_position = min(Position),
      last_position = max(Position),
      n_markers = .N,
      markers = paste(SNP, collapse = ";")
    ), by = .(seg_id = cumsum(breaks))]

    segs[, seg_id := NULL]
    segs[]
  }), use.names = TRUE, fill = TRUE)

  if (nrow(segments) == 0) return(empty)

  segments[, length := last_position - first_position]
  segments[, span := last_index - first_index + 1L]
  segments[, density := n_markers / pmax(span, 1L)]
  segments[, score := n_markers * 1000 + density * 100 - log1p(length)]

  segments = segments[n_markers >= min_block_snps]
  if (nrow(segments) == 0) return(empty)

  setcolorder(segments, names(empty))
  setorder(segments, first_index, last_index)
  segments[]
}


# select_nonoverlapping --------------------------------------------------------
# Picks a set of segments that do not overlap in marker index. Candidates are
# taken in descending score, then descending marker count, then ascending map
# length, then ascending first_index, and a candidate is kept when its whole
# [first_index, last_index] interval is still free.
#
# The interval is claimed entire, holes included, so a segment that skips markers
# still blocks anything that would have covered them. This is a greedy pass, not
# a search for the best overall set.
#
# segments: candidate segments from linearise_blocks()
#
# Returns the kept rows, same columns, ordered by first_index then last_index.
select_nonoverlapping = function(segments) {

  candidates = as.data.table(segments)
  if (nrow(candidates) == 0) return(candidates)

  setorder(candidates, -score, -n_markers, length, first_index)
  occupied = data.table(first_index = integer(), last_index = integer())
  selected = vector("list", nrow(candidates))
  keep_n = 0L

  overlaps_any = function(a1, a2, occ) {
    if (nrow(occ) == 0) return(FALSE)
    any(!(a2 < occ$first_index | a1 > occ$last_index))
  }

  for (i in seq_len(nrow(candidates))) {
    a1 = candidates$first_index[i]
    a2 = candidates$last_index[i]
    if (!overlaps_any(a1, a2, occupied)) {
      keep_n = keep_n + 1L
      selected[[keep_n]] = candidates[i]
      # The whole interval is claimed, holes included
      occupied = rbind(occupied, data.table(first_index = a1, last_index = a2))
    }
  }

  out = rbindlist(selected[seq_len(keep_n)], use.names = TRUE, fill = TRUE)
  setorder(out, first_index, last_index)
  out
}


##################################
###### Stage III: refill #########
##################################

# block_ld_support -------------------------------------------------------------
# Tests whether a marker is in strong enough LD with a block to belong to it.
# TRUE as soon as any one member reaches the threshold: this is a test of the
# best link, not of mean LD across the block, and not of how many members agree.
#
# LD is recomputed from the genotypes rather than read from the Stage I edge
# table, so neither window_ld nor ld_min_r2 constrains it and the marker is
# compared against every member however far away it sits.
#
# geno_chr  : genotype rows for the chromosome
# target_snp: name of the marker being tested
# block_snps: names of the block's current members
# threshold : minimum r^2 against any one member
#
# Returns TRUE or FALSE. FALSE when the marker is not in geno_chr, and when the
# block has no members other than the marker itself. r^2 is computed as in
# ld_func_c(min_obs = 3), so a pair with fewer than three shared observations or
# no variance does not count as support.
block_ld_support = function(geno_chr, target_snp, block_snps, threshold) {

  if (is.null(geno_chr) || nrow(geno_chr) == 0) return(FALSE)

  dosages = geno_matrix(geno_chr)
  if (!target_snp %in% colnames(dosages)) return(FALSE)

  members = intersect(block_snps, colnames(dosages))
  members = setdiff(members, target_snp)
  if (length(members) == 0) return(FALSE)

  member_cols = match(members, colnames(dosages))
  member_cols = member_cols[!is.na(member_cols)]
  if (length(member_cols) == 0) return(FALSE)

  has_strong_ld_to_block_cpp(
    geno = dosages,
    target_col = as.integer(match(target_snp, colnames(dosages))),
    member_cols = as.integer(member_cols),
    threshold = as.numeric(threshold)
  )
}


# refill_internal_markers ------------------------------------------------------
# Decides what happens to markers that no selected segment contains. A marker is
# internal when its index falls inside some selected segment's span, and external
# when it falls outside every one of them.
#
# An internal marker is tested against the block whose span encloses it, and only
# that block. It is added if block_ld_support() holds at theta_refill, and
# dropped otherwise - dropped meaning it appears in no block at all, so the
# blocking does not account for every marker on the chromosome. External markers
# are returned for the caller to keep as one-marker blocks.
#
# selected    : selected segments from select_nonoverlapping()
# map_chr     : map rows for the chromosome, ordered by position
# geno_chr    : genotype rows for the chromosome
# theta_refill: minimum r^2 against any member for an internal marker to be refilled
#
# Returns a list of:
#   blocks   - selected, with refilled markers added to markers and n_markers,
#              first_position, last_position, length, density recomputed. The index
#              bounds cannot change, since a refilled marker was inside the span.
#   refilled - data frame of SNP, block (its row in blocks), and index
#   dropped  - data frame of the same shape, for internal markers that failed
#   external - names of uncovered markers outside every selected span, in map order
refill_internal_markers = function(selected, map_chr, geno_chr, theta_refill) {

  blocks = as.data.table(selected)
  snp2idx = marker_index(map_chr)
  no_markers = data.frame(SNP = character(), block = integer(), index = integer(),
                          stringsAsFactors = FALSE)

  covered = unlist(strsplit(paste(blocks$markers, collapse = ";"), ";", fixed = TRUE),
                   use.names = FALSE)
  covered = covered[nzchar(covered)]
  uncovered = setdiff(as.character(map_chr$SNP), covered)

  if (length(uncovered) == 0) {
    return(list(blocks = blocks, refilled = no_markers, dropped = no_markers,
                external = character()))
  }

  refilled = list()
  dropped = list()
  external = character()

  for (s in uncovered) {
    s_idx = snp2idx[[s]]

    # Internal means inside some selected block's span. A marker is only ever
    # tested against the block enclosing it, never a neighbour, and the first
    # enclosing block wins - the spans do not overlap, so there is normally
    # only one.
    inside = which(blocks$first_index <= s_idx & blocks$last_index >= s_idx)
    if (length(inside) == 0) {
      external = c(external, s)
      next
    }

    row_i = inside[1]
    members = unlist(strsplit(blocks$markers[row_i], ";", fixed = TRUE), use.names = FALSE)
    decision = data.frame(SNP = s, block = row_i, index = s_idx,
                          stringsAsFactors = FALSE)

    if (!block_ld_support(geno_chr, s, members, theta_refill)) {
      # Dropped, meaning it ends up in no block at all
      dropped[[length(dropped) + 1L]] = decision
      next
    }

    refilled[[length(refilled) + 1L]] = decision

    new_members = c(members, s)
    new_members = new_members[order(snp2idx[new_members])]
    positions = map_chr$Position[match(new_members, as.character(map_chr$SNP))]

    # The index bounds cannot move, since the marker was inside the span already
    set(blocks, row_i, "markers", paste(new_members, collapse = ";"))
    set(blocks, row_i, "n_markers", length(new_members))
    set(blocks, row_i, "first_position", min(positions))
    set(blocks, row_i, "last_position", max(positions))
    set(blocks, row_i, "length", max(positions) - min(positions))
    set(blocks, row_i, "density", length(new_members) / max(blocks$span[row_i], 1L))
  }

  list(
    blocks = blocks,
    refilled = if (length(refilled) == 0) no_markers else do.call(rbind, refilled),
    dropped = if (length(dropped) == 0) no_markers else do.call(rbind, dropped),
    external = external
  )
}


##################################
######## Orchestration ###########
##################################

# perform_graph_blocking -------------------------------------------------------
# Top-level graph blocking. Splits the genotypes and map by chromosome, runs the
# three stages over each, and returns blocks in the same shape as
# perform_ld_blocking() and perform_window_blocking(), so the result can be
# handed to block_obj_to_df(), compute_local_GEBV() and the plotting functions.
#
# geno    : genotype data frame in the HapSelect layout, as read_vcf_geno() produces
# map     : marker map with SNP, Chromosome and Position, ordered by order_map()
# strategy: graph_strategy() object
#
# geno and map arrive as separate arguments, so this is where they have to be
# checked against each other. Required:
#   - geno must be supplied, and must satisfy check_ld_matrix()
#   - geno and map must hold the same markers in the same order
#   - every name in strategy$theta_core_by_chr must match a chromosome in map
#
# The second of those is the one that matters. A marker's index is its row in the
# map, and dosages are read by row from geno, so if the two disagree every r^2 is
# computed between the wrong pair of markers. Nothing downstream notices: the run
# finishes and returns a full set of blocks built from nonsense. order_map()
# reorders rows, so handing this an ordered map with an unordered geno is an easy
# mistake to make and an invisible one to debug.
#
# The third is a silent no-op rather than a wrong answer - an override naming a
# chromosome that is not there simply never fires - but it fails the same way,
# by looking like it worked.
#
# Returns a named list with one entry per chromosome, each a list of blocks, each
# block a character vector of marker names in map order. Blocks are ordered by
# their first marker. Both the multi-marker blocks that survived Stage III and
# the one-marker blocks for uncovered markers are included.
#
# Markers dropped by the refill step are in no block and so do not appear at all,
# which means the output is not guaranteed to account for every marker in the
# map. That is worth knowing when comparing against the other blocking methods,
# which are strict partitions. The counts are reported in the diagnostics.
#
# A "graph_diagnostics" attribute carries one row per chromosome, with columns:
#   chromosome, theta_core (the value used, after any theta_core_by_chr override),
#   n_blocks, n_multi_marker_blocks, n_single_marker_blocks,
#   n_markers_in, n_markers_out, n_refilled, n_dropped
perform_graph_blocking = function(geno, map, strategy) {

  if (is.null(geno)) {
    stop("geno must be supplied for a graph_strategy: the method computes LD from ",
         "the genotypes directly. Pass it to def_blocks(strategy, map, geno = ...).")
  }
  check_ld_matrix(geno)

  geno_snps = as.character(geno[[1]])
  map_snps = as.character(map$SNP)
  if (length(geno_snps) != length(map_snps) || !setequal(geno_snps, map_snps)) {
    stop("geno and map must hold the same markers.")
  }
  if (!identical(geno_snps, map_snps)) {
    stop("geno and map must hold the same markers in the same order, since a ",
         "marker's index is its row in the map and its dosages are read from the ",
         "matching row of geno. order_map() reorders the map, so reorder geno to ",
         "match it.")
  }

  chromosomes = unique(map$Chromosome)

  if (!is.null(strategy$theta_core_by_chr)) {
    unknown = setdiff(names(strategy$theta_core_by_chr), as.character(chromosomes))
    if (length(unknown) > 0) {
      stop("theta_core_by_chr names no chromosome in the map: ",
           paste(unknown, collapse = ", "), ". An override that matches nothing ",
           "never fires, and the global theta_core would be used instead.")
    }
  }

  blocks_by_chr = list()
  diagnostics = list()

  for (chromosome in chromosomes) {
    key = as.character(chromosome)
    rows = map$Chromosome == chromosome
    map_chr = map[rows, , drop = FALSE]
    geno_chr = geno[rows, , drop = FALSE]
    snp2idx = marker_index(map_chr)

    # A per-chromosome override stands in for theta_core on this chromosome only
    theta_core = strategy$theta_core
    if (!is.null(strategy$theta_core_by_chr) && key %in% names(strategy$theta_core_by_chr)) {
      theta_core = unname(strategy$theta_core_by_chr[[key]])
    }
    chr_strategy = strategy
    chr_strategy$theta_core = theta_core

    stage1 = graph_chromosome_blocks(geno_chr, map_chr, chr_strategy)

    segments = linearise_blocks(stage1, map_chr, strategy$max_gap_markers,
                                strategy$max_gap_position, strategy$min_block_snps)
    selected = select_nonoverlapping(segments)

    if (nrow(selected) == 0) {
      # Nothing survived linearisation, so every marker stands alone
      blocks = as.list(as.character(map_chr$SNP))
      n_refilled = 0L
      n_dropped = 0L
    } else {
      refill = refill_internal_markers(selected, map_chr, geno_chr, strategy$theta_refill)

      blocks = lapply(refill$blocks$markers, function(m) {
        unlist(strsplit(m, ";", fixed = TRUE), use.names = FALSE)
      })
      blocks = c(blocks, as.list(refill$external))
      blocks = lapply(blocks, function(b) b[order(snp2idx[b])])
      blocks = blocks[order(vapply(blocks, function(b) min(snp2idx[b]), numeric(1)))]

      n_refilled = nrow(refill$refilled)
      n_dropped = nrow(refill$dropped)
    }

    blocks_by_chr[[key]] = blocks

    sizes = lengths(blocks)
    diagnostics[[key]] = data.frame(
      chromosome = chromosome,
      theta_core = theta_core,
      n_blocks = length(blocks),
      n_multi_marker_blocks = sum(sizes >= 2),
      n_single_marker_blocks = sum(sizes == 1),
      n_markers_in = nrow(map_chr),
      n_markers_out = sum(sizes),
      n_refilled = n_refilled,
      n_dropped = n_dropped,
      stringsAsFactors = FALSE
    )
  }

  attr(blocks_by_chr, "graph_diagnostics") =
    do.call(rbind, c(diagnostics, list(make.row.names = FALSE)))

  blocks_by_chr
}
