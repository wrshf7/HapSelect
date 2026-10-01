# Tests: graph-based haploblocking ---------------------------------------------
#
# Covers R/blocking_graph.R. Every expectation here was taken from the original
# implementation this one was ported from, run over these same fixtures, so a
# failure means the behaviour has moved rather than that the expectation was
# guessed. Where that behaviour is surprising, the comment says so - the point of
# these tests is to notice when it changes, not to argue that it is right.
#
# Table comparisons go through as.data.frame() so the implementation is free to
# use data.table, dplyr or base R internally.


# Fixtures ---------------------------------------------------------------------
#
# Seven marker patterns over twelve individuals, chosen for the r2 they produce:
#
#          base  flip1 flip2  b3    b2    n1    n2
#   base   1.000 0.897 0.505 0.250 0.016 0.016 0.062
#   flip1  0.897 1.000 0.631 0.350 0.014 0.000 0.014
#   flip2  0.505 0.631 1.000 0.126 0.014 0.000 0.014
#   b3     0.250 0.350 0.126 1.000 0.016 0.141 0.016
#   b2     0.016 0.014 0.014 0.016 1.000 0.016 0.062
#   n1     0.016 0.000 0.000 0.141 0.016 1.000 0.141
#   n2     0.062 0.014 0.014 0.016 0.062 0.141 1.000
#
# base/flip1 clear a 0.8 core threshold, flip2 sits in between so it is only
# reachable by the refill step, and b2/n1/n2 are noise that nothing links to.

gb_base  <- c(0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 1, 2)
gb_flip1 <- c(0, 0, 0, 0, 1, 1, 2, 2, 2, 0, 1, 2)
gb_flip2 <- c(0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 1, 2)
gb_b3    <- c(0, 0, 1, 0, 2, 2, 1, 2, 1, 0, 2, 1)
gb_b2    <- c(1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2)
gb_n1    <- c(2, 1, 0, 2, 0, 1, 0, 2, 1, 1, 0, 2)
gb_n2    <- c(1, 2, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0)

# Builds a geno data frame from a list of list(name, chromosome, position, dosages)
gb_geno <- function(markers) {
  geno <- data.frame(
    SNP        = vapply(markers, function(m) m[[1]], character(1)),
    Chromosome = vapply(markers, function(m) m[[2]], numeric(1)),
    Position   = vapply(markers, function(m) m[[3]], numeric(1)),
    stringsAsFactors = FALSE
  )
  dosages <- do.call(rbind, lapply(markers, function(m) m[[4]]))
  colnames(dosages) <- paste0("Ind", seq_len(ncol(dosages)))
  cbind(geno, as.data.frame(dosages))
}

# The end-to-end fixture. Chromosome 1 exercises core calling, extension beyond
# window_core, bridging across a gap, and a refill that both keeps and drops a
# marker. Chromosome 2 exercises a bridge that is refused for want of a boundary
# edge, a refill that keeps a marker, and uncovered markers outside every block.
graph_geno_fixture <- function() {
  gb_geno(list(
    list("c1_01", 1, 100, gb_base),
    list("c1_02", 1, 200, gb_flip1),
    list("c1_03", 1, 300, gb_n1),
    list("c1_04", 1, 400, gb_flip2),
    list("c1_05", 1, 500, gb_flip1),
    list("c1_06", 1, 600, gb_n2),
    list("c1_07", 1, 700, gb_b3),
    list("c1_08", 1, 800, gb_b3),
    list("c2_01", 2, 100, gb_base),
    list("c2_02", 2, 200, gb_flip2),
    list("c2_03", 2, 300, gb_flip1),
    list("c2_04", 2, 400, gb_n1),
    list("c2_05", 2, 500, gb_b2),
    list("c2_06", 2, 600, gb_b2),
    list("c2_07", 2, 700, gb_n2)
  ))
}

graph_map_fixture <- function() graph_geno_fixture()[, 1:3]

# One chromosome of the fixture, as the per-chromosome functions take it
graph_chr_fixture <- function(chr, what = c("geno", "map")) {
  what <- match.arg(what)
  geno <- graph_geno_fixture()
  rows <- geno$Chromosome == chr
  if (what == "map") geno[rows, 1:3, drop = FALSE] else geno[rows, , drop = FALSE]
}

# window_core below window_ld is what leaves c1_05 to the extension step, and
# theta_refill below theta_extend is what leaves markers to the refill step
# rather than absorbing them in Stage I.
graph_test_strategy <- function(...) {
  args <- list(
    theta_core = 0.80, theta_extend = 0.70, theta_bridge = 0.20, theta_refill = 0.50,
    ld_min_r2 = 0.20, window_ld = 3, window_core = 2, window_extend = 3,
    min_links = 1, max_gap_snps = 2, max_gap_markers = 3, min_block_snps = 2
  )
  do.call(graph_strategy, utils::modifyList(args, list(...)))
}

# A plain map of n markers 100 apart on one chromosome
gb_map <- function(n, chrom = 1) {
  data.frame(SNP = sprintf("m%02d", seq_len(n)), Chromosome = chrom,
             Position = seq_len(n) * 100, stringsAsFactors = FALSE)
}

# An edge table from list(locus1, locus2, ld) triples
gb_edges <- function(map, pairs) {
  data.frame(
    Chrom  = map$Chromosome[1],
    Locus1 = vapply(pairs, function(p) p[1], numeric(1)),
    Locus2 = vapply(pairs, function(p) p[2], numeric(1)),
    Name1  = map$SNP[vapply(pairs, function(p) p[1], numeric(1))],
    Name2  = map$SNP[vapply(pairs, function(p) p[2], numeric(1))],
    LD     = vapply(pairs, function(p) p[3], numeric(1)),
    stringsAsFactors = FALSE
  )
}

# A geno data frame straight from a list of dosage vectors
gb_geno_rows <- function(rows, chrom = 1) {
  gb_geno(lapply(seq_along(rows), function(i) {
    list(sprintf("m%02d", i), chrom, i * 100, rows[[i]])
  }))
}


# graph_strategy ---------------------------------------------------------------

test_that("graph_strategy stores its arguments and carries the strategy classes", {
  strategy <- graph_test_strategy()

  expect_s3_class(strategy, "graph_strategy")
  expect_s3_class(strategy, "block_strategy")
  expect_equal(strategy$theta_core, 0.80)
  expect_equal(strategy$window_ld, 3)
  expect_equal(strategy$min_block_snps, 2)

  # a strategy holds how to block, not what to block: the genotypes reach
  # perform_graph_blocking() through def_blocks(), not through here
  expect_null(strategy$geno)
})


test_that("graph_strategy defaults match the reference implementation's parameters", {
  strategy <- graph_strategy()

  expect_equal(strategy$theta_core, 0.80)
  expect_equal(strategy$theta_extend, 0.20)
  expect_equal(strategy$theta_bridge, 0.20)
  expect_equal(strategy$theta_refill, 0.80)
  expect_equal(strategy$ld_min_r2, 0.20)
  expect_equal(strategy$window_ld, 20)
  expect_equal(strategy$window_core, 20)
  expect_equal(strategy$window_extend, 50)
  expect_equal(strategy$min_links, 1)
  expect_equal(strategy$max_gap_snps, 2)
  expect_equal(strategy$max_gap_markers, 2)
  expect_null(strategy$max_gap_position)
  expect_equal(strategy$min_block_snps, 2)
  expect_null(strategy$theta_core_by_chr)
})


test_that("graph_strategy rejects thresholds outside 0 to 1", {
  for (bad in list(list(theta_core = 1.5), list(theta_extend = -0.1),
                   list(theta_bridge = 2), list(theta_refill = -1),
                   list(ld_min_r2 = 1.2))) {
    expect_error(do.call(graph_test_strategy, bad), "between 0 and 1")
  }
})


test_that("graph_strategy rejects per-chromosome overrides it could not apply", {
  # Three distinct mistakes, so three distinct messages
  expect_error(graph_test_strategy(theta_core_by_chr = "0.9"), "numeric vector")
  expect_error(graph_test_strategy(theta_core_by_chr = numeric()), "numeric vector")

  # An unnamed vector has nothing to match chromosomes against, so it would be
  # ignored in full and the run would quietly use the global theta_core
  expect_error(graph_test_strategy(theta_core_by_chr = c(0.9, 0.8)), "named by chromosome")
  expect_error(graph_test_strategy(theta_core_by_chr = c("1" = 0.9, 0.8)),
               "named by chromosome")

  expect_error(graph_test_strategy(theta_core_by_chr = c("1" = 1.5)), "between 0 and 1")
  expect_error(graph_test_strategy(theta_core_by_chr = c("1" = NA_real_)),
               "between 0 and 1")

  expect_s3_class(graph_test_strategy(theta_core_by_chr = c("1" = 0.9, "2" = 0.8)),
                  "graph_strategy")
})


test_that("graph_strategy rejects an ld_min_r2 that would raise another threshold", {
  # Edges below ld_min_r2 are never computed, so a threshold underneath it silently
  # becomes ld_min_r2 - a wrong answer rather than a slow one
  expect_error(graph_test_strategy(ld_min_r2 = 0.5, theta_extend = 0.3), "ld_min_r2")
  expect_error(graph_test_strategy(ld_min_r2 = 0.9, theta_core = 0.8), "ld_min_r2")
  expect_error(
    graph_test_strategy(ld_min_r2 = 0.5, theta_core_by_chr = c("1" = 0.3)),
    "ld_min_r2"
  )

  # theta_refill reads from the genotypes, not the edge table, so it is exempt
  expect_silent(graph_test_strategy(ld_min_r2 = 0.6, theta_refill = 0.1,
                                    theta_core = 0.8, theta_extend = 0.7,
                                    theta_bridge = 0.7))
})


test_that("graph_strategy warns when window_core reaches past window_ld", {
  # Those pairs were never computed, so the extra reach cannot recover them
  expect_warning(graph_test_strategy(window_ld = 5, window_core = 10), "window_ld")
  expect_silent(graph_test_strategy(window_ld = 5, window_core = 5))
})


test_that("graph_strategy rejects window and count arguments that are not whole numbers", {
  expect_error(graph_test_strategy(window_ld = 0), "window_ld")
  expect_error(graph_test_strategy(window_ld = 2.5), "window_ld")
  expect_error(graph_test_strategy(min_links = 0), "min_links")
  expect_error(graph_test_strategy(min_block_snps = 0), "min_block_snps")
  expect_error(graph_test_strategy(max_gap_snps = -1), "max_gap_snps")
  expect_error(graph_test_strategy(max_gap_position = -5), "max_gap_position")
})


test_that("graph_strategy needs no data, so it can be built before any is loaded", {
  expect_s3_class(graph_strategy(theta_core = 0.5, ld_min_r2 = 0.1), "graph_strategy")
})


# ld_func_c, as the graph method calls it ---------------------------------------------------------------

test_that("ld_func_c returns the package LD column convention", {
  edges <- ld_func_c(graph_chr_fixture(1), window = 3, min_r2 = 0.20,
                          min_obs = 3L)

  expect_equal(colnames(edges), c("Chrom", "Locus1", "Locus2", "Name1", "Name2", "LD"))
  expect_true(all(edges$Locus1 < edges$Locus2))
  expect_true(all(edges$Chrom == 1))
})


test_that("ld_func_c reproduces the reference edge table", {
  edges <- ld_func_c(graph_chr_fixture(1), window = 3, min_r2 = 0.20,
                          min_obs = 3L)
  edges <- as.data.frame(edges)
  edges <- edges[order(edges$Locus1, edges$Locus2), ]

  expect_equal(edges$Locus1, c(1, 1, 2, 2, 4, 5, 5, 7))
  expect_equal(edges$Locus2, c(2, 4, 4, 5, 5, 7, 8, 8))
  expect_equal(round(edges$LD, 4),
               c(0.8972, 0.5047, 0.6311, 1.0000, 0.6311, 0.3505, 0.3505, 1.0000))
  expect_equal(edges$Name1[1], "c1_01")
  expect_equal(edges$Name2[1], "c1_02")
})


test_that("ld_func_c computes r2 as squared pairwise-complete correlation", {
  geno <- graph_chr_fixture(1)
  edges <- ld_func_c(geno, window = 3, min_r2 = 0, min_obs = 3L)

  for (i in seq_len(nrow(edges))) {
    x <- as.numeric(geno[edges$Locus1[i], -(1:3)])
    y <- as.numeric(geno[edges$Locus2[i], -(1:3)])
    expect_equal(edges$LD[i], cor(x, y, use = "pairwise.complete.obs")^2)
  }
})


test_that("ld_func_c only compares markers inside the window", {
  same <- c(0, 1, 2, 0, 1, 2)
  geno <- gb_geno_rows(list(same, same, same, same))

  edges <- ld_func_c(geno, window = 2, min_r2 = 0.1, min_obs = 3L)
  edges <- as.data.frame(edges)
  edges <- edges[order(edges$Locus1, edges$Locus2), ]

  # every pair is perfectly correlated, so only the window decides what appears:
  # 1-4 is three apart and never computed
  expect_equal(paste(edges$Locus1, edges$Locus2), c("1 2", "1 3", "2 3", "2 4", "3 4"))
  expect_true(all(edges$LD == 1))
})


test_that("ld_func_c drops pairs below the floor", {
  edges <- ld_func_c(graph_chr_fixture(1), window = 3, min_r2 = 0.64,
                          min_obs = 3L)

  expect_true(all(edges$LD >= 0.64))
  expect_equal(nrow(edges), 3)   # 0.8972, 1.0000 and 1.0000
})


test_that("ld_func_c skips a pair with fewer than three shared observations", {
  # cor() returns r2 = 1 for two complete observations, which would plant perfect LD
  # on a pair that has no evidence behind it at all
  partial <- c(0, 1, NA, NA, NA, NA)
  geno <- gb_geno_rows(list(partial, partial))

  expect_equal(cor(partial, partial, use = "pairwise.complete.obs")^2, 1)
  expect_equal(nrow(ld_func_c(geno, window = 2, min_r2 = 0.1, min_obs = 3L)), 0)
})


test_that("ld_func_c skips a marker with no variance", {
  geno <- gb_geno_rows(list(rep(1, 6), c(0, 1, 2, 0, 1, 2)))

  expect_equal(nrow(ld_func_c(geno, window = 2, min_r2 = 0.1, min_obs = 3L)), 0)
})


test_that("ld_func_c returns an empty table, with its columns, for a single marker", {
  geno <- gb_geno_rows(list(c(0, 1, 2, 0, 1, 2)))

  edges <- ld_func_c(geno, window = 2, min_r2 = 0.1, min_obs = 3L)

  expect_equal(nrow(edges), 0)
  expect_equal(colnames(edges), c("Chrom", "Locus1", "Locus2", "Name1", "Name2", "LD"))
})


# core_blocks ------------------------------------------------------------------

test_that("core_blocks takes components of at least two markers as blocks", {
  map <- graph_chr_fixture(1, "map")
  edges <- ld_func_c(graph_chr_fixture(1), window = 3, min_r2 = 0.20, min_obs = 3L)

  core <- core_blocks(edges, map, theta_core = 0.80, window_core = 2)

  expect_equal(core$blocks, list(c("c1_01", "c1_02"), c("c1_07", "c1_08")))
  expect_equal(core$unassigned, c("c1_03", "c1_04", "c1_05", "c1_06"))
})


test_that("core_blocks ignores edges below theta_core", {
  map <- gb_map(3)
  edges <- gb_edges(map, list(c(1, 2, 0.5), c(2, 3, 0.9)))

  expect_equal(core_blocks(edges, map, 0.8, 5)$blocks, list(c("m02", "m03")))
  expect_equal(core_blocks(edges, map, 0.4, 5)$blocks, list(c("m01", "m02", "m03")))
})


test_that("core_blocks ignores edges reaching past window_core", {
  map <- gb_map(4)
  edges <- gb_edges(map, list(c(1, 2, 0.9), c(2, 4, 0.9)))

  # the 2-4 edge is two apart, so it joins the component only at window_core 2
  expect_equal(core_blocks(edges, map, 0.8, 1)$blocks, list(c("m01", "m02")))
  expect_equal(core_blocks(edges, map, 0.8, 2)$blocks, list(c("m01", "m02", "m04")))
})


test_that("core_blocks returns everything as unassigned when no edge qualifies", {
  map <- gb_map(3)
  core <- core_blocks(gb_edges(map, list(c(1, 2, 0.5))), map, theta_core = 0.8,
                      window_core = 5)

  expect_equal(core$blocks, list())
  expect_equal(core$unassigned, c("m01", "m02", "m03"))
})


test_that("core_blocks can produce a block that is not contiguous in marker order", {
  # This is the whole reason Stage II exists: a component need not be an interval
  map <- gb_map(3)
  core <- core_blocks(gb_edges(map, list(c(1, 3, 0.9))), map, 0.8, 5)

  expect_equal(core$blocks, list(c("m01", "m03")))
  expect_equal(core$unassigned, "m02")
})


# extend_blocks ----------------------------------------------------------------

test_that("extend_blocks attaches a marker whose edge was too long for the core step", {
  map <- graph_chr_fixture(1, "map")
  edges <- ld_func_c(graph_chr_fixture(1), window = 3, min_r2 = 0.20, min_obs = 3L)
  core <- core_blocks(edges, map, theta_core = 0.80, window_core = 2)

  ext <- extend_blocks(core$blocks, core$unassigned, edges, map,
                       theta_extend = 0.70, window_extend = 3, min_links = 1)

  # c1_05 is a perfect match for c1_02 but three markers away, past window_core
  expect_equal(ext$blocks, list(c("c1_01", "c1_02", "c1_05"), c("c1_07", "c1_08")))
  expect_equal(ext$unassigned, c("c1_03", "c1_04", "c1_06"))
})


test_that("extend_blocks prefers the block with more supporting edges over stronger LD", {
  map <- gb_map(5)
  blocks <- list(c("m01", "m02"), c("m04", "m05"))
  edges <- gb_edges(map, list(c(2, 3, 0.9), c(3, 4, 0.5), c(3, 5, 0.5)))

  ext <- extend_blocks(blocks, "m03", edges, map, 0.3, 5, 1)

  # one edge at 0.9 loses to two edges at 0.5
  expect_equal(ext$blocks, list(c("m01", "m02"), c("m03", "m04", "m05")))
  expect_equal(ext$unassigned, character())
})


test_that("extend_blocks breaks a tie on edge count with mean LD", {
  map <- gb_map(5)
  blocks <- list(c("m01", "m02"), c("m04", "m05"))
  edges <- gb_edges(map, list(c(2, 3, 0.5), c(3, 4, 0.9)))

  ext <- extend_blocks(blocks, "m03", edges, map, 0.3, 5, 1)

  expect_equal(ext$blocks, list(c("m01", "m02"), c("m03", "m04", "m05")))
})


test_that("extend_blocks breaks a tie on edge count and LD with the nearer block", {
  map <- gb_map(8)
  blocks <- list(c("m01", "m02"), c("m06", "m07"))
  edges <- gb_edges(map, list(c(2, 3, 0.5), c(3, 6, 0.5)))

  ext <- extend_blocks(blocks, "m03", edges, map, 0.3, 5, 1)

  # one marker from block 1, three from block 2
  expect_equal(ext$blocks, list(c("m01", "m02", "m03"), c("m06", "m07")))
})


test_that("extend_blocks measures distance to the nearest member, not the block span", {
  map <- gb_map(8)
  edges <- gb_edges(map, list(c(3, 6, 0.5)))

  expect_equal(extend_blocks(list(c("m06", "m07")), "m03", edges, map, 0.3, 2, 1)$unassigned,
               "m03")
  expect_equal(extend_blocks(list(c("m06", "m07")), "m03", edges, map, 0.3, 3, 1)$blocks,
               list(c("m03", "m06", "m07")))
})


test_that("extend_blocks honours min_links", {
  map <- gb_map(5)
  blocks <- list(c("m01", "m02"), c("m04", "m05"))
  edges <- gb_edges(map, list(c(2, 3, 0.5), c(3, 4, 0.9)))

  ext <- extend_blocks(blocks, "m03", edges, map, 0.3, 5, min_links = 2)

  expect_equal(ext$blocks, blocks)
  expect_equal(ext$unassigned, "m03")
})


test_that("extend_blocks is sequential, so a marker it places can carry the next one", {
  map <- gb_map(5)
  edges <- gb_edges(map, list(c(2, 3, 0.9), c(3, 4, 0.9)))

  # m04 has no edge to the block at all, only to m03, which joins first
  ext <- extend_blocks(list(c("m01", "m02")), c("m03", "m04"), edges, map, 0.3, 5, 1)
  expect_equal(ext$blocks, list(c("m01", "m02", "m03", "m04")))

  # without m03 in the queue there is nothing for m04 to reach
  ext <- extend_blocks(list(c("m01", "m02")), "m04", edges, map, 0.3, 5, 1)
  expect_equal(ext$blocks, list(c("m01", "m02")))
  expect_equal(ext$unassigned, "m04")
})


test_that("extend_blocks leaves a marker with no qualifying edge alone", {
  map <- gb_map(4)
  edges <- gb_edges(map, list(c(1, 2, 0.9), c(2, 3, 0.2)))

  ext <- extend_blocks(list(c("m01", "m02")), c("m03", "m04"), edges, map, 0.5, 5, 1)

  expect_equal(ext$blocks, list(c("m01", "m02")))
  expect_equal(ext$unassigned, c("m03", "m04"))
})


# bridge_blocks ----------------------------------------------------------------

test_that("bridge_blocks merges two blocks across a small gap", {
  map <- graph_chr_fixture(1, "map")
  edges <- ld_func_c(graph_chr_fixture(1), window = 3, min_r2 = 0.20, min_obs = 3L)
  blocks <- list(c("c1_01", "c1_02", "c1_05"), c("c1_07", "c1_08"))

  bridged <- bridge_blocks(blocks, edges, map, theta_bridge = 0.20, max_gap_snps = 2)

  # c1_06 lies in the gap and is NOT absorbed, so the block skips it
  expect_equal(bridged, list(c("c1_01", "c1_02", "c1_05", "c1_07", "c1_08")))
})


test_that("bridge_blocks refuses a gap wider than max_gap_snps", {
  map <- gb_map(8)
  blocks <- list(c("m01", "m02"), c("m06", "m07"))
  edges <- gb_edges(map, list(c(2, 6, 0.9)))

  expect_equal(bridge_blocks(blocks, edges, map, 0.2, max_gap_snps = 2), blocks)
  expect_equal(bridge_blocks(blocks, edges, map, 0.2, max_gap_snps = 3),
               list(c("m01", "m02", "m06", "m07")))
})


test_that("bridge_blocks refuses a gap with no qualifying boundary edge", {
  map <- gb_map(8)
  blocks <- list(c("m01", "m02"), c("m06", "m07"))
  edges <- gb_edges(map, list(c(2, 6, 0.9)))

  expect_equal(bridge_blocks(blocks, edges, map, theta_bridge = 0.95, max_gap_snps = 3),
               blocks)
})


test_that("bridge_blocks only looks at the three markers facing the gap", {
  map <- gb_map(10)
  blocks <- list(c("m01", "m02", "m03", "m04"), c("m06", "m07"))

  # m01 is four markers back from the gap, so its edge does not count
  expect_equal(bridge_blocks(blocks, gb_edges(map, list(c(1, 6, 0.9))), map, 0.2, 2),
               blocks)
  # m04 faces the gap, so the same edge from there does
  expect_equal(bridge_blocks(blocks, gb_edges(map, list(c(4, 6, 0.9))), map, 0.2, 2),
               list(c("m01", "m02", "m03", "m04", "m06", "m07")))
})


test_that("bridge_blocks leaves a single block alone", {
  map <- gb_map(4)
  blocks <- list(c("m01", "m02"))

  expect_equal(bridge_blocks(blocks, gb_edges(map, list(c(1, 2, 0.9))), map, 0.2, 2),
               blocks)
})


# graph_chromosome_blocks ------------------------------------------------------

test_that("graph_chromosome_blocks runs Stage I and keeps every marker", {
  blocks <- graph_chromosome_blocks(graph_chr_fixture(1), graph_chr_fixture(1, "map"),
                                    graph_test_strategy())

  expect_equal(blocks, list(
    c("c1_01", "c1_02", "c1_05", "c1_07", "c1_08"),
    "c1_03",
    "c1_04",
    "c1_06"
  ))
  expect_setequal(unlist(blocks), graph_chr_fixture(1, "map")$SNP)
})


test_that("graph_chromosome_blocks refuses a bridge with no boundary edge", {
  blocks <- graph_chromosome_blocks(graph_chr_fixture(2), graph_chr_fixture(2, "map"),
                                    graph_test_strategy())

  # c2_01/c2_03 and c2_05/c2_06 are one marker apart but share no LD
  expect_equal(blocks, list(
    c("c2_01", "c2_03"),
    "c2_02",
    "c2_04",
    c("c2_05", "c2_06"),
    "c2_07"
  ))
})


# linearise_blocks -------------------------------------------------------------

test_that("linearise_blocks describes a segment with its span, density and score", {
  blocks <- list(c("m01", "m02", "m05"))

  segs <- as.data.frame(linearise_blocks(blocks, gb_map(10), max_gap_markers = 3,
                                         max_gap_position = NULL, min_block_snps = 2))

  expect_equal(nrow(segs), 1)
  expect_equal(segs$first_index, 1)
  expect_equal(segs$last_index, 5)
  expect_equal(segs$first_position, 100)
  expect_equal(segs$last_position, 500)
  expect_equal(segs$length, 400)
  expect_equal(segs$n_markers, 3)
  expect_equal(segs$span, 5)
  expect_equal(segs$density, 3 / 5)
  expect_equal(segs$score, 1000 * 3 + 100 * (3 / 5) - log1p(400))
  expect_equal(segs$markers, "m01;m02;m05")
})


test_that("linearise_blocks splits a block where its markers are too far apart", {
  blocks <- list(c("m01", "m02", "m05"))

  wide <- as.data.frame(linearise_blocks(blocks, gb_map(10), 3, NULL, 1))
  expect_equal(wide$markers, "m01;m02;m05")

  narrow <- as.data.frame(linearise_blocks(blocks, gb_map(10), 2, NULL, 1))
  expect_equal(narrow$markers, c("m01;m02", "m05"))
  expect_equal(narrow$first_index, c(1, 5))
  expect_equal(narrow$last_index, c(2, 5))
})


test_that("linearise_blocks splits on map distance when max_gap_position is set", {
  map <- data.frame(SNP = c("m01", "m02", "m03"), Chromosome = 1,
                    Position = c(100, 200, 5000), stringsAsFactors = FALSE)
  blocks <- list(c("m01", "m02", "m03"))

  expect_equal(as.data.frame(linearise_blocks(blocks, map, 5, NULL, 1))$markers,
               "m01;m02;m03")
  expect_equal(as.data.frame(linearise_blocks(blocks, map, 5, 1000, 1))$markers,
               c("m01;m02", "m03"))
})


test_that("linearise_blocks drops segments below min_block_snps", {
  blocks <- list(c("m01", "m02", "m05"))

  segs <- as.data.frame(linearise_blocks(blocks, gb_map(10), 2, NULL, min_block_snps = 2))

  expect_equal(segs$markers, "m01;m02")
})


test_that("linearise_blocks returns no rows when every segment is too small", {
  segs <- linearise_blocks(list("m01", "m02"), gb_map(10), 2, NULL, 2)

  expect_equal(nrow(segs), 0)
})


# select_nonoverlapping --------------------------------------------------------

test_that("select_nonoverlapping rejects a candidate overlapping one already taken", {
  segs <- data.frame(
    source_block = c("b1", "b2", "b3"), chromosome = 1,
    first_index = c(1, 4, 8), last_index = c(5, 7, 9),
    first_position = c(100, 400, 800), last_position = c(500, 700, 900),
    n_markers = c(4, 3, 2), markers = c("m01;m02;m03;m05", "m04;m06;m07", "m08;m09"),
    length = c(400, 300, 100), span = c(5, 4, 2), density = c(4 / 5, 3 / 4, 1),
    stringsAsFactors = FALSE
  )
  segs$score <- 1000 * segs$n_markers + 100 * segs$density - log1p(segs$length)

  kept <- as.data.frame(select_nonoverlapping(segs))

  # b1 scores highest and claims 1 to 5; b2 starts at 4 and is refused; b3 is clear
  expect_equal(kept$source_block, c("b1", "b3"))
  expect_equal(kept$first_index, c(1, 8))
})


test_that("select_nonoverlapping claims the whole span, holes included", {
  segs <- data.frame(
    source_block = c("wide", "inner"), chromosome = 1,
    first_index = c(1, 3), last_index = c(5, 4),
    first_position = c(100, 300), last_position = c(500, 400),
    n_markers = c(3, 2), markers = c("m01;m02;m05", "m03;m04"),
    length = c(400, 100), span = c(5, 2), density = c(3 / 5, 1),
    stringsAsFactors = FALSE
  )
  segs$score <- 1000 * segs$n_markers + 100 * segs$density - log1p(segs$length)

  kept <- as.data.frame(select_nonoverlapping(segs))

  # "inner" covers the markers "wide" skips, but sits inside the span it claimed
  expect_equal(kept$source_block, "wide")
})


test_that("select_nonoverlapping returns its result in map order", {
  segs <- data.frame(
    source_block = c("late", "early"), chromosome = 1,
    first_index = c(10, 1), last_index = c(12, 3),
    first_position = c(1000, 100), last_position = c(1200, 300),
    n_markers = c(3, 2), markers = c("m10;m11;m12", "m01;m03"),
    length = c(200, 200), span = c(3, 3), density = c(1, 2 / 3),
    stringsAsFactors = FALSE
  )
  segs$score <- 1000 * segs$n_markers + 100 * segs$density - log1p(segs$length)

  expect_equal(as.data.frame(select_nonoverlapping(segs))$source_block,
               c("early", "late"))
})


test_that("select_nonoverlapping handles an empty candidate table", {
  segs <- data.frame(
    source_block = character(), chromosome = numeric(),
    first_index = numeric(), last_index = numeric(),
    first_position = numeric(), last_position = numeric(),
    n_markers = numeric(), markers = character(), length = numeric(),
    span = numeric(), density = numeric(), score = numeric(),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(select_nonoverlapping(segs)), 0)
})


# block_ld_support -------------------------------------------------------------

test_that("block_ld_support is TRUE when any one member reaches the threshold", {
  geno <- gb_geno_rows(list(c(0, 0, 1, 1, 2, 2), c(0, 0, 1, 1, 2, 2), c(2, 1, 0, 2, 0, 1)))

  # m02 matches perfectly, m03 does not: one qualifying member is enough
  expect_true(block_ld_support(geno, "m01", c("m02", "m03"), 0.8))
  expect_false(block_ld_support(geno, "m01", "m03", 0.8))
})


test_that("block_ld_support ignores the marker itself", {
  geno <- gb_geno_rows(list(c(0, 0, 1, 1, 2, 2), c(2, 1, 0, 2, 0, 1)))

  expect_false(block_ld_support(geno, "m01", "m01", 0.8))
})


test_that("block_ld_support is FALSE when there is nothing to compare against", {
  geno <- gb_geno_rows(list(c(0, 0, 1, 1, 2, 2), c(0, 0, 1, 1, 2, 2)))

  expect_false(block_ld_support(geno, "m01", character(), 0.8))
  expect_false(block_ld_support(geno, "absent", c("m01", "m02"), 0.8))
})


test_that("block_ld_support reaches past the LD window", {
  # It recomputes from the genotypes, so a member ld_func_c() never compared
  # against still counts. This is what lets the refill step find support that
  # Stage I could not see.
  same <- c(0, 0, 1, 1, 2, 2)
  geno <- gb_geno_rows(list(same, c(2, 1, 0, 2, 0, 1), c(1, 0, 2, 1, 0, 2), same))

  expect_equal(nrow(ld_func_c(geno, window = 1, min_r2 = 0.8, min_obs = 3L)), 0)
  expect_true(block_ld_support(geno, "m01", "m04", 0.8))
})


# refill_internal_markers ------------------------------------------------------

test_that("refill_internal_markers keeps a marker with strong LD to the block", {
  map  <- graph_chr_fixture(1, "map")
  geno <- graph_chr_fixture(1)
  blocks <- graph_chromosome_blocks(geno, map, graph_test_strategy())
  selected <- select_nonoverlapping(linearise_blocks(blocks, map, 3, NULL, 2))

  res <- refill_internal_markers(selected, map, geno, theta_refill = 0.50)

  # c1_04 reaches 0.63 against c1_02; c1_03 and c1_06 are noise
  expect_equal(res$refilled$SNP, "c1_04")
  expect_equal(res$dropped$SNP, c("c1_03", "c1_06"))
  expect_equal(as.data.frame(res$blocks)$markers,
               "c1_01;c1_02;c1_04;c1_05;c1_07;c1_08")
  expect_equal(as.data.frame(res$blocks)$n_markers, 6)
  expect_equal(res$external, character())
})


test_that("refill_internal_markers leaves the index bounds alone", {
  map  <- graph_chr_fixture(1, "map")
  geno <- graph_chr_fixture(1)
  blocks <- graph_chromosome_blocks(geno, map, graph_test_strategy())
  selected <- select_nonoverlapping(linearise_blocks(blocks, map, 3, NULL, 2))

  res <- refill_internal_markers(selected, map, geno, 0.50)
  refilled <- as.data.frame(res$blocks)

  # a refilled marker was inside the span already, so only the density moves
  expect_equal(refilled$first_index, as.data.frame(selected)$first_index)
  expect_equal(refilled$last_index, as.data.frame(selected)$last_index)
  expect_equal(refilled$span, 8)
  expect_equal(refilled$density, 6 / 8)
})


test_that("refill_internal_markers reports markers outside every block as external", {
  map  <- graph_chr_fixture(2, "map")
  geno <- graph_chr_fixture(2)
  blocks <- graph_chromosome_blocks(geno, map, graph_test_strategy())
  selected <- select_nonoverlapping(linearise_blocks(blocks, map, 3, NULL, 2))

  res <- refill_internal_markers(selected, map, geno, theta_refill = 0.50)

  # c2_02 sits inside the 1 to 3 span, c2_04 and c2_07 sit outside both blocks
  expect_equal(res$refilled$SNP, "c2_02")
  expect_equal(res$external, c("c2_04", "c2_07"))
  expect_equal(nrow(res$dropped), 0)
  expect_equal(as.data.frame(res$blocks)$markers,
               c("c2_01;c2_02;c2_03", "c2_05;c2_06"))
})


test_that("refill_internal_markers tests a marker only against the block enclosing it", {
  map  <- graph_chr_fixture(2, "map")
  geno <- graph_chr_fixture(2)
  blocks <- graph_chromosome_blocks(geno, map, graph_test_strategy())
  selected <- select_nonoverlapping(linearise_blocks(blocks, map, 3, NULL, 2))

  res <- refill_internal_markers(selected, map, geno, theta_refill = 0.50)

  expect_equal(res$refilled$block, 1)
})


# perform_graph_blocking -------------------------------------------------------

test_that("perform_graph_blocking returns one entry per chromosome", {
  blocks <- perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(),
                                   graph_test_strategy())

  expect_type(blocks, "list")
  expect_named(blocks, c("1", "2"))
})


test_that("perform_graph_blocking reproduces the reference blocking", {
  blocks <- perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(),
                                   graph_test_strategy())

  expect_equal(blocks[["1"]],
               list(c("c1_01", "c1_02", "c1_04", "c1_05", "c1_07", "c1_08")))
  expect_equal(blocks[["2"]], list(
    c("c2_01", "c2_02", "c2_03"),
    "c2_04",
    c("c2_05", "c2_06"),
    "c2_07"
  ))
})


test_that("perform_graph_blocking drops markers the refill step refused", {
  # Unlike the window and LD methods, this one is not a partition of the map:
  # c1_03 and c1_06 sit inside a block's span with no LD to it, and are lost
  blocks <- perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(),
                                   graph_test_strategy())
  placed <- unlist(blocks, use.names = FALSE)

  expect_false(any(c("c1_03", "c1_06") %in% placed))
  expect_equal(length(placed), nrow(graph_map_fixture()) - 2)
})


test_that("perform_graph_blocking orders blocks and their markers by position", {
  blocks <- perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(),
                                   graph_test_strategy())
  map <- graph_map_fixture()

  for (chr in names(blocks)) {
    idx <- stats::setNames(seq_len(sum(map$Chromosome == chr)),
                           map$SNP[map$Chromosome == chr])
    firsts <- vapply(blocks[[chr]], function(b) idx[[b[1]]], numeric(1))
    expect_false(is.unsorted(firsts))
    for (b in blocks[[chr]]) expect_false(is.unsorted(idx[b]))
  }
})


test_that("perform_graph_blocking reports per-chromosome diagnostics", {
  blocks <- perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(),
                                   graph_test_strategy())
  diagnostics <- as.data.frame(attr(blocks, "graph_diagnostics"))

  expect_equal(diagnostics$chromosome, c(1, 2))
  expect_equal(diagnostics$theta_core, c(0.80, 0.80))
  expect_equal(diagnostics$n_blocks, c(1, 4))
  expect_equal(diagnostics$n_multi_marker_blocks, c(1, 2))
  expect_equal(diagnostics$n_single_marker_blocks, c(0, 2))
  expect_equal(diagnostics$n_markers_in, c(8, 7))
  expect_equal(diagnostics$n_markers_out, c(6, 7))
  expect_equal(diagnostics$n_refilled, c(1, 1))
  expect_equal(diagnostics$n_dropped, c(2, 0))
})


test_that("perform_graph_blocking applies a per-chromosome theta_core override", {
  # At 0.99 the only core edge left on chromosome 1 is the identical c1_07/c1_08
  # pair, so its single block falls apart into seven. Chromosome 2 is untouched.
  strategy <- graph_test_strategy(theta_core_by_chr = c("1" = 0.99))
  blocks <- perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(), strategy)

  expect_equal(blocks[["1"]], c(as.list(sprintf("c1_%02d", 1:6)),
                                list(c("c1_07", "c1_08"))))
  expect_equal(blocks[["2"]], list(
    c("c2_01", "c2_02", "c2_03"),
    "c2_04",
    c("c2_05", "c2_06"),
    "c2_07"
  ))

  diagnostics <- as.data.frame(attr(blocks, "graph_diagnostics"))
  expect_equal(diagnostics$theta_core, c(0.99, 0.80))
  expect_equal(diagnostics$n_blocks, c(7, 4))
})


test_that("perform_graph_blocking is deterministic", {
  first  <- perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(),
                                   graph_test_strategy())
  second <- perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(),
                                   graph_test_strategy())

  expect_identical(first, second)
})


test_that("perform_graph_blocking output feeds the shared block table helpers", {
  blocks <- perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(),
                                   graph_test_strategy())

  block_df <- block_obj_to_df(blocks, graph_map_fixture())

  expect_equal(nrow(block_df), 5)
  expect_equal(sum(block_df$Num_SNP), 13)
  expect_equal(block_summary(block_df)$Max_SNP_per_Block, 6)
  expect_equal(block_summary(block_df)$Singleton_Blocks, 2)
})


test_that("perform_graph_blocking leaves a chromosome with no LD as one-marker blocks", {
  geno <- gb_geno(list(
    list("x1", 1, 100, gb_base),
    list("x2", 1, 200, gb_n1),
    list("x3", 1, 300, gb_b2)
  ))

  strategy <- graph_strategy(window_ld = 2, window_core = 2)
  blocks <- perform_graph_blocking(geno, geno[, 1:3], strategy)

  expect_equal(blocks[["1"]], list("x1", "x2", "x3"))
  expect_equal(as.data.frame(attr(blocks, "graph_diagnostics"))$n_single_marker_blocks, 3)
})


# perform_graph_blocking: the geno/map boundary --------------------------------
#
# Genotypes and map arrive as separate arguments, so this is the only place they
# can be checked against each other. Each of these fails silently otherwise.

test_that("perform_graph_blocking needs genotypes", {
  expect_error(
    perform_graph_blocking(NULL, graph_map_fixture(), graph_test_strategy()),
    "geno"
  )
})


test_that("perform_graph_blocking rejects genotypes not in the HapSelect layout", {
  expect_error(
    perform_graph_blocking("not a data frame", graph_map_fixture(),
                           graph_test_strategy()),
    "data frame"
  )
  expect_error(
    perform_graph_blocking(data.frame(a = 1, b = 2, c = 3), graph_map_fixture(),
                           graph_test_strategy()),
    "data frame"
  )
})


test_that("perform_graph_blocking aligns genotypes that arrive in another order", {
  # A marker's index is its row in the map, and its dosages are read from the
  # matching row of geno. Rather than demand the caller line them up, the two are
  # matched by marker name here, which is what the reference implementation does
  # when it builds both from the same VCF.
  geno <- graph_geno_fixture()
  map  <- graph_map_fixture()

  shuffled <- geno[rev(seq_len(nrow(geno))), ]
  rownames(shuffled) <- NULL

  expect_equal(perform_graph_blocking(shuffled, map, graph_test_strategy()),
               perform_graph_blocking(geno, map, graph_test_strategy()))
})


test_that("perform_graph_blocking rejects a map that disagrees with the genotypes", {
  geno <- graph_geno_fixture()
  map  <- graph_map_fixture()

  # Matching by name resolves a different row order, but not a different marker set
  expect_error(perform_graph_blocking(geno, map[-1, ], graph_test_strategy()),
               "same markers")

  renamed <- map
  renamed$SNP[1] <- "somewhere_else"
  expect_error(perform_graph_blocking(geno, renamed, graph_test_strategy()),
               "same markers")
})


test_that("perform_graph_blocking rejects a map that is not in order_map order", {
  # Aligning geno to the map cannot rescue a map that is out of order: the index
  # every gap and window is measured in is the map row itself.
  geno <- graph_geno_fixture()
  map  <- graph_map_fixture()

  reordered <- map[c(2, 1, seq(3, nrow(map))), ]
  expect_error(perform_graph_blocking(geno, reordered, graph_test_strategy()),
               "order_map")
})


test_that("perform_graph_blocking rejects a theta_core_by_chr naming no chromosome", {
  # A name that matches nothing is a no-op: the override never fires and the run
  # quietly uses the global theta_core instead
  strategy <- graph_test_strategy(theta_core_by_chr = c("7" = 0.9))

  expect_error(
    perform_graph_blocking(graph_geno_fixture(), graph_map_fixture(), strategy),
    "theta_core_by_chr"
  )
})
# Reuse of the package's VCF reader --------------------------------------------

test_that("read_vcf_geno and order_map reproduce the prototype's VCF extraction", {
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path))
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tInd1\tInd2\tInd3",
    "chr1A\t300\tsnpB\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1",
    "chr1A\t100\tsnpA\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|1\t./.",
    "chr1B\t200\tsnpC\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0"
  ), path)

  geno <- read_vcf_geno(path)
  chr1a <- geno[geno$Chromosome == "chr1A", ]
  chr1a <- chr1a[order(chr1a$Position), ]

  # What extract_chr_geno_map() produced for chr1A: markers sorted by position,
  # dosages counting ALT alleles, a missing call left as NA
  expect_equal(chr1a$SNP, c("snpA", "snpB"))
  expect_equal(chr1a$Position, c(100, 300))
  expect_equal(unname(as.matrix(chr1a[, -(1:3)])),
               matrix(c(2, 1, NA, 0, 1, 2), nrow = 2, byrow = TRUE))

  # order_map() supplies the marker index the graph functions count on, and
  # numbers the chromosome labels on the way through
  map <- order_map(geno[, 1:3])
  expect_equal(map$SNP, c("snpA", "snpB", "snpC"))
  expect_equal(map$Chromosome, c(1, 1, 2))
})


test_that("a VCF reaches def_blocks through order_map and order_geno", {
  # The path the reference implementation had built in: it read a VCF and produced
  # the genotypes and the map together, from the one source. Here the reader, the
  # map and the alignment are three separate package functions, so this is the test
  # that they still compose into that path.
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path))

  set.seed(4)
  n_ind   <- 30
  samples <- paste0("I", seq_len(n_ind))
  calls   <- function(k) {
    g <- c("0|0", "0|1", "1|1")[sample(3, n_ind, TRUE)]
    if (k > 0) g[sample(n_ind, k)] <- "./."
    paste(g, collapse = "\t")
  }

  # Deliberately out of position order within chr1A, as a VCF may well be
  rows <- c(
    paste("chr1A", 500, "s5", "A", "G", ".", "PASS", ".", "GT", calls(1), sep = "\t"),
    paste("chr1A", 100, "s1", "A", "G", ".", "PASS", ".", "GT", calls(0), sep = "\t"),
    paste("chr1A", 300, "s3", "A", "G", ".", "PASS", ".", "GT", calls(0), sep = "\t"),
    paste("chr1A", 200, "s2", "A", "G", ".", "PASS", ".", "GT", calls(0), sep = "\t"),
    paste("chr1A", 400, "s4", "A", "G", ".", "PASS", ".", "GT", calls(2), sep = "\t"),
    paste("chr1B", 100, "t1", "A", "G", ".", "PASS", ".", "GT", calls(0), sep = "\t"),
    paste("chr1B", 200, "t2", "A", "G", ".", "PASS", ".", "GT", calls(0), sep = "\t"),
    paste("chr1B", 300, "t3", "A", "G", ".", "PASS", ".", "GT", calls(0), sep = "\t")
  )
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
            samples), collapse = "\t"),
    rows
  ), path)

  geno <- read_vcf_geno(path)

  # As read, the genotypes are in file order with the VCF's text chromosome labels,
  # so they are not yet something the blocking functions accept
  expect_equal(geno$SNP[1], "s5")
  expect_type(geno$Chromosome, "character")
  expect_error(check_ld_matrix(geno), "numeric chromosome")

  map     <- order_map(geno[, 1:3], verbose = FALSE)
  aligned <- order_geno(geno, map)

  expect_equal(map$SNP, c("s1", "s2", "s3", "s4", "s5", "t1", "t2", "t3"))
  expect_equal(aligned$SNP, map$SNP)
  expect_silent(check_ld_matrix(aligned))

  blocks <- def_blocks(graph_strategy(window_ld = 3, window_core = 3),
                       map, geno = aligned)

  expect_named(blocks, c("1", "2"))

  # Every marker is accounted for, either in a block or deliberately dropped by the
  # refill step, and nothing is invented
  placed <- unlist(blocks, use.names = FALSE)
  expect_true(all(placed %in% map$SNP))
  expect_false(anyDuplicated(placed) > 0)
})


test_that("def_blocks takes the genotypes as read, without pre-alignment", {
  # order_geno() is preparation a caller may reasonably forget, and forgetting it
  # used to be an error. perform_graph_blocking() aligns by marker name itself, so
  # the unaligned genotypes give the same answer as the aligned ones.
  geno <- graph_geno_fixture()
  map  <- graph_map_fixture()

  shuffled <- geno[sample(nrow(geno)), ]
  rownames(shuffled) <- NULL

  expect_equal(def_blocks(graph_test_strategy(), map, geno = shuffled),
               def_blocks(graph_test_strategy(), map, geno = geno))
})
