# Fixtures --------------------------------------------------------------------
#
# Marker spacing is uneven so that the two methods give different partitions,
# and so that an empty grid interval ([400,600) on chr 1) and a marker sitting
# exactly on a boundary (c6 at 800) are both covered.
#
#   Chr 1: 100, 250, 300, 700, 750, 800, 1500
#   Chr 2:  50,  60, 900

make_window_map_fixture <- function() {
  data.frame(
    SNP        = c("c1", "c2", "c3", "c4", "c5", "c6", "c7", "d1", "d2", "d3"),
    Chromosome = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L),
    Position   = c(100, 250, 300, 700, 750, 800, 1500, 50, 60, 900),
    stringsAsFactors = FALSE
  )
}


# Tests: structure ------------------------------------------------------------

test_that("perform_window_blocking returns a named list with one entry per chromosome", {
  map <- make_window_map_fixture()

  blocks <- perform_window_blocking(map, window = 3, method = "window_snp")

  expect_type(blocks, "list")
  expect_named(blocks, c("1", "2"))
  expect_length(blocks[["1"]], 3)   # 7 markers in runs of 3 -> 3, 3, 1
  expect_length(blocks[["2"]], 1)   # 3 markers in runs of 3 -> 1 full block
})


test_that("perform_window_blocking output is a partition of the map", {
  map <- make_window_map_fixture()

  for (args in list(list(window = 3,   method = "window_snp"),
                    list(window = 200, method = "window_map"))) {

    blocks <- perform_window_blocking(map, window = args$window, method = args$method)
    assigned <- unlist(blocks, use.names = FALSE)

    # every marker appears in exactly one block
    expect_setequal(assigned, map$SNP)
    expect_equal(length(assigned), nrow(map))
  }
})


# Tests: window_snp -----------------------------------------------------------

test_that("window_snp cuts consecutive runs of markers and keeps a short tail", {
  map <- make_window_map_fixture()

  blocks <- perform_window_blocking(map, window = 3, method = "window_snp")

  chr1 <- blocks[["1"]]
  chr2 <- blocks[["2"]]

  expect_equal(chr1[[1]], c("c1", "c2", "c3"))
  expect_equal(chr1[[2]], c("c4", "c5", "c6"))
  expect_equal(chr1[[3]], "c7")            # 7 %% 3 == 1, tail kept as its own block

  expect_equal(chr2[[1]], c("d1", "d2", "d3"))  # exact multiple, no tail
})


test_that("window_snp = 1 makes every marker a singleton block", {
  map <- make_window_map_fixture()

  blocks <- perform_window_blocking(map, window = 1, method = "window_snp")

  expect_length(blocks[["1"]], 7)
  expect_length(blocks[["2"]], 3)
  expect_true(all(lengths(blocks[["1"]]) == 1))
  expect_equal(blocks[["1"]][[1]], "c1")
})


test_that("window_snp larger than the chromosome gives a single block", {
  map <- make_window_map_fixture()

  blocks <- perform_window_blocking(map, window = 100, method = "window_snp")

  expect_length(blocks[["1"]], 1)
  expect_equal(blocks[["1"]][[1]], c("c1", "c2", "c3", "c4", "c5", "c6", "c7"))
  expect_equal(blocks[["2"]][[1]], c("d1", "d2", "d3"))
})


# Tests: window_map -----------------------------------------------------------

test_that("window_map assigns markers to a fixed grid anchored at 0", {
  map <- make_window_map_fixture()

  blocks <- perform_window_blocking(map, window = 200, method = "window_map")

  chr1 <- blocks[["1"]]
  chr2 <- blocks[["2"]]

  expect_length(chr1, 5)
  expect_equal(chr1[[1]], "c1")                # [0,200)
  expect_equal(chr1[[2]], c("c2", "c3"))       # [200,400)
  expect_equal(chr1[[3]], c("c4", "c5"))       # [600,800)
  expect_equal(chr1[[4]], "c6")                # [800,1000)
  expect_equal(chr1[[5]], "c7")                # [1400,1600)

  expect_length(chr2, 2)
  expect_equal(chr2[[1]], c("d1", "d2"))       # [0,200)
  expect_equal(chr2[[2]], "d3")                # [800,1000)
})


test_that("window_map skips grid intervals that contain no markers", {
  map <- make_window_map_fixture()

  blocks <- perform_window_blocking(map, window = 200, method = "window_map")

  # [400,600) is empty on chr 1 and must not appear as an empty block
  expect_false(any(lengths(blocks[["1"]]) == 0))

  # chr 2 spans 50..900, which is five grid intervals, but only two hold markers
  expect_length(blocks[["2"]], 2)
})


test_that("window_map intervals are half-open, so a marker on a boundary opens a new block", {
  map <- make_window_map_fixture()

  blocks <- perform_window_blocking(map, window = 200, method = "window_map")

  # c6 is at exactly 800; it belongs to [800,1000), not to [600,800) with c4/c5
  expect_equal(blocks[["1"]][[3]], c("c4", "c5"))
  expect_equal(blocks[["1"]][[4]], "c6")
})


test_that("window_map accepts a fractional window for cM maps", {
  map <- data.frame(
    SNP        = c("g1", "g2", "g3", "g4"),
    Chromosome = rep(1L, 4),
    Position   = c(0.1, 0.4, 0.6, 2.7),
    stringsAsFactors = FALSE
  )

  blocks <- perform_window_blocking(map, window = 0.5, method = "window_map")

  expect_length(blocks[["1"]], 3)
  expect_equal(blocks[["1"]][[1]], c("g1", "g2"))  # [0.0,0.5)
  expect_equal(blocks[["1"]][[2]], "g3")           # [0.5,1.0)
  expect_equal(blocks[["1"]][[3]], "g4")           # [2.5,3.0)
})


# Tests: input handling -------------------------------------------------------

test_that("perform_window_blocking warns and sorts a map that is not ordered by position", {
  map <- make_window_map_fixture()

  # An unsorted map is wrong rather than an error: window_snp would group
  # non-neighbours, and window_map would order markers within a block
  # arbitrarily.
  for (reordered in list(map[c(7, 2, 9, 4, 1, 10, 5, 3, 8, 6), ],
                         map[nrow(map):1, ],
                         map[c(2, 1, 3:10), ])) {

    expect_warning(snp <- perform_window_blocking(reordered, window = 3, method = "window_snp"),
                   "not sorted by Position")
    expect_warning(mp  <- perform_window_blocking(reordered, window = 200, method = "window_map"),
                   "not sorted by Position")

    # repaired output must match what the already-sorted map produces
    expect_equal(snp, perform_window_blocking(map, window = 3,   method = "window_snp"))
    expect_equal(mp,  perform_window_blocking(map, window = 200, method = "window_map"))
  }

  # the warning names the offending marker and points at the fix
  expect_warning(perform_window_blocking(map[c(2, 1, 3:10), ], window = 3, method = "window_snp"),
                 "'c1'")
  expect_warning(perform_window_blocking(map[c(2, 1, 3:10), ], window = 3, method = "window_snp"),
                 "order_map")

  # marker order within a block drives First_SNP/Last_SNP downstream
  suppressWarnings(
    blocks <- perform_window_blocking(map[nrow(map):1, ], window = 200, method = "window_map"))
  expect_true(all(block_obj_to_df(blocks, map)$Block_Size >= 0))
})


test_that("column order in the map does not matter", {
  # columns are located by name, never by position
  map <- make_window_map_fixture()
  reordered_cols <- map[nrow(map):1, c("Chromosome", "Position", "SNP")]

  expect_warning(blocks <- perform_window_blocking(reordered_cols, window = 3,
                                             method = "window_snp"),
                 "not sorted by Position")
  expect_equal(blocks, perform_window_blocking(map, window = 3, method = "window_snp"))
})


test_that("perform_window_blocking accepts a correctly sorted map", {
  map <- make_window_map_fixture()

  # equal positions are legitimate and must not trip the sort check
  tied <- map
  tied$Position[tied$SNP == "c3"] <- 250   # c2 and c3 now share a position
  expect_silent(perform_window_blocking(tied, window = 3, method = "window_snp"))

  # chromosomes need not be contiguous, since blocking splits on them anyway
  interleaved <- map[order(map$Position), ]
  expect_silent(perform_window_blocking(interleaved, window = 3, method = "window_snp"))

  # a single-marker map has no adjacent pairs to compare
  expect_silent(perform_window_blocking(map[1, ], window = 3, method = "window_snp"))

  # a missing position must not be mistaken for an ordering violation
  na_map <- map
  na_map$Position[na_map$SNP == "c3"] <- NA
  expect_warning(perform_window_blocking(na_map, window = 3, method = "window_snp"), "dropped")
})


test_that("perform_window_blocking requires the order_map() column names", {
  map <- make_window_map_fixture()

  # the map is the only source of chromosome structure here, so a missing or
  # differently-spelled column must fail up front
  alt <- map
  colnames(alt)[colnames(alt) == "Chromosome"] <- "Chrom"
  expect_error(perform_window_blocking(alt, window = 3, method = "window_snp"), "order_map")

  no_pos <- map
  colnames(no_pos)[colnames(no_pos) == "Position"] <- "pos"
  expect_error(perform_window_blocking(no_pos, window = 3, method = "window_snp"), "order_map")

  chr_txt <- map
  chr_txt$Position <- as.character(chr_txt$Position)
  expect_error(perform_window_blocking(chr_txt, window = 3, method = "window_snp"), "numeric")
})


test_that("perform_window_blocking drops markers with missing positions and warns", {
  map <- make_window_map_fixture()
  map$Position[map$SNP == "c3"] <- NA

  expect_warning(
    blocks <- perform_window_blocking(map, window = 3, method = "window_snp"),
    "dropped"
  )

  expect_false("c3" %in% unlist(blocks, use.names = FALSE))
  expect_equal(length(unlist(blocks, use.names = FALSE)), nrow(map) - 1)
  expect_equal(blocks[["1"]][[1]], c("c1", "c2", "c4"))
})


test_that("perform_window_blocking rejects invalid windows and maps", {
  map <- make_window_map_fixture()

  expect_error(perform_window_blocking(map, window = 0,    method = "window_snp"), "positive")
  expect_error(perform_window_blocking(map, window = -5,   method = "window_map"), "positive")
  expect_error(perform_window_blocking(map, window = c(1, 2), method = "window_snp"), "positive")
  expect_error(perform_window_blocking(map, window = "10", method = "window_snp"), "positive")

  # a fractional window is meaningless as a marker count, but fine as a distance
  expect_error(perform_window_blocking(map, window = 2.5, method = "window_snp"), "whole number")
  expect_silent(perform_window_blocking(map, window = 2.5, method = "window_map"))

  expect_error(perform_window_blocking(map, window = 3, method = "window_kb"), "arg")

  expect_error(perform_window_blocking(map[, c("SNP", "Position")], window = 3), "SNP")
  expect_error(perform_window_blocking("not a map", window = 3), "data frame")
})


# Tests: downstream compatibility ---------------------------------------------

test_that("window blocks feed into block_obj_to_df like LD blocks do", {
  map <- make_window_map_fixture()

  blocks   <- perform_window_blocking(map, window = 200, method = "window_map")
  block_df <- block_obj_to_df(blocks, map)

  expect_s3_class(block_df, "data.frame")
  expect_equal(nrow(block_df), 7)   # 5 blocks on chr 1 + 2 on chr 2
  expect_true(all(c("Block", "Block_ID", "Num_SNP", "First_SNP", "Last_SNP",
                    "Chrom", "Start_Pos", "End_Pos", "Block_Size") %in% colnames(block_df)))

  expect_equal(sum(block_df$Num_SNP), nrow(map))

  # the [200,400) block on chr 1 spans c2..c3
  two_snp <- block_df[block_df$First_SNP == "c2", ]
  expect_equal(two_snp$Last_SNP, "c3")
  expect_equal(two_snp$Start_Pos, 250)
  expect_equal(two_snp$End_Pos, 300)
  expect_equal(two_snp$Block_Size, 50)
})
