# Fixtures --------------------------------------------------------------------
#
# Two-chromosome dataset designed so expected block membership can be
# worked out by hand:
#
#   Chr 1 (5 markers):  m1-m2 (0.90), m2-m3 (0.85) are high-LD neighbours;
#                       m3-m4 (0.10) and m4-m5 (0.10) fall below threshold.
#                       Expected blocks: [m1,m2,m3]  [m4]  [m5]
#
#   Chr 2 (4 markers):  s1-s2 (0.90) and s3-s4 (0.80) are tight pairs
#                       separated by a low-LD gap s2-s3 (0.10).
#                       Expected blocks: [s1,s2]  [s3,s4]

make_ld_fixture <- function() {
  data.frame(
    Chrom  = c(rep(1L, 10), rep(2L, 6)),
    Locus1 = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4,
               1, 1, 1, 2, 2, 3),
    Locus2 = c(2, 3, 4, 5, 3, 4, 5, 4, 5, 5,
               2, 3, 4, 3, 4, 4),
    Name1  = c("m1","m1","m1","m1","m2","m2","m2","m3","m3","m4",
               "s1","s1","s1","s2","s2","s3"),
    Name2  = c("m2","m3","m4","m5","m3","m4","m5","m4","m5","m5",
               "s2","s3","s4","s3","s4","s4"),
    LD     = c(0.90, 0.80, 0.10, 0.05,
               0.85, 0.10, 0.05,
               0.10, 0.05, 0.10,
               0.90, 0.10, 0.05,
               0.10, 0.05,
               0.80),
    stringsAsFactors = FALSE
  )
}

make_map_fixture <- function() {
  data.frame(
    SNP      = c("m1","m2","m3","m4","m5","s1","s2","s3","s4"),
    Chrom    = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
    Position = c(100, 200, 300, 400, 500, 100, 200, 300, 400),
    stringsAsFactors = FALSE
  )
}

# Tolerance fixture -----------------------------------------------------------
#
# Single chromosome, 5 markers. A low-LD bridge marker (t4) sits between two
# high-LD clusters. t3 and t5 share high LD (0.80), so with tolerance >= 1
# the bridge is absorbed and all five markers form one block.
#
#   tolerance = 0 → [t1,t2,t3]  [t4,t5]
#   tolerance = 1 → [t1,t2,t3,t4,t5]

make_tol_ld_fixture <- function() {
  data.frame(
    Chrom  = rep(1L, 10),
    Locus1 = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4),
    Locus2 = c(2, 3, 4, 5, 3, 4, 5, 4, 5, 5),
    Name1  = c("t1","t1","t1","t1","t2","t2","t2","t3","t3","t4"),
    Name2  = c("t2","t3","t4","t5","t3","t4","t5","t4","t5","t5"),
    LD     = c(0.90, 0.85, 0.10, 0.05,
               0.85, 0.10, 0.05,
               0.10, 0.80,
               0.90),
    stringsAsFactors = FALSE
  )
}

make_tol_map_fixture <- function() {
  data.frame(
    SNP      = c("t1","t2","t3","t4","t5"),
    Chrom    = rep(1L, 5),
    Position = c(100, 200, 300, 400, 500),
    stringsAsFactors = FALSE
  )
}


# Tests: def_blocks -----------------------------------------------------------

test_that("def_blocks returns a named list with one entry per chromosome", {
  ld  <- make_ld_fixture()
  map <- make_map_fixture()

  blocks <- def_blocks(ld_strategy(ld,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 0,
                       tol_reset = FALSE,
                       start     = "LD",
                       parallel  = FALSE), map)

  expect_type(blocks, "list")
  expect_named(blocks, c("1", "2"))
  expect_length(blocks[["1"]], 3)   # [m1,m2,m3], [m4], [m5]
  expect_length(blocks[["2"]], 2)   # [s1,s2], [s3,s4]
})


test_that("def_blocks assigns correct markers to blocks (start = 'LD', flanking, tolerance = 0)", {
  ld  <- make_ld_fixture()
  map <- make_map_fixture()

  blocks <- def_blocks(ld_strategy(ld,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 0,
                       tol_reset = FALSE,
                       start     = "LD",
                       parallel  = FALSE), map)

  chr1 <- blocks[["1"]]
  chr2 <- blocks[["2"]]

  # Chr 1: high-LD run m1-m3 forms one block; m4 and m5 are singletons
  expect_equal(chr1[[1]], c("m1", "m2", "m3"))
  expect_equal(chr1[[2]], "m4")
  expect_equal(chr1[[3]], "m5")

  # Chr 2: two separated tight pairs
  expect_equal(chr2[[1]], c("s1", "s2"))
  expect_equal(chr2[[2]], c("s3", "s4"))
})


test_that("def_blocks produces the same blocks with start = 'beginning' on this fixture", {
  ld  <- make_ld_fixture()
  map <- make_map_fixture()

  blocks <- def_blocks(ld_strategy(ld,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 0,
                       tol_reset = FALSE,
                       start     = "beginning",
                       parallel  = FALSE), map)

  chr1 <- blocks[["1"]]
  chr2 <- blocks[["2"]]

  expect_equal(chr1[[1]], c("m1", "m2", "m3"))
  expect_equal(chr1[[2]], "m4")
  expect_equal(chr1[[3]], "m5")

  expect_equal(chr2[[1]], c("s1", "s2"))
  expect_equal(chr2[[2]], c("s3", "s4"))
})


test_that("tolerance = 0 stops extension at a low-LD bridge marker", {
  ld  <- make_tol_ld_fixture()
  map <- make_tol_map_fixture()

  blocks <- def_blocks(ld_strategy(ld,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 0,
                       tol_reset = FALSE,
                       start     = "LD",
                       parallel  = FALSE), map)

  chr1 <- blocks[["1"]]
  expect_length(chr1, 2)
  expect_equal(chr1[[1]], c("t1", "t2", "t3"))
  expect_equal(chr1[[2]], c("t4", "t5"))
})


test_that("tolerance = 1 absorbs a low-LD bridge marker into a single block", {
  ld  <- make_tol_ld_fixture()
  map <- make_tol_map_fixture()

  blocks <- def_blocks(ld_strategy(ld,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 1,
                       tol_reset = TRUE,
                       start     = "LD",
                       parallel  = FALSE), map)

  chr1 <- blocks[["1"]]
  expect_length(chr1, 1)
  expect_equal(chr1[[1]], c("t1", "t2", "t3", "t4", "t5"))
})


# Tests: the map is the marker universe ----------------------------------------
#
# An LD table is a table of pairs, so a marker with no qualifying pair does not
# appear in it at all. That is the normal result of an r^2 floor - pairwise_ld()
# takes min_r2, PLINK takes --ld-window-r2, and graph_strategy() defaults to
# ld_min_r2 = 0.20 - so a filtered table is ordinary input, not a malformed one.
#
# The marker list therefore has to come from the map, which is complete, rather
# than from the LD table, which is not. A marker the LD table never mentions is
# in LD with nothing above the floor, which makes it a single-marker block, not
# a marker to drop on the floor.

# Same markers as make_map_fixture(), with every pair involving m5 or s4 removed,
# as an r^2 floor would remove them.
make_filtered_ld_fixture <- function() {
  ld <- make_ld_fixture()
  ld[!(ld$Name1 %in% c("m5", "s4") | ld$Name2 %in% c("m5", "s4")), ]
}

test_that("a marker missing from the LD table is still blocked, as a singleton", {
  ld  <- make_filtered_ld_fixture()
  map <- make_map_fixture()

  # the fixture really does hide them: they are in the map and not in the LD table
  expect_false(any(c("m5", "s4") %in% c(ld$Name1, ld$Name2)))
  expect_true(all(c("m5", "s4") %in% map$SNP))

  blocks <- def_blocks(ld_strategy(ld, method = "flanking", threshold = 0.7,
                                   tolerance = 0, tol_reset = FALSE,
                                   start = "LD", parallel = FALSE), map)

  placed <- unlist(blocks, use.names = FALSE)

  # nothing in the map may go missing, whatever the LD table happens to contain
  expect_setequal(placed, map$SNP)
  expect_length(placed, nrow(map))

  # and the two hidden markers stand alone, since nothing links them to anything
  expect_true(list("m5") %in% blocks[["1"]])
  expect_true(list("s4") %in% blocks[["2"]])
})


test_that("markers dropped by an r-squared floor survive the round trip", {
  # The realistic route in: compute LD with a floor, then block with it. The
  # floor is what makes markers vanish from the table, so this is the path that
  # loses them.
  set.seed(11)
  n_ind <- 40
  base  <- sample(c(0, 1, 2), n_ind, TRUE)
  near  <- function(k) { x <- base; i <- sample(n_ind, k); x[i] <- sample(c(0,1,2), k, TRUE); x }

  geno <- data.frame(
    SNP        = c("a1", "a2", "a3", "a4"),
    Chromosome = 1,
    Position   = c(100, 200, 300, 400),
    rbind(near(2), near(2), near(2), sample(c(0, 1, 2), n_ind, TRUE)),
    stringsAsFactors = FALSE
  )
  names(geno)[-(1:3)] <- paste0("I", seq_len(n_ind))
  map <- geno[, 1:3]

  ld <- pairwise_ld(geno, parallelize = FALSE, min_r2 = 0.5)

  # a4 is independent of the rest, so the floor removes every pair it is in
  expect_false("a4" %in% c(ld$Name1, ld$Name2))

  blocks <- suppressMessages(
    def_blocks(ld_strategy(ld, method = "flanking", threshold = 0.7,
                           parallel = FALSE), map)
  )

  expect_setequal(unlist(blocks, use.names = FALSE), map$SNP)
})


test_that("filtering the LD table below the blocking threshold changes nothing", {
  # The sharpest statement of what the marker list being taken from the map buys:
  # pairs under the blocking threshold can never seed or extend a block, so
  # removing them is a no-op. Before the marker list came from the map, the same
  # filtering silently cost every marker whose pairs were all removed.
  #
  # This holds for "flanking", which tests one pair at a time. It does not hold
  # for "average", which averages LD over a block's members, so a removed
  # below-threshold pair still moves the mean.
  ld  <- make_ld_fixture()
  map <- make_map_fixture()

  blocking <- function(ld) {
    def_blocks(ld_strategy(ld, method = "flanking", threshold = 0.7, tolerance = 0,
                           tol_reset = FALSE, start = "LD", parallel = FALSE), map)
  }

  unfiltered <- blocking(ld)
  filtered   <- blocking(ld[ld$LD >= 0.2, ])

  # the filter really does remove pairs, and with them a marker's only mentions
  expect_lt(nrow(ld[ld$LD >= 0.2, ]), nrow(ld))
  expect_identical(filtered, unfiltered)
})
