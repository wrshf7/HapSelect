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

  blocks <- def_blocks(ld, map,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 0,
                       tol_reset = FALSE,
                       start     = "LD",
                       parallel  = FALSE)

  expect_type(blocks, "list")
  expect_named(blocks, c("1", "2"))
  expect_length(blocks[["1"]], 3)   # [m1,m2,m3], [m4], [m5]
  expect_length(blocks[["2"]], 2)   # [s1,s2], [s3,s4]
})


test_that("def_blocks assigns correct markers to blocks (start = 'LD', flanking, tolerance = 0)", {
  ld  <- make_ld_fixture()
  map <- make_map_fixture()

  blocks <- def_blocks(ld, map,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 0,
                       tol_reset = FALSE,
                       start     = "LD",
                       parallel  = FALSE)

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

  blocks <- def_blocks(ld, map,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 0,
                       tol_reset = FALSE,
                       start     = "beginning",
                       parallel  = FALSE)

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

  blocks <- def_blocks(ld, map,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 0,
                       tol_reset = FALSE,
                       start     = "LD",
                       parallel  = FALSE)

  chr1 <- blocks[["1"]]
  expect_length(chr1, 2)
  expect_equal(chr1[[1]], c("t1", "t2", "t3"))
  expect_equal(chr1[[2]], c("t4", "t5"))
})


test_that("tolerance = 1 absorbs a low-LD bridge marker into a single block", {
  ld  <- make_tol_ld_fixture()
  map <- make_tol_map_fixture()

  blocks <- def_blocks(ld, map,
                       method    = "flanking",
                       threshold = 0.7,
                       tolerance = 1,
                       tol_reset = TRUE,
                       start     = "LD",
                       parallel  = FALSE)

  chr1 <- blocks[["1"]]
  expect_length(chr1, 1)
  expect_equal(chr1[[1]], c("t1", "t2", "t3", "t4", "t5"))
})
