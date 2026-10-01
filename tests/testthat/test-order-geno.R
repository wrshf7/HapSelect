# Tests: order_geno ------------------------------------------------------------
#
# order_geno() is order_map() for a table that carries dosages alongside its map
# columns. The contract it has to keep is that geno and map end up describing the
# same markers in the same order, because blocking reads a marker's index from
# the map and its dosages from the matching row of geno.

geno_fixture <- function(chrom = c("chr1A", "chr1A", "chr1B"),
                         snp   = c("snpB", "snpA", "snpC"),
                         pos   = c(300, 100, 200)) {
  data.frame(
    SNP        = snp,
    Chromosome = chrom,
    Position   = pos,
    Ind1       = c(0, 2, 1),
    Ind2       = c(1, 1, 2),
    stringsAsFactors = FALSE
  )
}


# Ordering by its own map columns ----------------------------------------------

test_that("order_geno sorts by chromosome then position when no map is given", {
  ordered <- order_geno(geno_fixture())

  expect_equal(ordered$SNP, c("snpA", "snpB", "snpC"))
  expect_equal(ordered$Position, c(100, 300, 200))
  expect_equal(ordered$Chromosome, c(1, 1, 2))
})


test_that("the first three columns are exactly what order_map would give", {
  # the guarantee that keeps the two from drifting apart
  geno <- geno_fixture()

  expect_equal(order_geno(geno)[, 1:3],
               order_map(geno[, 1:3], verbose = FALSE))
})


test_that("dosages travel with their marker, not with their row number", {
  geno    <- geno_fixture()
  ordered <- order_geno(geno)

  dosages <- function(d, snp) unlist(d[d$SNP == snp, c("Ind1", "Ind2")], use.names = FALSE)

  for (snp in geno$SNP) {
    expect_equal(dosages(ordered, snp), dosages(geno, snp), info = snp)
  }
})


test_that("the result does not depend on the order the rows arrived in", {
  geno <- geno_fixture()

  expect_equal(order_geno(geno[c(3, 1, 2), ]), order_geno(geno))
  expect_equal(order_geno(geno[c(2, 3, 1), ]), order_geno(geno))
})


test_that("an already-ordered genotype table is left alone", {
  once  <- order_geno(geno_fixture())
  twice <- order_geno(once)

  expect_equal(twice, once)
})


# Ordering to follow a supplied map ---------------------------------------------

test_that("order_geno puts geno into the map's row order", {
  geno <- geno_fixture()
  map  <- order_map(geno[, 1:3], verbose = FALSE)

  aligned <- order_geno(geno, map)

  expect_equal(aligned$SNP, map$SNP)
  expect_equal(aligned$Position, map$Position)
  expect_equal(aligned$Chromosome, map$Chromosome)
})


test_that("the map is the authority for chromosome and position", {
  # geno keeps the VCF's text labels; the map carries the numbering, and that is
  # what the result has to end up with
  geno <- geno_fixture()
  map  <- order_map(geno[, 1:3], verbose = FALSE)

  expect_type(geno$Chromosome, "character")
  expect_true(is.numeric(order_geno(geno, map)$Chromosome))
})


test_that("the aligned genotypes satisfy check_ld_matrix", {
  # the whole point: a VCF read cannot be passed to the LD or blocking functions
  # until it has been through here
  geno <- geno_fixture()
  map  <- order_map(geno[, 1:3], verbose = FALSE)

  expect_error(check_ld_matrix(geno), "numeric chromosome")
  expect_silent(check_ld_matrix(order_geno(geno, map)))
})


test_that("dosages follow the marker when the map reorders it", {
  geno <- geno_fixture()
  map  <- order_map(geno[, 1:3], verbose = FALSE)

  aligned <- order_geno(geno, map)

  # snpA was row 2 of geno and is row 1 of the map
  expect_equal(aligned$Ind1[aligned$SNP == "snpA"], geno$Ind1[geno$SNP == "snpA"])
  expect_equal(aligned$Ind2[aligned$SNP == "snpC"], geno$Ind2[geno$SNP == "snpC"])
})


# What it refuses ----------------------------------------------------------------

test_that("a genotype table with no dosage columns is rejected", {
  expect_error(order_geno(geno_fixture()[, 1:3]), "at least 4 columns")
  expect_error(order_geno("not a data frame"),    "at least 4 columns")
})


test_that("duplicated SNP IDs are refused on either side", {
  dup_geno <- geno_fixture(snp = c("snpA", "snpA", "snpC"))
  expect_error(order_geno(dup_geno), "geno has duplicated SNP IDs")

  geno <- geno_fixture()
  map  <- order_map(geno[, 1:3], verbose = FALSE)
  map$SNP[2] <- map$SNP[1]
  expect_error(order_geno(geno, map), "map has duplicated SNP IDs")
})


test_that("a marker on only one side is refused, and named", {
  geno <- geno_fixture()
  map  <- order_map(geno[, 1:3], verbose = FALSE)

  expect_error(order_geno(geno, map[-1, ]),       "same markers")
  expect_error(order_geno(geno[-1, ], map),       "same markers")
  expect_error(order_geno(geno[-1, ], map),       "snpB")
})


test_that("a map whose chromosomes were never numbered is refused", {
  geno <- geno_fixture()
  raw  <- geno[, 1:3]   # still carries "chr1A" labels

  expect_error(order_geno(geno, raw), "order_map")
})
