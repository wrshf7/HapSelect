# Tests: parse_chromosome ------------------------------------------------------
#
# Tests that parse_chromosome() correctly parses chromosome labels to numeric values, and that it issues warnings when appropriate.

pc <- function(x) suppressWarnings(parse_chromosome(x))


test_that("numeric chromosomes are returned untouched and without warning", {
  expect_equal(parse_chromosome(c(1, 2, 10)), c(1, 2, 10))
  expect_silent(parse_chromosome(c(1, 2, 10)))
  expect_equal(parse_chromosome(c(3L, 1L, 2L)), c(3L, 1L, 2L))
})


test_that("numeric-looking chromosomes keep their own number", {
  expect_equal(pc(c("1", "2", "10")), c(1, 2, 10))
  expect_equal(pc(as.character(1:20)), as.numeric(1:20))
  expect_equal(pc(factor(c("1", "2", "10"))), c(1, 2, 10))
  expect_equal(pc(c(" 1", "2 ", " 10 ")), c(1, 2, 10))
})


test_that("a leading non-numeric prefix is stripped", {
  expect_equal(pc(c("chr1", "chr2", "chr10")), c(1, 2, 10))
  expect_equal(pc(c("ch9", "Chr02", "CHR11")), c(9, 2, 11))
  expect_equal(pc(c("LG1", "LG2", "LG12")),    c(1, 2, 12))

  # prefixed and bare forms of the same chromosome resolve to the same number
  expect_equal(pc(c("chr1", "2", "chr10")), c(1, 2, 10))
  expect_equal(pc(c("chr1", "1", "Chr1")),  c(1, 1, 1))
})


test_that("labels with no number of their own are assigned one above the real chromosomes", {
  expect_equal(pc(c("1", "2", "10", "X", "Y")), c(1, 2, 10, 11, 12))
  expect_equal(pc(c("chr1", "chr2", "chrX", "chrMT")), c(1, 2, 4, 3))
  expect_equal(pc(c("X", "Y", "MT")), c(2, 3, 1))
  expect_equal(pc(c("1A", "1B", "2A", "2D")), c(1, 2, 3, 4))
})


test_that("assigned numbers cannot collide with real chromosomes and are deterministic", {
  x <- c("7", "3", "22", "X", "Y", "MT")
  got <- pc(x)

  real      <- got[x %in% c("7", "3", "22")]
  synthetic <- got[x %in% c("X", "Y", "MT")]

  expect_equal(real, c(7, 3, 22))
  expect_true(all(synthetic > max(real)))
  expect_equal(length(unique(got)), length(unique(x)))

  # label order, not order of appearance, so repeated runs agree
  expect_equal(pc(c("X", "Y", "MT", "1")), pc(c("1", "MT", "Y", "X"))[c(4, 3, 2, 1)])
})


test_that("missing chromosomes stay missing", {
  expect_equal(pc(c("1", "2", NA, "", "10")), c(1, 2, NA, NA, 10))
  expect_true(is.na(pc(c(NA_character_))))

  # an NA must not be treated as a label needing a number
  expect_equal(pc(c("1", NA, "X")), c(1, NA, 2))
})


test_that("parse_chromosome warns, and names the labels it numbered", {
  expect_warning(parse_chromosome(c("chr1", "chr2")), "not numeric")
  expect_warning(parse_chromosome(c("1", "X", "Y")), "X -> 2, Y -> 3")
  expect_warning(parse_chromosome(c("1", "2")), "not numeric")

  # nothing was assigned, so no assignment list in the message
  w <- tryCatch(parse_chromosome(c("1", "2")), warning = conditionMessage)
  expect_false(grepl("assigned one", w))
})


# Tests: check_file ------------------------------------------------------------

test_that("check_file reads factor positions as values, not level codes", {
  map <- data.frame(SNP = c("s1", "s2"), Chromosome = c(1, 1),
                    Position = factor(c("100", "200")))

  expect_equal(suppressWarnings(check_file(map))[, 3], c(100, 200))
  expect_warning(check_file(map), "Positions were not numeric")
})


test_that("check_file coerces factor SNP IDs to character", {
  map <- data.frame(SNP = factor(c("s1", "s2")), Chromosome = c(1, 1),
                    Position = c(100, 200))

  expect_type(suppressWarnings(check_file(map))[, 1], "character")
  expect_equal(suppressWarnings(check_file(map))[, 1], c("s1", "s2"))
  expect_warning(check_file(map), "SNP ID were not characters")
})


test_that("check_file leaves a well-formed map alone", {
  map <- data.frame(SNP = c("s1", "s2"), Chromosome = c(1, 10),
                    Position = c(100, 200), stringsAsFactors = FALSE)

  expect_silent(check_file(map))
  expect_equal(check_file(map), map)
})


test_that("check_file still rejects a map that is not a 3-column data frame", {
  expect_error(check_file("not a map"), "at least 3 columns")
  expect_error(check_file(data.frame(SNP = "s1", Chromosome = 1)), "at least 3 columns")
})


# Tests: order_map -------------------------------------------------------------

test_that("order_map preserves chromosome identity end to end", {
  map <- data.frame(SNP        = c("a", "b", "c", "d", "e", "f"),
                    Chromosome = c("chr1", "chr1", "chr10", "chr10", "chrX", "chrX"),
                    Position   = c(20, 10, 20, 10, 20, 10),
                    stringsAsFactors = FALSE)

  out <- suppressWarnings(order_map(map))

  # chromosome 10 stays 10; chrX becomes 11, above the real chromosomes
  expect_equal(out$Chromosome[out$SNP %in% c("c", "d")], c(10, 10))
  expect_equal(out$Chromosome[out$SNP %in% c("e", "f")], c(11, 11))

  # and positions are sorted within each chromosome
  expect_equal(out$Position, c(10, 20, 10, 20, 10, 20))
  expect_equal(out$SNP,      c("b", "a", "d", "c", "f", "e"))
})
