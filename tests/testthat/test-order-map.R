# Tests: number_chromosomes ----------------------------------------------------
#

# Run without warnings to avoid cluttering the test output with expected warnings.
number_chromosomes_suppressed <- function(x) suppressWarnings(number_chromosomes(x))

# Core contract ----------------------------------------------------------------

test_that("a numeric column is returned untouched and in silence", {
  expect_equal(number_chromosomes(c(1, 2, 10)), c(1, 2, 10))
  expect_equal(number_chromosomes(c(3L, 1L, 2L)), c(3L, 1L, 2L))
  expect_silent(number_chromosomes(c(1, 2, 10)))
})


test_that("distinct labels always get distinct integers", {
  # nothing is inferred to be the same chromosome as anything else
  for (x in list(c("chr1", "1", "Chr1"),
                 c("lg1", "chrom1", "1"),
                 c("chr1", "Chr1", "CHR1"),
                 c("11", "12", "13", "scaffold_12"),
                 c("8", "9", "-9"),
                 c("1", "1A", "1B", "1D"),
                 c("chr1", "chr1_random", "chr1_alt"),
                 c(as.character(1:10), "scaffold_3", "-9", "X"))) {
    got <- number_chromosomes_suppressed(x)
    expect_equal(length(unique(got)), length(unique(x)))
    expect_false(anyNA(got))
  }
})


test_that("the assigned integers are whole numbers", {
  got <- number_chromosomes_suppressed(c("chr1", "X", "scaffold_9"))
  expect_true(all(got == as.integer(got)))
  expect_true(all(got > 0))
})


test_that("missing chromosomes stay missing", {
  expect_equal(number_chromosomes_suppressed(c("1", "2", NA, "", "10"))[c(3, 4)], c(NA_real_, NA_real_))
  expect_true(is.na(number_chromosomes_suppressed(NA_character_)))

  # an NA must not be handed an integer of its own
  got <- number_chromosomes_suppressed(c("A", NA, "B"))
  expect_true(is.na(got[2]))
  expect_equal(length(unique(got[!is.na(got)])), 2)
})


test_that("factors and whitespace are handled before mapping", {
  expect_equal(number_chromosomes_suppressed(factor(c("1", "2", "10"))), number_chromosomes_suppressed(c("1", "2", "10")))
  expect_equal(number_chromosomes_suppressed(c(" 1", "2 ", " 10 ")), number_chromosomes_suppressed(c("1", "2", "10")))

  # padding is not a difference, so these are one chromosome, not two
  expect_equal(length(unique(number_chromosomes_suppressed(c("chr1", " chr1 ")))), 1)
})


# Stability --------------------------------------------------------------------

test_that("the mapping does not depend on row order", {
  x <- c("1", "X", "Y", "10")
  expect_equal(number_chromosomes_suppressed(x), rev(number_chromosomes_suppressed(rev(x))))
})


test_that("the mapping does not depend on how often a label appears", {
  expect_equal(unique(number_chromosomes_suppressed(c("X", "1", "1", "1"))), unique(number_chromosomes_suppressed(c("X", "1"))))
})


test_that("mapping an already-mapped column changes nothing", {
  once <- number_chromosomes_suppressed(c("chr1", "X", "chr10"))
  expect_silent(twice <- number_chromosomes(once))
  expect_equal(twice, once)
})


# Numbering follows the labels --------------------------------------------------

test_that("labels are numbered 1 to n in natural order", {
  expect_equal(number_chromosomes_suppressed(c("1", "2", "10")), c(1, 2, 3))
  expect_equal(number_chromosomes_suppressed(as.character(1:20)), as.numeric(1:20))
  expect_equal(number_chromosomes_suppressed(c("chr1", "chr2", "chr10")), c(1, 2, 3))

  # leading zeros make a different label, so these are three chromosomes
  expect_equal(length(unique(number_chromosomes_suppressed(c("01", "1", "001")))), 3)
})


test_that("the numbers always run in the same order as the labels", {
  for (x in list(c("7", "3", "22", "X", "Y", "MT"),
                 c("1", "2_1", "2_2", "3"),
                 c("1.1", "1.2", "2", "10"),
                 c("chr1", "chr2", "chr10", "chrX"),
                 c("11", "12", "13", "scaffold_12"))) {
    got <- number_chromosomes_suppressed(x)

    # a dense 1 to n with no gaps, whatever the labels looked like
    expect_equal(sort(unique(got)), as.numeric(seq_along(unique(x))))

    # putting the rows in label order puts the numbers in ascending order too
    expect_false(is.unsorted(got[order(natural_key(x))]))
  }
})


test_that("named chromosomes are numbered after the numbered ones", {
  x   <- c("7", "3", "22", "X", "Y", "MT")
  got <- number_chromosomes_suppressed(x)

  # digits sort before letters, so 3 < 7 < 22 are numbered first, then MT < X < Y
  expect_equal(got[x %in% c("7", "3", "22")], c(2, 1, 3))
  expect_true(all(got[x %in% c("X", "Y", "MT")] > 3))
})


test_that("a name containing a number is still its own chromosome", {
  got <- number_chromosomes_suppressed(c("11", "12", "13", "scaffold_12"))
  expect_equal(length(unique(got)), 4)
  expect_false(got[4] == got[2])          # scaffold_12 is not chromosome 12

  got <- number_chromosomes_suppressed(c("8", "9", "-9"))
  expect_equal(length(unique(got)), 3)
  expect_false(got[3] == got[2])          # "-9" is not chromosome 9
})


# Awkward labels ---------------------------------------------------------------

test_that("sub-chromosome labels keep their place among the chromosomes", {
  expect_equal(number_chromosomes_suppressed(c("1", "2_1", "2_2", "3")), c(1, 2, 3, 4))
})


test_that("decimal labels keep their place among the chromosomes", {
  x <- c("1.1", "1.2", "2", "3", "4", "5", "6", "7", "8")
  expect_equal(number_chromosomes_suppressed(x), as.numeric(1:9))
})


test_that("a sub-chromosome label colliding with a real chromosome is flagged", {
  # "2_1" and "21" are the same once punctuation is stripped
  x <- c("1", "2_1", "2_2", "3", "21")

  expect_warning(number_chromosomes(x), "2_1")
  expect_equal(number_chromosomes_suppressed(x), c(1, 2, 3, 4, 5))
})


# Ordering ---------------------------------------------------------------------

test_that("labels are ordered with their digits compared as numbers", {
  # chr10 must not sort between chr1 and chr2
  got <- number_chromosomes_suppressed(c("chr1", "chr2", "chr10"))
  expect_equal(got, sort(got))
  expect_true(got[2] < got[3])

  got <- number_chromosomes_suppressed(c("scaffold_1", "scaffold_2", "scaffold_10", "scaffold_20"))
  expect_equal(got, sort(got))
})


test_that("a single-scheme column keeps the order its names imply", {
  for (x in list(paste0("chr", 1:12), paste0("LG", 1:12), paste0("scaffold_", 1:12))) {
    expect_equal(number_chromosomes_suppressed(x), as.numeric(1:12))
  }
})


# Reporting --------------------------------------------------------------------

test_that("verbose reports the label to integer mapping", {
  expect_message(number_chromosomes(c("chr1", "chr2", "X"), verbose = TRUE), "chr1")
  expect_message(number_chromosomes(c("chr1", "chr2", "X"), verbose = TRUE), "X")

  # nothing to report for a column that was already numeric
  expect_silent(number_chromosomes(c(1, 2, 10), verbose = TRUE))
})


test_that("verbose can be switched off", {
  expect_silent(suppressWarnings(number_chromosomes(c("chr1", "chr2"), verbose = FALSE)))
})


test_that("labels differing only in case or spacing are flagged", {
  expect_warning(number_chromosomes(c("chr1", "Chr1")), "chr1")
  expect_warning(number_chromosomes(c("chr 1", "chr1")), "chr1")

  # flagged, but still numbered separately
  expect_equal(length(unique(number_chromosomes_suppressed(c("chr1", "Chr1", "CHR1")))), 3)
})


test_that("labels that look like missing data are flagged", {
  for (sentinel in c("NA", "N/A", ".", "-", "?", "-9")) {
    expect_warning(number_chromosomes(c("1", "2", sentinel)), "missing")
  }
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

test_that("order_map sorts within chromosome and keeps chromosomes distinct", {
  map <- data.frame(SNP        = c("a", "b", "c", "d", "e", "f"),
                    Chromosome = c("chr1", "chr1", "chr2", "chr2", "chrX", "chrX"),
                    Position   = c(20, 10, 20, 10, 20, 10),
                    stringsAsFactors = FALSE)

  out <- suppressWarnings(order_map(map))

  expect_equal(length(unique(out$Chromosome)), 3)
  expect_equal(out$Position, c(10, 20, 10, 20, 10, 20))
  expect_equal(out$SNP,      c("b", "a", "d", "c", "f", "e"))
})
