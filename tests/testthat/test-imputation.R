# Tests: Beagle imputation and phasing wrappers -------------------------------
#
# Covers the geno-level wrappers in R/imputation.R. Only their input validation
# is exercised here, since running the real thing needs a Beagle install. The
# VCF reading and writing they delegate to is covered in test-vcf.R.

test_that("beagle_impute errors when the input VCF is missing", {
  expect_error(
    beagle_impute("does_not_exist.vcf", tempfile()),
    "Beagle input VCF not found"
  )
})

test_that("beagle_impute_geno errors when geno is not a valid data frame", {
  expect_error(
    beagle_impute_geno("not_a_data_frame"),
    "geno must be a data frame"
  )
  expect_error(
    beagle_impute_geno(data.frame(a = 1, b = 2, c = 3)),
    "geno must be a data frame"
  )
})

test_that("beagle_phase_geno errors when geno is not a valid data frame", {
  expect_error(
    beagle_phase_geno("not_a_data_frame"),
    "geno must be a data frame"
  )
  expect_error(
    beagle_phase_geno(data.frame(a = 1, b = 2, c = 3)),
    "geno must be a data frame"
  )
})
