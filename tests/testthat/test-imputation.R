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

test_that("the Beagle wrappers reject non-diploid dosages before running Beagle", {
  # A dosage of 3 is valid for a polyploid, but Beagle only handles diploids
  geno <- data.frame(SNP = "snp1", Chromosome = 1, Position = 100, Ind1 = 3, stringsAsFactors = FALSE)

  expect_error(beagle_impute_geno(geno), "Beagle only supports diploid genotypes")
  expect_error(beagle_phase_geno(geno), "Beagle only supports diploid genotypes")
})

# prefer_ds --------------------------------------------------------------------
# Beagle writes a DS field only when it imputes against a reference panel. These
# tests pin down that beagle_impute_geno() asks the reader for the field the
# caller wanted, and says so when the run cannot produce one.

test_that("beagle_impute_geno asks the reader for the field prefer_ds names", {
  geno <- data.frame(SNP = "snp1", Chromosome = 1, Position = 100, Ind1 = 1,
                     stringsAsFactors = FALSE)
  asked <- NULL

  # Capture the reader run_beagle_on_geno() would have been given, rather than
  # running Beagle, then ask it what it would read.
  local_mocked_bindings(
    run_beagle_on_geno = function(geno, reader, ...) {
      asked <<- reader
      invisible(NULL)
    },
    read_vcf_geno = function(path, prefer_ds = FALSE, ...) prefer_ds
  )

  beagle_impute_geno(geno, ref = "panel.vcf", prefer_ds = TRUE)
  expect_true(asked("ignored.vcf"))

  beagle_impute_geno(geno, ref = "panel.vcf", prefer_ds = FALSE)
  expect_false(asked("ignored.vcf"))
})

test_that("beagle_impute_geno defaults to hard calls, so its dosages stay whole", {
  geno <- data.frame(SNP = "snp1", Chromosome = 1, Position = 100, Ind1 = 1,
                     stringsAsFactors = FALSE)
  asked <- NULL

  local_mocked_bindings(
    run_beagle_on_geno = function(geno, reader, ...) {
      asked <<- reader
      invisible(NULL)
    },
    read_vcf_geno = function(path, prefer_ds = FALSE, ...) prefer_ds
  )

  beagle_impute_geno(geno, ref = "panel.vcf")
  expect_false(asked("ignored.vcf"))
})

test_that("beagle_impute_geno warns that prefer_ds needs a reference panel", {
  geno <- data.frame(SNP = "snp1", Chromosome = 1, Position = 100, Ind1 = 1,
                     stringsAsFactors = FALSE)

  local_mocked_bindings(run_beagle_on_geno = function(...) invisible(NULL))

  expect_warning(
    beagle_impute_geno(geno, prefer_ds = TRUE),
    "no effect without a reference panel"
  )
  expect_silent(beagle_impute_geno(geno, prefer_ds = FALSE))
  expect_silent(beagle_impute_geno(geno, ref = "panel.vcf", prefer_ds = TRUE))
})
