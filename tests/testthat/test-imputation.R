test_that("write_vcf_geno and read_vcf_geno round-trip a genotype data frame", {
  geno <- data.frame(
    Marker = c("snp1", "snp2", "snp3"),
    Chrom = c(1, 1, 2),
    Position = c(100, 200, 100),
    Ind1 = c(0, 1, NA),
    Ind2 = c(2, NA, 1),
    stringsAsFactors = FALSE
  )

  vcf_path <- tempfile(fileext = ".vcf")
  on.exit(unlink(vcf_path))

  HapSelect:::write_vcf_geno(geno, vcf_path)
  observed <- read_vcf_geno(vcf_path)

  expected <- data.frame(
    Marker = geno$Marker,
    Chrom = geno$Chrom,
    Position = geno$Position,
    Ind1 = geno$Ind1,
    Ind2 = geno$Ind2,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_equal(observed, expected)
})

test_that("read_vcf_geno errors when the file is missing or malformed", {
  expect_error(read_vcf_geno("does_not_exist.vcf"), "VCF file not found")

  bad_vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(bad_vcf))
  writeLines(c("##fileformat=VCFv4.2", "1\t100\tsnp1\tA\tG\t.\tPASS\t.\tGT\t0/0"), bad_vcf)

  expect_error(read_vcf_geno(bad_vcf), "#CHROM header line")
})

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
