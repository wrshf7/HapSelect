test_that("write_vcf_geno and read_vcf_geno round-trip a genotype data frame", {
  geno <- data.frame(
    SNP = c("snp1", "snp2", "snp3"),
    Chromosome = c(1, 1, 2),
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
    SNP = geno$SNP,
    Chromosome = geno$Chromosome,
    Position = geno$Position,
    Ind1 = geno$Ind1,
    Ind2 = geno$Ind2,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_equal(observed, expected)

  # Pin the geno column convention, not just round-trip self-consistency
  expect_true(all(c("SNP", "Chromosome", "Position") %in% colnames(observed)))
})

test_that("read_vcf_geno errors when the file is missing or malformed", {
  expect_error(read_vcf_geno("does_not_exist.vcf"), "VCF file not found")

  bad_vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(bad_vcf))
  writeLines(c("##fileformat=VCFv4.2", "1\t100\tsnp1\tA\tG\t.\tPASS\t.\tGT\t0/0"), bad_vcf)

  expect_error(read_vcf_geno(bad_vcf), "#CHROM header line")
})

test_that("read_vcf_geno errors when the VCF is not tab-delimited", {
  space_vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(space_vcf))
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Ind1",
    "1 100 snp1 A G . PASS . GT 0|0"
  ), space_vcf)

  expect_error(read_vcf_geno(space_vcf), "tab-delimited")
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

# Minimal phased VCF fixture, as Beagle would emit it
phased_vcf_fixture <- function() {
  path <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tInd1\tInd2",
    "1\t100\tsnp1\tA\tG\t.\tPASS\t.\tGT\t0|0\t1|1",
    "1\t200\tsnp2\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0",
    "2\t100\tsnp3\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|1"
  ), path)
  path
}

test_that("gt_to_haplotypes splits phased calls into haplotype columns", {
  expect_equal(
    gt_to_haplotypes(c("0|0", "0|1", "1|0", "1|1")),
    matrix(c(0L, 0L, 0L, 1L, 1L, 0L, 1L, 1L), ncol = 2, byrow = TRUE)
  )

  # A valid vector must not warn: a coercion-based implementation would warn here
  expect_silent(gt_to_haplotypes(c("0|0", "1|1")))
})

test_that("gt_to_haplotypes rejects unphased and missing calls", {
  expect_error(gt_to_haplotypes(c("0|1", "0/1")), "phased diploid GT calls")
  expect_error(gt_to_haplotypes(c("0|1", "./.")), "phased diploid GT calls")

  # The phased missing form splits into two alleles, so a count-only check would let it through
  expect_error(gt_to_haplotypes(c("0|1", ".|.")), "phased diploid GT calls")
})

test_that("gt_to_haplotypes errors on calls a 0/1 haplotype cannot express", {
  # Valid VCF at a multi-allelic site, but not a biallelic 0/1 allele presence
  expect_error(gt_to_haplotypes(c("0|1", "0|2")), "Unrecognised GT call")
  expect_error(gt_to_haplotypes(c("0|1", "A|T")), "Unrecognised GT call")
  expect_error(gt_to_haplotypes(c("0|1", NA)), "Unrecognised GT call")
  expect_error(gt_to_haplotypes(c("0|1", "0|1|1")), "Unrecognised GT call")
})

test_that("read_vcf_phased returns the geno_phased layout", {
  vcf_path <- phased_vcf_fixture()
  on.exit(unlink(vcf_path))

  observed <- read_vcf_phased(vcf_path)

  expect_equal(colnames(observed),
               c("SNP", "Chromosome", "Position", "Ind1_1", "Ind1_2", "Ind2_1", "Ind2_2"))
  expect_equal(observed$SNP, c("snp1", "snp2", "snp3"))
  expect_true(all(unlist(observed[, -(1:3)]) %in% c(0L, 1L)))
  expect_equal(observed$Ind1_1, c(0L, 0L, 1L))
  expect_equal(observed$Ind1_2, c(0L, 1L, 1L))
})

test_that("phased haplotypes sum to the dosages read from the same VCF", {
  vcf_path <- phased_vcf_fixture()
  on.exit(unlink(vcf_path))

  phased <- read_vcf_phased(vcf_path)
  dosage <- read_vcf_geno(vcf_path)

  for (s in c("Ind1", "Ind2")) {
    expect_equal(as.numeric(phased[[paste0(s, "_1")]] + phased[[paste0(s, "_2")]]),
                 as.numeric(dosage[[s]]))
  }
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

test_that("gt_to_dosage maps phased and unphased diploid calls, and missing to NA", {
  expect_equal(gt_to_dosage(c("0/0", "0/1", "1/0", "1/1", "./.")), c(0, 1, 1, 2, NA))
  expect_equal(gt_to_dosage(c("0|0", "0|1", "1|0", "1|1", ".|.")), c(0, 1, 1, 2, NA))
  expect_equal(gt_to_dosage(character()), numeric())
})

test_that("gt_to_dosage errors on calls a dosage cannot express", {
  expect_error(gt_to_dosage(c("0/0", "1/2")), "Unrecognised GT call")
  expect_error(gt_to_dosage(c("0/0", "0/0/1/1")), "diploid biallelic")
})
