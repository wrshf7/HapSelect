# Tests: VCF reading and writing ----------------------------------------------
#
# Covers R/vcf.R: the readers and writer themselves, and the GT call converters
# they are built on. The Beagle wrappers that call them are covered separately
# in test-imputation.R.

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

# Minimal phased VCF fixture, as phasing software would emit it
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
  expect_error(gt_to_haplotypes(c("0|1", "0/1")), "Expected phased calls of 2 alleles")
  expect_error(gt_to_haplotypes(c("0|1", "./.")), "Expected phased calls of 2 alleles")

  # The phased missing form splits into two alleles, so a count-only check would let it through
  expect_error(gt_to_haplotypes(c("0|1", ".|.")), "Expected phased calls of 2 alleles")
})

test_that("gt_to_haplotypes errors on calls a 0/1 haplotype cannot express", {
  # Valid VCF at a multi-allelic site, but not a biallelic 0/1 allele presence
  expect_error(gt_to_haplotypes(c("0|1", "0|2")), "Unrecognised GT call")
  expect_error(gt_to_haplotypes(c("0|1", "A|T")), "Unrecognised GT call")
  expect_error(gt_to_haplotypes(c("0|1", NA)), "Unrecognised GT call")
  expect_error(gt_to_haplotypes(c("0|1", "0|1|1")), "Unrecognised GT call")
})

# A VCF whose GT calls are whatever is passed in, for the reader error tests
gt_vcf_fixture <- function(calls) {
  path <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tInd1",
    paste0("1\t100\tsnp1\tA\tG\t.\tPASS\t.\tGT\t", calls)
  ), path)
  path
}

test_that("read_vcf_phased explains what to do when the VCF is not phased", {
  # write_vcf_geno() emits unphased calls, so a round-tripped geno file is exactly this mistake
  geno <- data.frame(
    SNP = "snp1", Chromosome = 1, Position = 100, Ind1 = 1, Ind2 = 0,
    stringsAsFactors = FALSE
  )
  vcf_path <- tempfile(fileext = ".vcf")
  on.exit(unlink(vcf_path))
  HapSelect:::write_vcf_geno(geno, vcf_path)

  # fixed = TRUE throughout: "0/1" as a regex is an alternation that would match anything
  expect_error(read_vcf_phased(vcf_path), "0/1", fixed = TRUE)
  expect_error(read_vcf_phased(vcf_path), vcf_path, fixed = TRUE)
  expect_error(read_vcf_phased(vcf_path), "Ind1", fixed = TRUE)
  expect_error(read_vcf_phased(vcf_path), "Phase the genotypes first", fixed = TRUE)
})

test_that("read_vcf_phased surfaces the offending call for missing and multi-allelic VCFs", {
  missing_vcf <- gt_vcf_fixture(".|.")
  multi_vcf <- gt_vcf_fixture("0|2")
  on.exit(unlink(c(missing_vcf, multi_vcf)))

  expect_error(read_vcf_phased(missing_vcf), ".|.", fixed = TRUE)
  expect_error(read_vcf_phased(multi_vcf), "0|2", fixed = TRUE)
  expect_error(read_vcf_phased(multi_vcf), multi_vcf, fixed = TRUE)
})

test_that("read_vcf_geno names the file when a call cannot be a dosage", {
  multi_vcf <- gt_vcf_fixture("1/2")
  on.exit(unlink(multi_vcf))

  expect_error(read_vcf_geno(multi_vcf), "1/2", fixed = TRUE)
  expect_error(read_vcf_geno(multi_vcf), multi_vcf, fixed = TRUE)
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

test_that("gt_to_dosage maps phased and unphased diploid calls, and missing to NA", {
  expect_equal(gt_to_dosage(c("0/0", "0/1", "1/0", "1/1", "./.")), c(0, 1, 1, 2, NA))
  expect_equal(gt_to_dosage(c("0|0", "0|1", "1|0", "1|1", ".|.")), c(0, 1, 1, 2, NA))
  expect_equal(gt_to_dosage(character()), numeric())
})

test_that("gt_to_dosage errors on calls a dosage cannot express", {
  expect_error(gt_to_dosage(c("0/0", "1/2")), "Unrecognised GT call")
  expect_error(gt_to_dosage(c("0/0", "0/0/1/1")), "wrong number of alleles for ploidy 2")
})

# A VCF carrying DS alongside GT, as imputation software emits. ds_calls holds one
# "GT:DS" string per record for the single sample, or a bare GT to leave DS off that record.
ds_vcf_fixture <- function(ds_calls, format = "GT:DS") {
  path <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated ALT dose\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tInd1",
    paste0("1\t", seq_along(ds_calls) * 100, "\tsnp", seq_along(ds_calls),
           "\tA\tG\t.\tPASS\t.\t", format, "\t", ds_calls)
  ), path)
  path
}

test_that("read_vcf_geno reads fractional dosages from DS when prefer_ds is TRUE", {
  vcf_path <- ds_vcf_fixture(c("0|0:0.02", "0|1:0.83", "1|1:1.96"))
  on.exit(unlink(vcf_path))

  observed <- read_vcf_geno(vcf_path, prefer_ds = TRUE)

  # The fractions survive: rounding to the GT hard calls would give 0, 1, 2
  expect_equal(observed$Ind1, c(0.02, 0.83, 1.96))
})

test_that("read_vcf_geno ignores DS unless asked, so existing callers are unaffected", {
  vcf_path <- ds_vcf_fixture(c("0|0:0.02", "0|1:0.83", "1|1:1.96"))
  on.exit(unlink(vcf_path))

  expect_equal(read_vcf_geno(vcf_path)$Ind1, c(0, 1, 2))
  expect_equal(read_vcf_geno(vcf_path, prefer_ds = FALSE)$Ind1, c(0, 1, 2))
})

test_that("read_vcf_geno falls back to GT when any record is missing DS", {
  # Middle record carries GT only, so the file is read as GT throughout rather than half and half
  vcf_path <- ds_vcf_fixture(c("0|0:0.02", "0|1", "1|1:1.96"), format = c("GT:DS", "GT", "GT:DS"))
  on.exit(unlink(vcf_path))

  expect_equal(read_vcf_geno(vcf_path, prefer_ds = TRUE)$Ind1, c(0, 1, 2))
})

test_that("read_vcf_geno reads a missing DS as NA", {
  vcf_path <- ds_vcf_fixture(c("0|0:0.02", ".|.:."))
  on.exit(unlink(vcf_path))

  expect_equal(read_vcf_geno(vcf_path, prefer_ds = TRUE)$Ind1, c(0.02, NA))
})

test_that("ds_to_dosage rejects values that are not ALT dosages at the default ploidy", {
  expect_equal(ds_to_dosage(c("0.00", "0.83", ".")), c(0, 0.83, NA))

  expect_error(ds_to_dosage(c("0.5", "0/1")), "Unrecognised DS value")
  expect_error(ds_to_dosage(c("0.5", "2.4")), "0 to 2 dosage range")
  expect_error(ds_to_dosage(c("0.5", "-1")), "0 to 2 dosage range")
})


test_that("a bare . is read as a missing genotype", {
  # "./." is the diploid spelling, "." the VCF spelling for a genotype that was
  # not called at all. Both mean the same thing to a dosage.
  expect_equal(gt_to_dosage(c("0/0", ".", "1/1", "./.")), c(0, NA, 2, NA))

  vcf_path <- gt_vcf_fixture(".")
  on.exit(unlink(vcf_path))
  expect_equal(read_vcf_geno(vcf_path)$Ind1, NA_real_)
})


test_that("a VCF with Windows line endings reads the same as one without", {
  # readLines() accepts LF, CRLF or CR as the terminator whatever mode the
  # connection is opened in, so the carriage returns never reach the parser. This
  # pins that, because the last field of every line is what would break first.
  body <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tInd1\tInd2",
    "1\t100\tsnp1\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1",
    "1\t200\tsnp2\tA\tG\t.\tPASS\t.\tGT\t0/1\t./."
  )

  unix_path <- tempfile(fileext = ".vcf")
  crlf_path <- tempfile(fileext = ".vcf")
  on.exit(unlink(c(unix_path, crlf_path)))

  writeLines(body, unix_path)
  # writeBin rather than writeLines, which would use the platform's own ending
  con <- file(crlf_path, "wb")
  writeBin(charToRaw(paste0(paste(body, collapse = "\r\n"), "\r\n")), con)
  close(con)

  expect_silent(crlf <- read_vcf_geno(crlf_path))
  expect_equal(crlf, read_vcf_geno(unix_path))

  # the final sample column is the one a stray carriage return would land on
  expect_equal(crlf$Ind2, c(2, NA))
})


# Tests: ploidy ----------------------------------------------------------------
#
# Every reader and the writer take ploidy as declared. Dosages run from 0 to ploidy,
# GT calls must carry exactly ploidy alleles, and phased tables get ploidy haplotype
# columns per sample. The default of 2L keeps the diploid behaviour above unchanged.

# A single-marker geno data frame holding each dosage 0 to ploidy, then a missing one
ploidy_geno_fixture <- function(ploidy) {
  dosages <- c(0:ploidy, NA)
  geno <- data.frame(SNP = "snp1", Chromosome = 1, Position = 100, stringsAsFactors = FALSE)
  for (i in seq_along(dosages)) geno[[paste0("Ind", i)]] <- dosages[i]
  geno
}

# The GT calls a written VCF holds for its first record, in sample order
written_calls <- function(path) {
  records <- HapSelect:::read_vcf_records(path)
  unname(records$tab[1, records$sample_cols])
}

test_that("check_ploidy accepts whole numbers of at least 1 and rejects anything else", {
  expect_identical(HapSelect:::check_ploidy(1L), 1L)
  expect_identical(HapSelect:::check_ploidy(2), 2L)
  expect_identical(HapSelect:::check_ploidy(4L), 4L)
  expect_identical(HapSelect:::check_ploidy(6L), 6L)

  for (bad in list(0L, -1L, 2.5, NA, NA_integer_, Inf, c(2L, 4L), "2", integer())) {
    expect_error(HapSelect:::check_ploidy(bad), "ploidy must be a single whole number")
  }
})

test_that("write_vcf_geno writes REF alleles then ALT alleles at every ploidy", {
  expected <- list(
    "1" = c("0", "1", "."),
    "2" = c("0/0", "0/1", "1/1", "./."),
    "3" = c("0/0/0", "0/0/1", "0/1/1", "1/1/1", "././."),
    "4" = c("0/0/0/0", "0/0/0/1", "0/0/1/1", "0/1/1/1", "1/1/1/1", "./././.")
  )

  for (ploidy in 1:4) {
    path <- tempfile(fileext = ".vcf")
    HapSelect:::write_vcf_geno(ploidy_geno_fixture(ploidy), path, ploidy = ploidy)
    expect_equal(written_calls(path), expected[[as.character(ploidy)]], info = paste("ploidy", ploidy))
    unlink(path)
  }
})

test_that("write_vcf_geno rejects dosages the declared ploidy cannot hold", {
  geno <- data.frame(SNP = "snp1", Chromosome = 1, Position = 100, Ind1 = 0, stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path))

  geno$Ind1 <- 3
  expect_error(HapSelect:::write_vcf_geno(geno, path), "0 to 2, or NA, for ploidy 2")
  geno$Ind1 <- 5
  expect_error(HapSelect:::write_vcf_geno(geno, path, ploidy = 4L), "0 to 4, or NA, for ploidy 4")
  geno$Ind1 <- 1.5
  expect_error(HapSelect:::write_vcf_geno(geno, path, ploidy = 4L), "whole numbers")
  geno$Ind1 <- -1
  expect_error(HapSelect:::write_vcf_geno(geno, path, ploidy = 4L), "whole numbers")
})

test_that("write_vcf_geno and read_vcf_geno round-trip dosages at every ploidy", {
  for (ploidy in c(1L, 2L, 3L, 4L, 6L)) {
    geno <- ploidy_geno_fixture(ploidy)
    path <- tempfile(fileext = ".vcf")
    HapSelect:::write_vcf_geno(geno, path, ploidy = ploidy)

    observed <- read_vcf_geno(path, ploidy = ploidy)
    expect_equal(unlist(observed[, -(1:3)]), unlist(geno[, -(1:3)]), info = paste("ploidy", ploidy))
    unlink(path)
  }
})

test_that("gt_to_dosage counts ALT alleles in triploid calls", {
  expect_equal(gt_to_dosage(c("0/0/0", "0/0/1", "1|1|0", "1/1/1"), ploidy = 3L), c(0, 1, 2, 3))
})

test_that("gt_to_dosage counts ALT alleles in tetraploid calls, whatever the separators", {
  expect_equal(
    gt_to_dosage(c("0/0/0/0", "0/0/0/1", "0/1/0/1", "1|1|1|0", "1/1/1/1"), ploidy = 4L),
    c(0, 1, 2, 3, 4)
  )
  # A partly phased call is still a count of ALT alleles
  expect_equal(gt_to_dosage("0/1|1/0", ploidy = 4L), 2)
})

test_that("gt_to_dosage reads a call with any unknown allele as NA", {
  expect_equal(gt_to_dosage(c("0/.", "./1", "1|."), ploidy = 2L), c(NA_real_, NA_real_, NA_real_))
  expect_equal(gt_to_dosage(c("0/0/./1", "./././.", ".", "0/0/1/1"), ploidy = 4L),
               c(NA, NA, NA, 2))
  expect_equal(gt_to_dosage(c(".", "0/0/1"), ploidy = 3L), c(NA, 1))
})

test_that("gt_to_dosage errors when a call does not carry ploidy alleles", {
  expect_error(gt_to_dosage("0/0/1/1", ploidy = 2L), "wrong number of alleles for ploidy 2")
  expect_error(gt_to_dosage("0/1", ploidy = 4L), "wrong number of alleles for ploidy 4")
  expect_error(gt_to_dosage("0/0/1", ploidy = 4L), "0/0/1", fixed = TRUE)
  # The allele count applies to a partly missing call too
  expect_error(gt_to_dosage("./.", ploidy = 4L), "wrong number of alleles for ploidy 4")
})

test_that("gt_to_dosage errors on multi-allelic and malformed calls at any ploidy", {
  expect_error(gt_to_dosage("0/0/1/2", ploidy = 4L), "multi-allelic calls are not supported")
  expect_error(gt_to_dosage("0/2/1", ploidy = 3L), "Unrecognised GT call")
  expect_error(gt_to_dosage(c("0/0", NA)), "Unrecognised GT call")
  expect_error(gt_to_dosage(c("0/0", "")), "Unrecognised GT call")
})

test_that("read_vcf_geno names the file when its calls do not match the declared ploidy", {
  vcf_path <- gt_vcf_fixture("0/0/1/1")
  on.exit(unlink(vcf_path))

  expect_error(read_vcf_geno(vcf_path), "wrong number of alleles for ploidy 2")
  expect_error(read_vcf_geno(vcf_path), vcf_path, fixed = TRUE)
  expect_equal(read_vcf_geno(vcf_path, ploidy = 4L)$Ind1, 2)
})

test_that("ds_to_dosage bounds DS by the declared ploidy", {
  expect_equal(ds_to_dosage(c("0.2", "3.6", "4"), ploidy = 4L), c(0.2, 3.6, 4))
  expect_error(ds_to_dosage("4.2", ploidy = 4L), "0 to 4 dosage range for ploidy 4")
  expect_error(ds_to_dosage("2.4"), "0 to 2 dosage range for ploidy 2")
  expect_error(ds_to_dosage("-0.1", ploidy = 4L), "0 to 4 dosage range")
})

test_that("read_vcf_geno passes ploidy to the DS range check", {
  vcf_path <- ds_vcf_fixture(c("0/0/0/1:0.9", "0/1/1/1:3.4"))
  on.exit(unlink(vcf_path))

  expect_equal(read_vcf_geno(vcf_path, prefer_ds = TRUE, ploidy = 4L)$Ind1, c(0.9, 3.4))
  expect_error(read_vcf_geno(vcf_path, prefer_ds = TRUE), "0 to 2 dosage range")
})

test_that("gt_to_haplotypes splits triploid and tetraploid phased calls", {
  expect_equal(
    gt_to_haplotypes(c("0|1|1", "1|0|0"), ploidy = 3L),
    matrix(c(0L, 1L, 1L, 1L, 0L, 0L), ncol = 3, byrow = TRUE)
  )
  expect_equal(
    gt_to_haplotypes(c("1|0|0|1", "0|0|0|0", "1|0|0|1"), ploidy = 4L),
    matrix(c(1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L), ncol = 4, byrow = TRUE)
  )
})

test_that("gt_to_haplotypes reads haploid calls, which have nothing to phase", {
  expect_equal(gt_to_haplotypes(c("0", "1", "1"), ploidy = 1L), matrix(c(0L, 1L, 1L), ncol = 1))
})

test_that("gt_to_haplotypes rejects polyploid calls that are not fully phased and complete", {
  # Unphased, partly phased, and unphased homozygous calls are all rejected
  expect_error(gt_to_haplotypes("0/1/1", ploidy = 3L), "Expected phased calls of 3 alleles")
  expect_error(gt_to_haplotypes("0|1/1|0", ploidy = 4L), "Expected phased calls of 4 alleles")
  expect_error(gt_to_haplotypes("0/0/0/0", ploidy = 4L), "Expected phased calls of 4 alleles")

  expect_error(gt_to_haplotypes("0|.|1|1", ploidy = 4L), "0|.|1|1", fixed = TRUE)
  expect_error(gt_to_haplotypes("0|1|1", ploidy = 4L), "Expected phased calls of 4 alleles")
  expect_error(gt_to_haplotypes("0|1|2|1", ploidy = 4L), "Unrecognised GT call")
})

# A phased tetraploid VCF with two samples and three records
tetraploid_phased_fixture <- function() {
  path <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tInd1\tInd2",
    "1\t100\tsnp1\tA\tG\t.\tPASS\t.\tGT\t0|0|0|1\t1|1|1|1",
    "1\t200\tsnp2\tA\tG\t.\tPASS\t.\tGT\t0|1|1|0\t0|0|0|0",
    "2\t100\tsnp3\tA\tG\t.\tPASS\t.\tGT\t1|1|1|0\t1|0|1|0"
  ), path)
  path
}

test_that("read_vcf_phased gives ploidy haplotype columns per sample", {
  vcf_path <- tetraploid_phased_fixture()
  on.exit(unlink(vcf_path))

  observed <- read_vcf_phased(vcf_path, ploidy = 4L)

  expect_equal(colnames(observed),
               c("SNP", "Chromosome", "Position", paste0("Ind1_", 1:4), paste0("Ind2_", 1:4)))
  expect_equal(observed$Ind1_4, c(1L, 0L, 0L))
  expect_equal(observed$Ind2_3, c(1L, 0L, 1L))

  # Read at the default ploidy, the same file names the mismatch and the file
  expect_error(read_vcf_phased(vcf_path), "Expected phased calls of 2 alleles")
  expect_error(read_vcf_phased(vcf_path), vcf_path, fixed = TRUE)
})

test_that("polyploid haplotypes sum to the dosages read from the same VCF", {
  vcf_path <- tetraploid_phased_fixture()
  on.exit(unlink(vcf_path))

  phased <- read_vcf_phased(vcf_path, ploidy = 4L)
  dosage <- read_vcf_geno(vcf_path, ploidy = 4L)

  for (s in c("Ind1", "Ind2")) {
    hap_sum <- Reduce(`+`, phased[paste0(s, "_", 1:4)])
    expect_equal(as.numeric(hap_sum), as.numeric(dosage[[s]]), info = s)
  }
})

test_that("a triploid VCF written as dosages reads as dosages but not as haplotypes", {
  geno <- ploidy_geno_fixture(3L)
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path))
  HapSelect:::write_vcf_geno(geno, path, ploidy = 3L)

  # Written unphased, so the dosage reader accepts it and the phased reader refuses it
  expect_equal(unlist(read_vcf_geno(path, ploidy = 3L)[, -(1:3)]), unlist(geno[, -(1:3)]))
  expect_error(read_vcf_phased(path, ploidy = 3L), "Phase the genotypes first")
})
