test_that("format_plink_ld converts PLINK output to FastStack LD format", {
  tmp_dir <- tempfile("plink_ld_test_")
  dir.create(tmp_dir)
  prefix <- file.path(tmp_dir, "toy")
  bim_path <- paste0(prefix, ".bim")
  ld_path <- file.path(tmp_dir, "toy.ld")

  writeLines(
    c(
      "1 snp1 0 100 A G",
      "1 snp2 0 200 A G",
      "1 snp3 0 300 A G",
      "2 snp4 0 100 A G",
      "2 snp5 0 200 A G"
    ),
    bim_path
  )

  writeLines(
    c(
      "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2",
      "1 100 snp1 1 200 snp2 0.81",
      "1 100 snp1 1 300 snp3 0.25",
      "2 100 snp4 2 200 snp5 0.64"
    ),
    ld_path
  )

  bim <- FastStack:::read_plink_bim(bim_path)
  observed <- format_plink_ld(ld_path, bim_path)

  expected <- data.frame(
    Chrom = c(1, 1, 2),
    Locus1 = c(1, 1, 1),
    Locus2 = c(2, 3, 2),
    Name1 = c("snp1", "snp1", "snp4"),
    Name2 = c("snp2", "snp3", "snp5"),
    LD = c(0.81, 0.25, 0.64),
    stringsAsFactors = FALSE
  )

  expect_equal(bim$Locus, c(1, 2, 3, 1, 2))
  expect_equal(observed, expected)
})


test_that("plink_pairwise_ld errors when the PLINK fileset is incomplete", {
  expect_error(
    plink_pairwise_ld("missing_prefix"),
    "Missing required PLINK input files"
  )
})
