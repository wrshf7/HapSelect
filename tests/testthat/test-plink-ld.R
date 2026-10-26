test_that("format_plink_ld converts PLINK output to HapSelect LD format", {
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

  bim <- HapSelect:::read_plink_bim(bim_path)
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

test_that("plink_pairwise_ld_geno errors when geno is not a valid data frame", {
  expect_error(
    plink_pairwise_ld_geno("not_a_data_frame"),
    "geno must be a data frame"
  )
  expect_error(
    plink_pairwise_ld_geno(data.frame(a = 1, b = 2, c = 3)),
    "geno must be a data frame"
  )
})


# Progress reporting ------------------------------------------------------------
#
# Both entry points report the one plink_ld step, so a call shows one step whichever
# input it starts from, and the genotype path's own preparation falls inside it.

test_that("plink_pairwise_ld_geno reports a single plink_ld step", {
  skip_if(
    inherits(tryCatch(HapSelect:::find_plink(), error = function(e) e), "error"),
    "PLINK executable not found"
  )

  genotypes <- data.frame(
    marker = paste0("m", 1:4),
    chrom = rep(1, 4),
    pos = seq(100, 400, by = 100),
    ind1 = c(0, 0, 2, 0),
    ind2 = c(0, 1, 2, 0),
    ind3 = c(1, 1, 1, 0),
    ind4 = c(1, 2, 1, 1),
    ind5 = c(2, 2, 0, 1),
    stringsAsFactors = FALSE
  )

  events <- with_machine_logging(function() plink_pairwise_ld_geno(genotypes))

  expect_equal(kinds(events), c("step_started", "step_finished"))
  expect_equal(steps(events), c("plink_ld", "plink_ld"))
  expect_equal(events[[2]]$status, "succeeded")
})


test_that("plink_pairwise_ld_geno reports a rejected input as a failed plink_ld step", {
  events <- with_machine_logging(function() {
    expect_error(plink_pairwise_ld_geno("not_a_data_frame"), "geno must be a data frame")
  })

  expect_equal(kinds(events), c("step_started", "step_finished"))
  expect_equal(steps(events), c("plink_ld", "plink_ld"))
  expect_equal(events[[2]]$status, "failed")
  expect_match(events[[2]]$message, "geno must be a data frame")
})
