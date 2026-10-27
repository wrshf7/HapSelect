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


test_that("format_plink_ld keeps marker names that look like numbers exactly", {
  # Read as numbers, 007 would become 7 and 1e5 would become 1e+05, and no longer match
  tmp_dir <- tempfile("plink_ld_names_")
  dir.create(tmp_dir)
  bim_path <- file.path(tmp_dir, "toy.bim")
  ld_path <- file.path(tmp_dir, "toy.ld")

  writeLines(c("1 007 0 100 A G", "1 1e5 0 200 A G", "1 snp3 0 300 A G"), bim_path)
  writeLines(
    c(
      " CHR_A         BP_A SNP_A  CHR_B         BP_B SNP_B           R2 ",
      "     1          100   007      1          200   1e5         0.81 ",
      "     1          200   1e5      1          300  snp3         0.25 "
    ),
    ld_path
  )

  observed <- format_plink_ld(ld_path, bim_path)

  expect_identical(observed$Name1, c("007", "1e5"))
  expect_identical(observed$Name2, c("1e5", "snp3"))
  expect_equal(observed$Locus1, c(1, 2))
  expect_equal(observed$LD, c(0.81, 0.25))
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

skip_if_no_plink <- function() {
  skip_if(
    inherits(tryCatch(HapSelect:::find_plink(), error = function(e) e), "error"),
    "PLINK executable not found"
  )
}

plink_progress_genotypes <- function() {
  data.frame(
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
}

test_that("plink_pairwise_ld_geno reports a single plink_ld step, ticking at each stage", {
  skip_if_no_plink()

  n_ticks <- NULL
  events <- with_machine_logging(function() {
    n_ticks <<- count_ticks(function() plink_pairwise_ld_geno(plink_progress_genotypes()))
  })
  progress <- Filter(function(e) identical(e$event, "progress"), events)

  expect_equal(kinds(events)[c(1, length(events))], c("step_started", "step_finished"))
  expect_true(all(steps(events) == "plink_ld"))
  expect_equal(events[[length(events)]]$status, "succeeded")

  # Writing the input, converting it, running PLINK and reading its output
  expect_equal(progress[[length(progress)]]$total, 4)
  expect_equal(n_ticks, 4)
})


test_that("plink_pairwise_ld ticks once PLINK has run and once its output is read", {
  skip_if_no_plink()

  prefix <- tempfile("plink_progress_")
  HapSelect:::write_plink_ped_map(plink_progress_genotypes(), prefix)
  HapSelect:::run_plink_command(c("--file", prefix, "--make-bed", "--out", prefix,
                                  "--chr-set", "1"))

  n_ticks <- NULL
  events <- with_machine_logging(function() {
    n_ticks <<- count_ticks(function() {
      plink_pairwise_ld(prefix, extra_args = c("--chr-set", "1"))
    })
  })

  expect_true(all(steps(events) == "plink_ld"))
  expect_equal(events[[length(events)]]$status, "succeeded")
  expect_equal(n_ticks, 2)
})


test_that("a failed PLINK run stops with PLINK's output and a failed step", {
  skip_if_no_plink()

  events <- with_machine_logging(function() {
    expect_error(
      plink_pairwise_ld_geno(plink_progress_genotypes(), extra_args = "--not-a-real-flag"),
      "PLINK LD calculation failed[\\s\\S]*not-a-real-flag", perl = TRUE
    )
  })

  expect_equal(events[[length(events)]]$step, "plink_ld")
  expect_equal(events[[length(events)]]$status, "failed")
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
