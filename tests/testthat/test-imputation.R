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


# Progress reporting ------------------------------------------------------------
#
# Beagle announces each window it works on, with the chromosome and position range, but never
# how many windows there will be. Progress is measured in markers instead: once the next window
# is announced, every marker up to the end of the previous one is done.

# A progressor stand-in that records each amount it is ticked by
recording_progressor <- function() {
  amounts <- integer()
  p <- function(amount = 1, ...) amounts <<- c(amounts, as.integer(amount))
  list(p = p, amounts = function() amounts)
}

beagle_test_markers <- function() {
  data.frame(Chromosome = c(rep(1, 5), rep(2, 3)),
             Position = c(100, 200, 300, 400, 500, 100, 200, 300))
}

test_that("Beagle progress completes a window once the next one is announced", {
  rec <- recording_progressor()
  parse <- beagle_progress_parser(beagle_test_markers(), rec$p)

  parse("Window 1 [1:100-300]")
  expect_equal(rec$amounts(), integer())

  # Window 1 ends at the third marker, so the second window's header completes three
  parse("Burnin  iteration 1:           0 seconds")
  parse("Window 2 [1:250-500]")
  expect_equal(rec$amounts(), 3L)

  parse("Window 3 [2:100-300]")
  expect_equal(rec$amounts(), c(3L, 2L))

  # The closing statistics complete the last window, and everything with it
  parse("Cumulative Statistics:")
  expect_equal(sum(rec$amounts()), 8L)
})


test_that("Beagle progress never steps back over an overlapping window", {
  rec <- recording_progressor()
  parse <- beagle_progress_parser(beagle_test_markers(), rec$p)

  parse("Window 1 [1:100-400]")
  parse("Window 2 [1:350-380]")
  # Window 2 ends inside window 1, so finishing it moves nothing
  parse("Window 3 [1:380-500]")
  expect_equal(rec$amounts(), 4L)
})


test_that("Beagle progress reads a window end lying between markers", {
  rec <- recording_progressor()
  parse <- beagle_progress_parser(beagle_test_markers(), rec$p)

  # 250 falls between the second and third markers, so two are done
  parse("Window 1 [1:100-250]")
  parse("Window 2 [1:240-500]")
  expect_equal(rec$amounts(), 2L)
})


test_that("Beagle progress ignores other output and unknown chromosomes", {
  rec <- recording_progressor()
  parse <- beagle_progress_parser(beagle_test_markers(), rec$p)

  parse("beagle.27Feb25.75f.jar (version 5.5)")
  parse("Study     markers:                  251")
  parse("Window 1 [9:100-300]")
  parse("Window 2 [1:100-300]")
  expect_equal(rec$amounts(), integer())
})


test_that("Beagle progress parses a chromosome name containing a colon", {
  rec <- recording_progressor()
  markers <- data.frame(Chromosome = "HLA:1", Position = c(10, 20, 30))
  parse <- beagle_progress_parser(markers, rec$p)

  parse("Window 1 [HLA:1:10-20]")
  parse("Window 2 [HLA:1:15-30]")
  expect_equal(rec$amounts(), 2L)
})


test_that("the Beagle wrappers report a rejected input as their own failed step", {
  events <- with_machine_logging(function() {
    expect_error(beagle_impute_geno("not_a_data_frame"), "geno must be a data frame")
  })
  expect_equal(steps(events), c("beagle_impute", "beagle_impute"))
  expect_equal(events[[2]]$status, "failed")

  events <- with_machine_logging(function() {
    expect_error(beagle_phase_geno("not_a_data_frame"), "geno must be a data frame")
  })
  expect_equal(steps(events), c("beagle_phase", "beagle_phase"))
  expect_equal(events[[2]]$status, "failed")

  events <- with_machine_logging(function() {
    expect_error(beagle_impute("does_not_exist.vcf", tempfile()), "Beagle input VCF not found")
  })
  expect_equal(steps(events), c("beagle_impute", "beagle_impute"))
  expect_equal(events[[2]]$status, "failed")
})


# Real Beagle runs ---------------------------------------------------------------

skip_if_no_beagle <- function() {
  skip_if(
    inherits(tryCatch(HapSelect:::find_beagle(), error = function(e) e), "error"),
    "Beagle not found"
  )
}

# Two chromosomes of 40 markers 0.5 Mb apart, which Beagle's default of 1 cM = 1 Mb and a
# 4 cM window splits into several windows each. The positions are round, which write_vcf_geno()
# once wrote as 5e+05 and Beagle could not read.
beagle_run_fixture <- function() {
  set.seed(42)
  n <- 40
  geno <- data.frame(SNP = paste0("s", seq_len(2 * n)), Chromosome = rep(1:2, each = n),
                     Position = rep(seq(500000, by = 500000, length.out = n), 2))
  for (i in 1:30) {
    geno[[paste0("ind", i)]] <- sample(c(0, 1, 2, NA), 2 * n, replace = TRUE,
                                       prob = c(0.3, 0.4, 0.25, 0.05))
  }
  geno
}

test_that("beagle_phase_geno reports its progress as Beagle runs", {
  skip_if_no_beagle()
  geno <- beagle_run_fixture()

  # Every tick reported rather than throttled, so the intermediate ones can be seen
  old <- options(HapSelect.progress_interval = 0)
  on.exit(options(old), add = TRUE)

  events <- with_machine_logging(function() {
    beagle_phase_geno(geno, extra_args = c("window=4", "overlap=1"))
  })
  ticks <- Filter(function(e) identical(e$event, "progress"), events)
  current <- vapply(ticks, function(e) e$current, numeric(1))

  expect_equal(kinds(events)[c(1, length(events))], c("step_started", "step_finished"))
  expect_true(all(steps(events) == "beagle_phase"))
  expect_equal(events[[length(events)]]$status, "succeeded")

  expect_equal(ticks[[length(ticks)]]$total, nrow(geno))
  # Progress arrived while Beagle ran, rather than all at once when it finished
  expect_true(any(current > 0 & current < nrow(geno)))
  expect_false(is.unsorted(current))
})


test_that("a failed Beagle run stops with Beagle's output and a failed step", {
  skip_if_no_beagle()
  geno <- beagle_run_fixture()
  # Beagle refuses a marker repeated at the same position
  geno <- rbind(geno[1, ], geno)

  events <- with_machine_logging(function() {
    expect_error(beagle_phase_geno(geno), "Beagle imputation/phasing failed[\\s\\S]*Duplicate marker",
                 perl = TRUE)
  })

  expect_equal(events[[length(events)]]$step, "beagle_phase")
  expect_equal(events[[length(events)]]$status, "failed")
})
