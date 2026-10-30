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

# The counted progress events, in order
counted_progress <- function(events) {
  Filter(function(e) identical(e$event, "progress") && !is.null(e$current), events)
}

test_that("plink_pairwise_ld_geno reports a single plink_ld step, moving through each stage", {
  skip_if_no_plink()
  old <- options(HapSelect.progress_interval = 0)
  on.exit(options(old), add = TRUE)

  events <- with_machine_logging(function() plink_pairwise_ld_geno(plink_progress_genotypes()))
  progress <- counted_progress(events)
  currents <- vapply(progress, function(e) e$current, numeric(1))

  expect_equal(kinds(events)[c(1, length(events))], c("step_started", "step_finished"))
  expect_true(all(steps(events) == "plink_ld"))
  expect_equal(events[[length(events)]]$status, "succeeded")

  # Writing the input, converting it, running PLINK and reading its output, 100 steps each
  expect_equal(progress[[length(progress)]]$total, 400)
  expect_equal(currents[length(currents)], 400)
  expect_false(is.unsorted(currents))
  messages <- unlist(lapply(events, function(e) e$message))
  expect_true(all(c("Writing genotypes for PLINK", "Converting to PLINK format",
                    "Reading LD results") %in% messages))
  expect_true(any(grepl("^Calculating LD [0-9]+%$", messages)))
})


test_that("plink_pairwise_ld moves through running PLINK and reading its output", {
  skip_if_no_plink()

  prefix <- tempfile("plink_progress_")
  HapSelect:::write_plink_ped_map(plink_progress_genotypes(), prefix)
  HapSelect:::run_plink_command(c("--file", prefix, "--make-bed", "--out", prefix,
                                  "--chr-set", "1"))

  events <- with_machine_logging(function() {
    plink_pairwise_ld(prefix, extra_args = c("--chr-set", "1"))
  })
  progress <- counted_progress(events)

  expect_true(all(steps(events) == "plink_ld"))
  expect_equal(events[[length(events)]]$status, "succeeded")
  expect_equal(progress[[length(progress)]]$total, 200)
  expect_equal(progress[[length(progress)]]$current, 200)
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


# Following PLINK's own progress -------------------------------------------------
#
# PLINK redraws a percentage in place for each pass of a command. These feed the parser the
# pieces run_tool_command() hands over for output captured from PLINK 1.9.

# Stages that record where they were moved to, in place of a progressor
recording_stages <- function() {
  seen <- new.env()
  seen$at <- numeric()
  seen$messages <- character()
  list(
    seen = seen,
    within = function(fraction, message = NULL) {
      seen$at <- c(seen$at, fraction)
      seen$messages <- c(seen$messages, message)
    },
    finish = function() NULL
  )
}

# Hand the parser each piece of some output, as run_tool_command() would
feed_plink <- function(output, passes) {
  stages <- recording_stages()
  parse <- HapSelect:::plink_progress_parser(stages, passes = passes)
  pieces <- HapSelect:::complete_pieces(output, "[\r\n\b]+", final = TRUE)$done
  for (piece in pieces[nzchar(pieces)]) parse(piece)
  stages$seen
}

# Captured from --make-bed on a .ped fileset, with most of each pass's percentages left out
plink_make_bed_output <- paste0(
  "Scanning .ped file... 0%\b\b1%\b\b50%\b\b\b100%.ped scan complete (for binary autoconversion).\n",
  "Performing single-pass .bed write (6070 variants, 444 people).\n",
  "0%\b\b50%\b\b\b94%--file: x-temporary.bed + x-temporary.bim + x-temporary.fam written.\n",
  "6070 variants loaded from .bim file.\n",
  "Calculating allele frequencies... 0%\b\b50%\b\b\b done.\n",
  "--make-bed to x.bed + x.bim + x.fam ... 0%\b\b50%\b\b\bdone.\n"
)

test_that("each pass of a PLINK command gets an equal share of its stage", {
  seen <- feed_plink(plink_make_bed_output, passes = 4)

  # A pass starts after a "...", or when the percentage falls back to 0
  expect_equal(seen$at, c(0, 0.01, 0.5, 1, 1, 1.5, 1.94, 2, 2.5, 3, 3.5) / 4)
  # Each pass is named for what it does, rather than in PLINK's words
  expect_equal(unique(sub(" [0-9]+%$", "", seen$messages)),
               c("Reading genotypes", "Converting genotypes", "Calculating allele frequencies",
                 "Writing PLINK files"))
})


test_that("the --r2 pass is followed through its processing and writing redraws", {
  output <- paste0(
    "Calculating allele frequencies... 0%\b\b100%\b\b\b\b done.\n",
    "--r2 to /tmp/x.ld ... 0% [processing]\b\b\b\b\b\b\b\b\b\b\bwriting]   ",
    "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b 37% [processing]\b\b\b\b\b\b\b\b\b\b\bwriting]   ",
    "\b\b\b\b\b\b\b\b\b\b\b\b\b          \b\b\b\b\b\b\b\b\b\b\b\bdone.\n"
  )
  seen <- feed_plink(output, passes = 2)

  expect_equal(seen$at, c(0, 1, 1, 1.37) / 2)
  expect_equal(seen$messages[3:4], c("Calculating LD 0%", "Calculating LD 37%"))
})


test_that("a pass PLINK names in a way it is not known by still reads as PLINK running", {
  seen <- feed_plink("Pruning variants... 0%\b\b40%\n", passes = 1)
  expect_equal(seen$messages, c("Running PLINK 0%", "Running PLINK 40%"))
})


test_that("output that is not a percentage moves nothing", {
  seen <- feed_plink("PLINK v1.9.0-rc2 64-bit\nLogging to x.log.\nError: no such file\n", passes = 2)
  expect_length(seen$at, 0)
})


test_that("PLINK stages only move forward, and never past their share", {
  ticks <- numeric()
  messages <- character()
  p <- function(amount = 1, message = NULL) {
    ticks <<- c(ticks, amount)
    messages <<- c(messages, message)
  }
  stages <- HapSelect:::plink_stages(p)

  stages$within(0.5, "half")
  stages$within(0.25, "back")  # a step back only updates the message
  stages$within(3)             # more passes than expected stops at the stage's end
  stages$finish()
  stages$within(0.1)
  stages$finish()

  expect_equal(ticks, c(50, 0, 50, 10, 90))
  expect_equal(sum(ticks), 2 * HapSelect:::PLINK_STAGE_STEPS)
  expect_equal(messages, c("half", "back"))
})


test_that("the .ped file is written a chunk of individuals at a time, reporting as it goes", {
  geno <- plink_progress_genotypes()
  wide <- cbind(geno[, 1:3], do.call(cbind, rep(list(geno[, -(1:3)]), 24)))
  names(wide)[-(1:3)] <- paste0("ind", seq_len(ncol(wide) - 3))
  prefix <- tempfile("ped_chunks_")
  on.exit(unlink(paste0(prefix, c(".ped", ".map"))), add = TRUE)

  fractions <- numeric()
  HapSelect:::write_plink_ped_map(wide, prefix, progress = function(f) fractions <<- c(fractions, f))

  ped <- utils::read.table(paste0(prefix, ".ped"), colClasses = "character")
  expect_equal(nrow(ped), 120)
  expect_equal(ped$V2, paste0("I", 1:120))
  # Individual 6 is ind1's dosages again: 0, 0, 2, 0
  expect_equal(unlist(ped[6, 7:14], use.names = FALSE), c("A", "A", "A", "A", "G", "G", "A", "A"))
  expect_gt(length(fractions), 1)
  expect_false(is.unsorted(fractions))
  expect_equal(fractions[length(fractions)], 1)
})
