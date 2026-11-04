# Tests: what each chromosome is doing, beside the bar and as each one finishes -------------
#
# The contract these pin down:
#   - the status lists the chromosomes being worked on, spelt out where there is room and cut
#     back to fit where there is not, and always says how many are done
#   - every stage a chromosome signals is seen, and a line is printed for each one that finishes

# --- the status ----------------------------------------------------------------------------

# One done, two being worked on and one not started, listed in the order given
stages_fixture <- function() c("1" = "done", "7" = "refilling", "12" = "choosing", "3" = NA)

test_that("the status spells out each chromosome being worked on when there is room", {
  expect_equal(chromosome_status(stages_fixture()),
               "Chromosome 7 refilling, Chromosome 12 choosing | 1/4 done")
})

test_that("the status is cut back step by step to fit the room it has", {
  cut_back <- c(
    "Chromosome 7 refilling, Chromosome 12 choosing | 1/4 done",
    "Chr 7 refilling, Chr 12 choosing | 1/4 done",
    "Chr 7 refilling, +1 more | 1/4 done",
    "2 active | 1/4 done",
    "1/4 done"
  )
  # Each exactly fits the width it is given, and a character less gives the next one down
  for (i in seq_along(cut_back)) {
    expect_equal(chromosome_status(stages_fixture(), width = nchar(cut_back[[i]])), cut_back[[i]])
  }
  for (i in seq_len(length(cut_back) - 1L)) {
    expect_equal(chromosome_status(stages_fixture(), width = nchar(cut_back[[i]]) - 1L),
                 cut_back[[i + 1L]])
  }
  # The count of done chromosomes is kept even when it does not fit
  expect_equal(chromosome_status(stages_fixture(), width = 0), "1/4 done")
})

test_that("the status is only the count when nothing is being worked on", {
  expect_equal(chromosome_status(c("1" = NA, "2" = NA)), "0/2 done")
  expect_equal(chromosome_status(c("1" = "done", "2" = "done")), "2/2 done")
})

# --- tracking ------------------------------------------------------------------------------

test_that("every stage is tracked and a line is printed for each finished chromosome", {
  run <- function() {
    report_step("graph_blocking")
    progressr::with_progress({
      p <- progressr::progressor(steps = 2)
      track_chromosome_stages({
        signal_chromosome_stage(1, "building")
        signal_chromosome_stage(2, "building")
        signal_chromosome_stage(1, "done", n_blocks = 5)
        p()
        signal_chromosome_stage(2, "done", n_blocks = 7)
        p()
        "result"
      }, chromosomes = 1:2, describe_done = function(chromosome, d) {
        sprintf("Chromosome %s: %d blocks", chromosome, d$n_blocks)
      }, p = p)
    })
  }
  # Every status rather than the throttled few
  old <- options(HapSelect.progress_interval = 0)
  on.exit(options(old), add = TRUE)
  lines <- NULL
  value <- NULL
  events <- with_machine_logging(function() lines <<- capture_messages(value <<- run()))

  expect_equal(value, "result")
  expect_equal(trimws(lines), c("Chromosome 1: 5 blocks", "Chromosome 2: 7 blocks"))

  statuses <- vapply(Filter(function(e) identical(e$event, "progress") && !is.null(e$message),
                            events), function(e) e$message, character(1))
  # A tick carries the latest status again, so compare each change of status
  changes <- statuses[c(TRUE, statuses[-1] != statuses[-length(statuses)])]
  expect_equal(changes, c(
    "0/2 done",
    "Chromosome 1 building | 0/2 done",
    "Chromosome 1 building, Chromosome 2 building | 0/2 done",
    "Chromosome 2 building | 1/2 done",
    "2/2 done"
  ))

  # A platform gets each finished chromosome's numbers as a metric
  metrics <- Filter(function(e) identical(e$event, "metric"), events)
  expect_equal(vapply(metrics, function(e) e$values$chromosome, character(1)), c("1", "2"))
  expect_equal(vapply(metrics, function(e) e$values$n_blocks, numeric(1)), c(5, 7))
})

test_that("a stage signalled with nothing tracking it is ignored", {
  expect_silent(signal_chromosome_stage(1, "building"))
})

# --- LD blocking ---------------------------------------------------------------------------

test_that("LD blocking prints a line per chromosome and finishes with all of them done", {
  set.seed(3)
  n <- 20L
  g <- matrix(sample(0:2, n * 15L, replace = TRUE), n)
  geno <- data.frame(SNP = paste0("s", seq_len(n)), Chromosome = rep(1:2, each = n / 2),
                     Position = rep(seq_len(n / 2), 2), g, stringsAsFactors = FALSE)
  ld <- pairwise_ld(geno, parallelize = FALSE)

  lines <- NULL
  events <- with_machine_logging(function() {
    lines <<- capture_messages(perform_ld_blocking(ld, geno[, 1:3], threshold = 0.1))
  })

  expect_length(lines, 2)
  expect_match(lines, "^\\S+ Chromosome [12]: \\d+ blocks \\(\\d+ with 2\\+ markers\\)\n$")
  statuses <- Filter(function(e) identical(e$event, "progress") && !is.null(e$message), events)
  expect_equal(statuses[[length(statuses)]]$message, "2/2 done")
  expect_true(any(grepl("Chromosome 1 blocking", vapply(statuses, `[[`, "", "message"))))
})
