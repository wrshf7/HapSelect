# Tests: progress and the future plan across parallel workers ------------------
#
# The contract these pin down:
#   - ticks made on workers reach the main session once each, under the step that made them
#   - a parallel run puts back whatever future plan was in place before it, even when it fails

# Three chromosomes of random dosages, laid out as pairwise_ld() takes them
parallel_ld_fixture <- function(chromosomes = 3L, markers = 15L, individuals = 12L) {
  set.seed(7)
  n <- chromosomes * markers
  g <- matrix(sample(0:2, n * individuals, replace = TRUE), n,
              dimnames = list(NULL, paste0("i", seq_len(individuals))))
  data.frame(SNP = paste0("s", seq_len(n)), Chromosome = rep(seq_len(chromosomes), each = markers),
             Position = rep(seq_len(markers), chromosomes), g, stringsAsFactors = FALSE)
}

# Run body with two workers allowed and a recognisable plan in place, putting both back after.
# The plan is sequential, tweaked only so that it can be told apart from the plain sequential
# plan the package used to reset to
with_caller_plan <- function(body) {
  old_cores <- options(HapSelect.cpu_cores = 2)
  caller_plan <- future::tweak(future::sequential, split = TRUE)
  old_plan <- future::plan(caller_plan)
  on.exit({
    future::plan(old_plan)
    options(old_cores)
  }, add = TRUE)
  force(body)
}

caller_plan_in_place <- function() inherits(future::plan("next"), "tweaked")

test_that("ticks from parallel workers arrive once each, under the step that made them", {
  skip_on_cran()
  geno <- parallel_ld_fixture()

  with_caller_plan({
    old <- options(HapSelect.progress_interval = 0)
    on.exit(options(old), add = TRUE)
    events <- with_machine_logging(function() pairwise_ld(geno, parallelize = TRUE))
  })

  ticks <- Filter(function(e) identical(e$event, "progress"), events)
  expect_gt(length(ticks), 0)
  for (tick in ticks) expect_equal(tick$step, "pairwise_ld")
  # The count is in markers, and never runs past the total or goes backwards
  counted <- Filter(function(e) !is.null(e$current), ticks)
  currents <- vapply(counted, function(e) e$current, numeric(1))
  expect_false(is.unsorted(currents))
  expect_equal(counted[[length(counted)]]$current, nrow(geno))
  expect_equal(counted[[length(counted)]]$total, nrow(geno))
  # The workers' status gives way to what they are doing once they start counting
  expect_equal(ticks[[1]]$message, "Starting 2 workers")
  expect_equal(counted[[length(counted)]]$message, "Calculating LD")
  expect_equal(kinds(events)[length(events)], "step_finished")
})

test_that("a parallel run puts back the plan that was in place before it", {
  skip_on_cran()
  geno <- parallel_ld_fixture(chromosomes = 2L)

  with_caller_plan({
    progressr::with_progress(pairwise_ld(geno, parallelize = TRUE), handlers = progressr::handler_void())
    expect_true(caller_plan_in_place())
  })
})

test_that("a parallel LD blocking run that fails still puts back the plan", {
  skip_on_cran()
  ld <- data.frame(Chrom = c(1, 2), Locus1 = 1, Locus2 = 2, Name1 = c("a", "c"),
                   Name2 = c("b", "d"), LD = 0.9)
  # A map naming none of the LD's markers, so blocking fails on the workers
  map <- data.frame(SNP = "z", Chromosome = 1, Position = 1)

  with_caller_plan({
    expect_error(suppressWarnings(perform_ld_blocking(ld, map, parallel = TRUE)))
    expect_true(caller_plan_in_place())
  })
})

test_that("a progressor sent to workers does not carry HapSelect's imports with it", {
  # future looks up the names inside a progressor from the frame of the HapSelect function
  # sending it to the workers. A name found among HapSelect's imports, such as progressr's own
  # progression(), would send the whole imports environment with every future: tens of
  # megabytes of other packages' functions, which stalled parallel local GEBVs for 20 seconds
  old <- options(progressr.enable = TRUE)
  on.exit(options(old), add = TRUE)
  send_progressor <- function() {
    p <- progressr::progressor(steps = 1)
    # As furrr gathers the globals of the arguments it passes on to each worker
    future::getGlobalsAndPackages(list(p = p), envir = environment())$globals
  }
  environment(send_progressor) <- asNamespace("HapSelect")

  expect_lt(length(serialize(send_progressor(), NULL)), 1e6)
})
