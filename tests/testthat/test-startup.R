# Tests: package startup -------------------------------------------------------
#
# The startup banner belongs to the session that attached the package. A parallel worker
# attaches it too, and future relays the worker's output back, so the banner has to know
# when it is running in one.

test_that("in_parallel_worker is FALSE in the session that runs the tests", {
  expect_false(in_parallel_worker())
})


test_that("in_parallel_worker is TRUE in a socket worker", {
  cl <- parallel::makePSOCKcluster(1)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # A copy outside the namespace, so the worker needs no installed HapSelect to run it
  check <- in_parallel_worker
  environment(check) <- globalenv()

  expect_true(parallel::clusterCall(cl, check)[[1]])
})


test_that("the whole banner is one startup message, so it can be suppressed", {
  # Before, the cli half was printed directly and survived suppressPackageStartupMessages()
  expect_silent(suppressPackageStartupMessages(HapSelect:::.onAttach("", "HapSelect")))

  expect_message(HapSelect:::.onAttach("", "HapSelect"), "Welcome to HapSelect", class = "packageStartupMessage")
  expect_output(suppressMessages(HapSelect:::.onAttach("", "HapSelect")), NA)
})
