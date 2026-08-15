# Tests: def_blocks strategy dispatch ------------------------------------------
#
# These tests check def_blocks() itself: given a strategy object, does it call
# the right backend with the right arguments, and only that backend? The
# backends' own blocking logic is covered separately (test-def-haploblocks.R
# for perform_ld_blocking, test-window-blocks.R for perform_window_blocking),
# so here perform_ld_blocking() and perform_window_blocking() are mocked out
# rather than run for real.

make_dispatch_map_fixture <- function() {
  data.frame(
    SNP        = c("m1", "m2"),
    Chromosome = c(1L, 1L),
    Position   = c(100, 200),
    stringsAsFactors = FALSE
  )
}


# ld_strategy dispatch ----------------------------------------------------------

test_that("def_blocks dispatches an ld_strategy to perform_ld_blocking", {
  captured <- NULL

  local_mocked_bindings(
    perform_ld_blocking = function(...) {
      captured <<- list(...)
      "LD_BACKEND_RESULT"
    }
  )

  map      <- make_dispatch_map_fixture()
  ld       <- "fake-ld"   # opaque placeholder; the mock never inspects it
  strategy <- ld_strategy(ld, method = "average", threshold = 0.6, tolerance = 2,
                          tol_reset = FALSE, start = "beginning", parallel = TRUE)

  result <- def_blocks(map, strategy)

  expect_identical(result, "LD_BACKEND_RESULT")
  expect_identical(captured$ld, ld)
  expect_identical(captured$map, map)
  expect_equal(captured$method, "average")
  expect_equal(captured$threshold, 0.6)
  expect_equal(captured$tolerance, 2)
  expect_equal(captured$tol_reset, FALSE)
  expect_equal(captured$start, "beginning")
  expect_equal(captured$parallel, TRUE)
})


test_that("def_blocks calls only the LD backend for an ld_strategy, never the window backend", {
  ld_called     <- FALSE
  window_called <- FALSE

  local_mocked_bindings(
    perform_ld_blocking     = function(...) { ld_called     <<- TRUE; list() },
    perform_window_blocking = function(...) { window_called <<- TRUE; list() }
  )

  def_blocks(make_dispatch_map_fixture(), ld_strategy("fake-ld"))

  expect_true(ld_called)
  expect_false(window_called)
})


# window_strategy dispatch -------------------------------------------------------

test_that("def_blocks dispatches a window_strategy to perform_window_blocking", {
  captured <- NULL

  local_mocked_bindings(
    perform_window_blocking = function(...) {
      captured <<- list(...)
      "WINDOW_BACKEND_RESULT"
    }
  )

  map      <- make_dispatch_map_fixture()
  strategy <- window_strategy(window = 5, method = "window_map")

  result <- def_blocks(map, strategy)

  expect_identical(result, "WINDOW_BACKEND_RESULT")
  expect_identical(captured$map, map)
  expect_equal(captured$window, 5)
  expect_equal(captured$method, "window_map")
})


test_that("def_blocks calls only the window backend for a window_strategy, never the LD backend", {
  ld_called     <- FALSE
  window_called <- FALSE

  local_mocked_bindings(
    perform_ld_blocking     = function(...) { ld_called     <<- TRUE; list() },
    perform_window_blocking = function(...) { window_called <<- TRUE; list() }
  )

  def_blocks(make_dispatch_map_fixture(), window_strategy(window = 3, method = "window_snp"))

  expect_false(ld_called)
  expect_true(window_called)
})


# Invalid / unrecognised strategies ----------------------------------------------

test_that("def_blocks rejects a strategy not built by a strategy constructor", {
  map <- make_dispatch_map_fixture()

  expect_error(def_blocks(map, list(method = "flanking")), "strategy constructor")
  expect_error(def_blocks(map, "ld_strategy"),             "strategy constructor")
  expect_error(def_blocks(map, NULL),                      "strategy constructor")
})


test_that("def_blocks errors clearly on a strategy subclass with no matching backend", {
  map              <- make_dispatch_map_fixture()
  mystery_strategy <- structure(list(), class = c("mystery_strategy", "block_strategy"))

  expect_error(def_blocks(map, mystery_strategy), "mystery_strategy")
})
