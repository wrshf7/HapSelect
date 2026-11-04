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

  result <- def_blocks(strategy, map)

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

  def_blocks(ld_strategy("fake-ld"), make_dispatch_map_fixture())

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

  result <- def_blocks(strategy, map)

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

  def_blocks(window_strategy(window = 3, method = "window_snp"), make_dispatch_map_fixture())

  expect_false(ld_called)
  expect_true(window_called)
})


# graph_strategy dispatch --------------------------------------------------------

test_that("def_blocks dispatches a graph_strategy to perform_graph_blocking", {
  captured <- NULL

  local_mocked_bindings(
    perform_graph_blocking = function(...) {
      captured <<- list(...)
      "GRAPH_BACKEND_RESULT"
    }
  )

  map      <- make_dispatch_map_fixture()
  geno     <- data.frame(SNP = c("m1", "m2"), Chromosome = 1, Position = c(100, 200),
                         Ind1 = c(0, 1), Ind2 = c(2, 1), stringsAsFactors = FALSE)
  strategy <- graph_strategy(theta_core = 0.6, window_ld = 5, window_core = 5)

  result <- def_blocks(strategy, map, geno = geno)

  expect_identical(result, "GRAPH_BACKEND_RESULT")
  expect_identical(captured$map, map)
  expect_identical(captured$geno, geno)

  # The graph backend takes the whole strategy rather than unpacked arguments, so
  # a parameter added to graph_strategy() reaches it without touching def_blocks()
  expect_identical(captured$strategy, strategy)
  expect_equal(captured$strategy$theta_core, 0.6)
  expect_equal(captured$strategy$window_ld, 5)
  expect_equal(captured$strategy$window_core, 5)
})


test_that("def_blocks calls only the graph backend for a graph_strategy", {
  graph_called  <- FALSE
  ld_called     <- FALSE
  window_called <- FALSE

  local_mocked_bindings(
    perform_graph_blocking  = function(...) { graph_called  <<- TRUE; list() },
    perform_ld_blocking     = function(...) { ld_called     <<- TRUE; list() },
    perform_window_blocking = function(...) { window_called <<- TRUE; list() }
  )

  geno <- data.frame(SNP = c("m1", "m2"), Chromosome = 1, Position = c(100, 200),
                     Ind1 = c(0, 1), Ind2 = c(2, 1), stringsAsFactors = FALSE)
  def_blocks(graph_strategy(), make_dispatch_map_fixture(), geno = geno)

  expect_true(graph_called)
  expect_false(ld_called)
  expect_false(window_called)
})

test_that("def_blocks warns when geno is passed to a strategy that ignores it", {
  # Accepting genotypes a strategy never reads would hide the mistake until the
  # blocks came back with no relation to the data that was handed over
  local_mocked_bindings(
    perform_window_blocking = function(...) list(),
    perform_ld_blocking     = function(...) list()
  )

  geno <- data.frame(SNP = c("m1", "m2"), Chromosome = 1, Position = c(100, 200),
                     Ind1 = c(0, 1), Ind2 = c(2, 1), stringsAsFactors = FALSE)
  map  <- make_dispatch_map_fixture()

  expect_warning(def_blocks(window_strategy(window = 2), map, geno = geno),
                 "only used by graph_strategy")
  expect_warning(def_blocks(ld_strategy("fake-ld"), map, geno = geno),
                 "only used by graph_strategy")

  # and says nothing when it is left out
  expect_silent(def_blocks(window_strategy(window = 2), map))
})
# Invalid / unrecognised strategies ----------------------------------------------

test_that("def_blocks rejects a strategy not built by a strategy constructor", {
  map <- make_dispatch_map_fixture()

  expect_error(def_blocks(list(method = "flanking"), map), "strategy constructor")
  expect_error(def_blocks("ld_strategy", map),             "strategy constructor")
  expect_error(def_blocks(NULL, map),                      "strategy constructor")
})


test_that("def_blocks errors clearly on a strategy subclass with no matching backend", {
  map              <- make_dispatch_map_fixture()
  mystery_strategy <- structure(list(), class = c("mystery_strategy", "block_strategy"))

  expect_error(def_blocks(mystery_strategy, map), "mystery_strategy")
})
