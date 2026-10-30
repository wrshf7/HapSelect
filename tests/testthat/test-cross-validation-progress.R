# Tests: cross-validation progress reporting -----------------------------------
#
# Each fold is a tick of the cross-validation's own step, rather than a marker effects step of
# its own, so a reader sees one bar across the folds.

cv_fixture <- function(n = 40L, markers = 60L) {
  set.seed(3)
  g <- matrix(sample(0:2, n * markers, replace = TRUE), markers,
              dimnames = list(NULL, paste0("i", seq_len(n))))
  list(
    geno = data.frame(SNP = paste0("s", seq_len(markers)), Chromosome = 1L,
                      Position = seq_len(markers), g, stringsAsFactors = FALSE),
    BLUE = data.frame(id = paste0("i", seq_len(n)), y = rnorm(n), stringsAsFactors = FALSE)
  )
}

expect_one_bar_over_folds <- function(events, step, folds) {
  ticks <- Filter(function(e) identical(e$event, "progress"), events)

  expect_equal(unique(steps(events)), step)
  expect_equal(kinds(events)[c(1, length(events))], c("step_started", "step_finished"))
  expect_equal(ticks[[length(ticks)]]$current, folds)
  expect_equal(ticks[[length(ticks)]]$total, folds)
}

test_that("n_fold_cross_validation reports one bar over its folds, with no step per fold", {
  d <- cv_fixture()
  events <- with_machine_logging(function() {
    capture.output(suppressWarnings(n_fold_cross_validation(d$geno, d$BLUE, nfold = 4L)))
  })

  expect_one_bar_over_folds(events, "n_fold_cross_validation", 4)
})

test_that("cross_validation reports one bar over its folds, with no step per fold", {
  d <- cv_fixture()
  events <- with_machine_logging(function() {
    capture.output(suppressWarnings(cross_validation(d$geno, d$BLUE, fold = 3L)))
  })

  expect_one_bar_over_folds(events, "cross_validation", 3)
})

test_that("each fold's tick names the fold", {
  d <- cv_fixture()
  messages <- character(0)
  with_machine_logging(function() {
    withCallingHandlers(
      capture.output(suppressWarnings(n_fold_cross_validation(d$geno, d$BLUE, nfold = 3L))),
      progression = function(c) if (identical(c$type, "update")) messages <<- c(messages, c$message)
    )
  })

  expect_setequal(messages, c("Fold 1/3", "Fold 2/3", "Fold 3/3"))
})

test_that("solve_marker_effects still reports its own step when called directly", {
  d <- cv_fixture()
  events <- with_machine_logging(function() {
    capture.output(solve_marker_effects(d$geno, d$BLUE, "VanRaden", 2L, TRUE))
  })

  expect_equal(unique(steps(events)), "marker_effects")
})

test_that("each fold's accuracy, the mean and every fold's variance components are reported as metrics", {
  d <- cv_fixture()
  cv <- NULL
  events <- with_machine_logging(function() {
    capture.output(cv <<- suppressWarnings(n_fold_cross_validation(d$geno, d$BLUE, nfold = 3L)))
  })
  metrics <- Filter(function(e) identical(e$event, "metric"), events)
  names <- vapply(metrics, function(e) e$name, character(1))

  expect_equal(sum(names == "variance_components"), 3)
  folds <- metrics[names == "cv_fold"]
  expect_equal(vapply(folds, function(e) e$values$fold, numeric(1)), 1:3)
  expect_equal(vapply(folds, function(e) e$values$accuracy, numeric(1)), as.numeric(cv))
  summary <- metrics[names == "cv_summary"][[1]]$values
  expect_equal(summary$mean_accuracy, mean(cv))
  for (m in metrics) expect_equal(m$step, "n_fold_cross_validation")
})
