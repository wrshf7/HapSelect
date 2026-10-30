# Tests: named step progress reporting ------------------------------------------
#
# The contract these pin down:
#   - every step the package declares has a name, in both directions
#   - machine mode emits parseable lines that reflect what actually ran
#   - a step is always closed, including when the function it names fails
#   - human mode puts nothing machine readable on stderr


# The registry -----------------------------------------------------------------

test_that("every step declared in the package has a label", {
  # A typo would otherwise ship a step the platform cannot map and nobody can name.
  sources <- list.files(test_path("..", ".."), pattern = "[.]R$",
                        recursive = TRUE, full.names = TRUE)
  sources <- sources[grepl("/R/", sources, fixed = TRUE)]
  skip_if(length(sources) == 0, "package sources not reachable from the test directory")

  declared <- unlist(lapply(sources, function(f) {
    lines <- readLines(f, warn = FALSE)
    hits <- regmatches(lines, regexpr('report_step\\("[^"]+"\\)', lines))
    gsub('report_step\\("|"\\)', "", hits)
  }))
  declared <- setdiff(declared, "")
  expect_gt(length(declared), 0)

  for (id in unique(declared)) expect_true(id %in% names(HS_STEPS), info = id)
})


test_that("an unlabelled step falls back to its id rather than failing", {
  # A missing label must never take down the analysis it is reporting on; the registry test
  # above is what stops one shipping.
  expect_equal(step_label("not_a_real_step"), "not_a_real_step")
  expect_equal(step_label("pairwise_ld"), "Pairwise LD")
})


test_that("every label is a non-empty human readable name", {
  expect_true(all(nzchar(HS_STEPS)))
  expect_equal(anyDuplicated(names(HS_STEPS)), 0L)
})


# Machine mode -----------------------------------------------------------------

test_that("a step is opened and closed around the function that declares it", {
  run <- function() { report_step("pairwise_ld"); "result" }

  events <- with_machine_logging(function() expect_equal(run(), "result"))

  expect_equal(kinds(events), c("step_started", "step_finished"))
  expect_equal(steps(events), c("pairwise_ld", "pairwise_ld"))
  expect_equal(events[[2]]$status, "succeeded")
})


test_that("a step that fails is closed and reported as failed, with the reason", {
  run <- function() { report_step("local_gebv"); stop("no markers on chromosome 7") }

  events <- with_machine_logging(function() expect_error(run(), "no markers"))

  # Left open, this would read as a job stalled forever.
  expect_equal(kinds(events), c("step_started", "step_finished"))
  expect_equal(events[[2]]$status, "failed")
  # The reason travels on the event, so a reader never has to match it to a nearby log line.
  expect_match(events[[2]]$message, "no markers on chromosome 7")
})


test_that("a step that succeeds carries no message at all", {
  run <- function() { report_step("pairwise_ld"); "ok" }

  events <- with_machine_logging(run)

  expect_equal(events[[2]]$status, "succeeded")
  expect_false("message" %in% names(events[[2]]))
})


test_that("an error message with JSON metacharacters survives the wire", {
  nasty <- 'bad "arg" = c(1,2) or {x}'
  run <- function() { report_step("block_var_test"); stop(nasty) }

  events <- with_machine_logging(function() expect_error(run()))

  expect_match(events[[2]]$message, 'bad "arg" = c\\(1,2\\) or \\{x\\}')
})


test_that("steps nest, and unwind innermost first", {
  inner <- function() { report_step("gebv_chunk_prep"); invisible(NULL) }
  outer <- function() { report_step("local_gebv"); inner(); invisible(NULL) }

  events <- with_machine_logging(outer)

  expect_equal(steps(events),
               c("local_gebv", "gebv_chunk_prep", "gebv_chunk_prep", "local_gebv"))
  expect_equal(kinds(events),
               c("step_started", "step_started", "step_finished", "step_finished"))
})


test_that("progress ticks are attributed to the innermost running step", {
  run <- function() {
    report_step("ld_blocking")
    progressr::with_progress({
      p <- progressr::progressor(steps = 4)
      for (i in 1:4) p()
    })
  }

  events <- with_machine_logging(run)
  ticks <- Filter(function(e) identical(e$event, "progress"), events)

  expect_gt(length(ticks), 0)
  for (tick in ticks) expect_equal(tick$step, "ld_blocking")
  expect_equal(ticks[[length(ticks)]]$current, 4)
  expect_equal(ticks[[length(ticks)]]$total, 4)
})


test_that("ticks arriving outside any step are dropped rather than reported anonymously", {
  run <- function() {
    progressr::with_progress({
      p <- progressr::progressor(steps = 3)
      for (i in 1:3) p()
    })
  }

  events <- with_machine_logging(run)
  expect_false("progress" %in% kinds(events))
})


test_that("declaring a step leaves the caller's own on.exit alone", {
  ran <- FALSE
  run <- function() {
    on.exit(ran <<- TRUE, add = TRUE)
    report_step("block_var_test")
    invisible(NULL)
  }

  events <- with_machine_logging(run)

  expect_true(ran)
  expect_equal(kinds(events), c("step_started", "step_finished"))
})


test_that("every emitted line is marked, parseable and self describing", {
  run <- function() { report_step("marker_effects"); invisible(NULL) }

  events <- with_machine_logging(run)

  expect_gt(length(events), 0)
  for (event in events) {
    expect_equal(event$v, as.character(utils::packageVersion("HapSelect")))
    expect_match(event$at, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}")
    expect_true(nzchar(event$event))
  }
})


test_that("ordinary R output sharing stderr is not mistaken for progress", {
  # stderr also carries message(), warnings and anything a shelled out tool prints, so the
  # marker is the whole of what separates progress from logs.
  lines <- c(
    "loading genotypes",
    "Warning message: 3 markers dropped",
    '@@HS {"v":1,"event":"step_started","step":"ld_blocking"}',
    "PLINK v1.90b6.21 64-bit",
    "  @@HS not marked, because the marker must start the line"
  )

  events <- parse_events(lines)
  expect_equal(length(events), 1)
  expect_equal(events[[1]]$step, "ld_blocking")
})


# Human mode -------------------------------------------------------------------

test_that("with machine logging off nothing machine readable reaches stderr", {
  old_env <- Sys.getenv("HAPSELECT_MACHINE_LOGGING", unset = NA)
  on.exit({
    if (is.na(old_env)) {
      Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")
    } else {
      Sys.setenv(HAPSELECT_MACHINE_LOGGING = old_env)
    }
    .progress_state$stack <- NULL
  }, add = TRUE)
  Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")

  run <- function() { report_step("window_blocking"); invisible(NULL) }
  emitted <- capture.output(type = "message", invisible(run()))

  expect_false(machine_logging())
  expect_equal(length(parse_events(emitted)), 0)
})


test_that("the machine logging switch reads the environment", {
  old_env <- Sys.getenv("HAPSELECT_MACHINE_LOGGING", unset = NA)
  on.exit({
    if (is.na(old_env)) {
      Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")
    } else {
      Sys.setenv(HAPSELECT_MACHINE_LOGGING = old_env)
    }
  }, add = TRUE)

  Sys.setenv(HAPSELECT_MACHINE_LOGGING = "1");     expect_true(machine_logging())
  # Anything but exactly "1" is off, so a typo fails to the human readable default rather
  # than to some half configured state.
  Sys.setenv(HAPSELECT_MACHINE_LOGGING = "true");  expect_false(machine_logging())
  Sys.setenv(HAPSELECT_MACHINE_LOGGING = "TRUE");  expect_false(machine_logging())
  Sys.setenv(HAPSELECT_MACHINE_LOGGING = "0");     expect_false(machine_logging())
  Sys.setenv(HAPSELECT_MACHINE_LOGGING = "");      expect_false(machine_logging())
  Sys.unsetenv("HAPSELECT_MACHINE_LOGGING");       expect_false(machine_logging())
})


# Registration -----------------------------------------------------------------

test_that("a handler the caller chose is never overridden", {
  old_handlers <- progressr::handlers()
  on.exit(progressr::handlers(old_handlers), add = TRUE)

  progressr::handlers("txtprogressbar")
  chosen <- getOption("progressr.handlers")

  # Choosing a handler belongs to whoever called the package.
  expect_false(register_handler())
  expect_identical(getOption("progressr.handlers"), chosen)
})


test_that("the handler is installed when nothing is registered", {
  old_handlers <- progressr::handlers()
  on.exit(progressr::handlers(old_handlers), add = TRUE)

  options(progressr.handlers = NULL)
  expect_true(register_handler())
  expect_false(is.null(getOption("progressr.handlers")))
})


test_that("a nested step's bar leaves the bar of the step around it open", {
  # cli closes a bar when another is started in the same environment, which left the outer
  # step's bar closed for good once a nested step began
  on.exit({
    bar_done("gebv_chunk_prep")
    bar_done("local_gebv")
  }, add = TRUE)

  bar_start("local_gebv")
  outer <- .progress_state$bars$local_gebv$id
  bar_start("gebv_chunk_prep")

  expect_false(is.null(cli:::clienv$progress[[outer]]))
})


# What events carry ------------------------------------------------------------

test_that("a started step carries the name to show a person", {
  run <- function() { report_step("ld_blocking"); invisible(NULL) }

  events <- with_machine_logging(run)

  expect_equal(events[[1]]$label, "LD Blocking")
})


test_that("a tick's message travels on the progress event", {
  run <- function() {
    report_step("cross_validation")
    progressr::with_progress({
      p <- progressr::progressor(steps = 2)
      p(message = "Fold 1/2")
      p(message = "Fold 2/2")
    })
  }

  old <- options(HapSelect.progress_interval = 0)
  on.exit(options(old), add = TRUE)
  events <- with_machine_logging(run)
  ticks <- Filter(function(e) identical(e$event, "progress"), events)

  expect_true("Fold 2/2" %in% vapply(ticks, function(e) e$message %||% NA_character_, character(1)))
})


test_that("a metric carries its step, its name and every value in full", {
  run <- function() {
    report_step("marker_effects")
    report_metric("variance_components", marker_variance = 4.0552e-05, h2 = 0.3212345678)
  }

  events <- with_machine_logging(run)
  metric <- Filter(function(e) identical(e$event, "metric"), events)[[1]]

  expect_equal(metric$step, "marker_effects")
  expect_equal(metric$name, "variance_components")
  # Rounded to a few decimal places, a variance this small would read as zero
  expect_equal(metric$values$marker_variance, 4.0552e-05)
  expect_equal(metric$values$h2, 0.3212345678)
})


test_that("a metric is not reported when nobody is reading the machine stream", {
  old_env <- Sys.getenv("HAPSELECT_MACHINE_LOGGING", unset = NA)
  on.exit(if (is.na(old_env)) Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")
          else Sys.setenv(HAPSELECT_MACHINE_LOGGING = old_env), add = TRUE)
  Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")

  emitted <- capture.output(type = "message", reported <- report_metric("cv_fold", fold = 1))

  expect_false(reported)
  expect_length(emitted, 0)
})


test_that("a status before the step has counted anything carries only its message", {
  run <- function() { report_step("beagle_impute"); report_status("Writing VCF") }

  events <- with_machine_logging(run)
  status <- Filter(function(e) identical(e$event, "progress"), events)

  expect_length(status, 1)
  expect_equal(status[[1]]$step, "beagle_impute")
  expect_equal(status[[1]]$message, "Writing VCF")
  expect_null(status[[1]]$current)
})


test_that("a status after the step has counted repeats its last counts", {
  run <- function() {
    report_step("read_vcf")
    progressr::with_progress({
      p <- progressr::progressor(steps = 4)
      for (i in 1:4) p()
    })
    report_status("Converting genotypes")
  }

  events <- with_machine_logging(run)
  status <- Filter(function(e) identical(e$message, "Converting genotypes"), events)[[1]]

  expect_equal(status$current, 4)
  expect_equal(status$total, 4)
})


test_that("a status outside any step is dropped", {
  events <- with_machine_logging(function() report_status("Writing VCF"))

  expect_length(events, 0)
})


test_that("a failed step's bar is closed as failed, a finished one as done", {
  old_env <- Sys.getenv("HAPSELECT_MACHINE_LOGGING", unset = NA)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")
    else Sys.setenv(HAPSELECT_MACHINE_LOGGING = old_env)
    .progress_state$stack <- NULL
  }, add = TRUE)
  Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")

  closed <- list()
  local_mocked_bindings(bar_done = function(step, failed = FALSE) {
    closed[[step]] <<- failed
  })

  ok   <- function() { report_step("pairwise_ld"); invisible(NULL) }
  fail <- function() { report_step("local_gebv"); stop("no markers") }
  ok()
  expect_error(fail(), "no markers")

  expect_false(closed$pairwise_ld)
  expect_true(closed$local_gebv)
})
