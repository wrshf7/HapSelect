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


# Run body with a bar drawn, returning each line cli draws. cli is made not dynamic, so it
# sends each draw as a whole line, by signalling a cliMessage that testthat would otherwise
# take for itself. Ticks are sent straight to report_tick(), as progressr would send them
draw_lines <- function(body) {
  old_env <- Sys.getenv("HAPSELECT_MACHINE_LOGGING", unset = NA)
  old <- options(progressr.enable = TRUE, cli.dynamic = FALSE, cli.ansi = FALSE,
                 cli.progress_show_after = 0, cli.width = 200)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")
    else Sys.setenv(HAPSELECT_MACHINE_LOGGING = old_env)
    options(old)
    .progress_state$stack <- NULL
    .progress_state$bar <- NULL
  }, add = TRUE)
  Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")

  lines <- character()
  withCallingHandlers(body(), cliMessage = function(cnd) {
    lines <<- c(lines, trimws(conditionMessage(cnd)))
  })
  lines[nzchar(lines)]
}

tick <- function(current, total, message = NULL) {
  report_tick(list(max_steps = total), list(step = current, message = message))
}


test_that("a step inside another shows in the bar's status, rather than as a bar of its own", {
  # cli draws one bar at a time, so a bar for the inner step drew over the outer step's, and
  # the finished bars were left on screen in the wrong order
  inner <- function() { report_step("read_vcf"); tick(3, 5); report_status("Converting genotypes") }
  outer <- function() {
    report_step("beagle_impute")
    report_status("Running Beagle")
    tick(4, 4)
    inner()
  }

  lines <- draw_lines(outer)

  expect_false(any(grepl("Reading VCF [^0-9]", lines) & !grepl("^Imputation", lines)))
  expect_true(any(grepl("^Imputation \\(Beagle\\) .*100% \\| Reading VCF 60% \\| Converting genotypes \\|", lines)))
  # The outer bar's count is its own, never the inner step's
  expect_false(any(grepl("^Imputation \\(Beagle\\) [^%]* 60% \\|", lines)))
  # Once the inner step is done, the bar says neither what it nor the outer step was doing
  last <- utils::tail(lines, 1)
  expect_match(last, "^Imputation \\(Beagle\\) .*100%")
  expect_no_match(last, "Reading VCF|Running Beagle")
  expect_length(.progress_state$stack, 0)
})


test_that("a step inside another can share its id, and each unwinds on its own", {
  inner <- function() { report_step("read_vcf"); tick(1, 2) }
  outer <- function() { report_step("read_vcf"); tick(1, 4); inner(); tick(4, 4) }

  lines <- draw_lines(outer)

  expect_true(any(grepl("^Reading VCF .*25% \\| Reading VCF 50% \\|", lines)))
  expect_match(utils::tail(lines, 1), "^Reading VCF .*100%")
  expect_null(.progress_state$bar)
})


test_that("a failed step's bar says it failed, rather than looking like it is still running", {
  # cli only uses format_failed when one is given, and otherwise leaves the running format
  run <- function() { report_step("local_gebv"); tick(2, 4, "Fold 2/4"); stop("no markers") }

  lines <- draw_lines(function() try(run(), silent = TRUE))

  expect_match(utils::tail(lines, 1), "^\u2716 Local GEBV .* 50% \\| Fold 2/4 \\| Failed after ")
})


test_that("a failed inner step leaves its status on the bar, which fails with it", {
  inner <- function() { report_step("read_vcf"); tick(3, 5); stop("bad VCF") }
  outer <- function() { report_step("beagle_impute"); inner() }

  lines <- draw_lines(function() try(outer(), silent = TRUE))

  # The outer step never counted, so its failed line has no bar
  expect_match(utils::tail(lines, 1),
               "^\u2716 Imputation \\(Beagle\\) \\| Reading VCF 60% \\| Failed after ")
  expect_length(.progress_state$stack, 0)
})


test_that("an inner step's status names it, with how far through it is", {
  expect_equal(nested_status(list(id = "read_vcf")), "Reading VCF")
  expect_equal(nested_status(list(id = "read_vcf", current = 1L, total = 3L)), "Reading VCF 33%")
  expect_equal(nested_status(list(id = "read_vcf", current = 3L, total = 3L,
                                  status = "Converting genotypes")),
               "Reading VCF 100% | Converting genotypes")
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

  closed <- logical()
  local_mocked_bindings(bar_done = function(failed = FALSE) {
    closed <<- c(closed, failed)
  })

  ok   <- function() { report_step("pairwise_ld"); invisible(NULL) }
  fail <- function() { report_step("local_gebv"); stop("no markers") }
  ok()
  expect_error(fail(), "no markers")

  expect_equal(closed, c(FALSE, TRUE))
})

test_that("the reporter can hide its bar, so held back output is printed above it", {
  # Without hide and unhide, progressr holds every cat() and message() inside with_progress()
  # until it exits, then prints them from the start of the bar's line, over the bar
  reporter <- progress_reporter()
  expect_true(all(c("hide", "unhide") %in% names(reporter)))

  # With no bar drawn there is nothing to clear or redraw
  old <- options(progressr.enable = FALSE)
  on.exit(options(old), add = TRUE)
  expect_false(reporter$hide())
  expect_false(reporter$unhide())
})


test_that("the bar is only cleared when there is held back output to print", {
  # with_progress() held back its own ticks too, so it cleared and redrew the bar on every tick,
  # throttled or not, which drew as the bar flashing
  old_env <- Sys.getenv("HAPSELECT_MACHINE_LOGGING", unset = NA)
  old_handlers <- progressr::handlers()
  old <- options(progressr.enable = TRUE)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")
    else Sys.setenv(HAPSELECT_MACHINE_LOGGING = old_env)
    progressr::handlers(old_handlers)
    options(old)
  }, add = TRUE)
  Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")
  progressr::handlers(handler_hapselect())

  hides <- 0L
  local_mocked_bindings(bar_hide = function() { hides <<- hides + 1L; invisible(TRUE) },
                        bar_unhide = function() invisible(TRUE))

  messages <- character()
  withCallingHandlers(
    hs_with_progress({
      p <- progressr::progressor(steps = 100)
      for (i in 1:100) {
        if (i %% 50 == 0) message("printed ", i)
        p()
      }
    }),
    message = function(m) {
      messages <<- c(messages, trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    }
  )

  # Messages are still held back and printed, but the bar is cleared only to print them
  expect_equal(messages, c("printed 50", "printed 100"))
  expect_lte(hides, 2L)
})
