# Helpers for tests that read machine readable progress events, shared by
# test-progress.R, test-graph-blocking.R and test-plink-ld.R

# Run body with machine logging on, then put every global it touched back. Returns the
# parsed events.
with_machine_logging <- function(body) {
  old_env <- Sys.getenv("HAPSELECT_MACHINE_LOGGING", unset = NA)
  old_enable <- getOption("progressr.enable")
  old_handlers <- progressr::handlers()
  on.exit({
    if (is.na(old_env)) {
      Sys.unsetenv("HAPSELECT_MACHINE_LOGGING")
    } else {
      Sys.setenv(HAPSELECT_MACHINE_LOGGING = old_env)
    }
    options(progressr.enable = old_enable)
    progressr::handlers(old_handlers)
    .progress_state$stack <- NULL
  }, add = TRUE)

  # The sequence a spawned process goes through: the platform sets the variable, then
  # .onLoad() switches progressr on and installs the handler.
  Sys.setenv(HAPSELECT_MACHINE_LOGGING = "1")
  options(progressr.handlers = NULL)
  options(progressr.enable = machine_logging() || interactive() || isatty(stderr()))
  register_handler()

  parse_events(capture.output(type = "message", invisible(body())))
}

# The reader a platform would write: take the marked lines, parse the rest as JSON, leave
# everything else alone as ordinary log output.
parse_events <- function(lines) {
  marked <- lines[startsWith(lines, "@@HS ")]
  lapply(marked, function(line) {
    jsonlite::fromJSON(substring(line, nchar("@@HS ") + 1L), simplifyVector = FALSE)
  })
}

kinds <- function(events) vapply(events, function(e) e$event, character(1))
steps <- function(events) vapply(events, function(e) e$step %||% NA_character_, character(1))
`%||%` <- function(x, y) if (is.null(x)) y else x

# Run body and count the progress ticks it makes. Counted from the progression
# conditions themselves, since the events cannot show it: progressr fills the count up
# to the total when a progressor finishes, ticked or not.
count_ticks <- function(body) {
  n <- 0L
  withCallingHandlers(body(), progression = function(c) {
    if (identical(c$type, "update")) n <<- n + 1L
  })
  n
}
