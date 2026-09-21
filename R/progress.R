# R/progress.R
#
# Progress reporting. A progressr handler only sees counts, never which analysis is counting,
# so report_step() names the step and both backends read the name from there.
#
# HAPSELECT_MACHINE_LOGGING picks the backend: unset draws a cli bar per step, "1" writes
# one JSON line per event to stderr. Steps nest, so a reader needs a stack, not a slot.
#
#   @@HS {"v":1,"at":"2026-09-10T04:21:09.114Z","event":"step_started","step":"ld_blocking"}
#
#   step_started   step | progress  step, current, total | step_finished  step, status


# --- the steps this package reports ---------------------------------------------------------

##### Every step HapSelect can report, and the name to show a person #####
# The ids go on the wire and are keyed on by whatever is reading, so treat them as interface
HS_STEPS = c(
  pairwise_ld              = "Pairwise LD",
  plink_ld                 = "PLINK Pairwise LD",
  plink_ld_geno            = "PLINK Pairwise LD",
  ld_blocking              = "LD Blocking",
  window_blocking          = "Window Blocking",
  marker_effects           = "Solving Marker Effects",
  cross_validation         = "Cross Validation",
  n_fold_cross_validation  = "N-Fold Cross Validation",
  gebv_chunk_prep          = "Preparing Chunks",
  local_gebv               = "Local GEBV",
  haplotype_effects        = "Haplotype Effects",
  block_var_test           = "Haploblock Variance Test",
  haplotype_simulation     = "Haplotype vs Truncation Simulation",
  local_gebv_simulation    = "Local GEBV vs Truncation Simulation"
)

# Machine readable progress logging #

# Return true if machine reporting is enabled
# Machine reporting is used for consuming platforms
machine_logging = function() {
  Sys.getenv("HAPSELECT_MACHINE_LOGGING", unset = "") == "1"
}

# Holds the current progress state
# This is held in a env var so it persists across package calls
.progress_state = new.env(parent = emptyenv())

# Returns the current runnings step at the top of the stack
current_step = function() {
  stack = .progress_state$stack
  if (is.null(stack) || length(stack) == 0L) return(NULL)
  stack[[length(stack)]]
}

# A prefix placed at the start of each machine readable progress state 
# Separates legitimate progress logs from warning and other tool output sharing stderr
HS_PREFIX = "@@HS "
# The version of the machine readable logging format
PROTOCOL_VERSION = 1L

# Create an event timestamp in ISO8601 format for parsing by consumers
event_timestamp = function(time = Sys.time()) {
  format(as.POSIXct(time), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

# Emit an event log with the event, step, current step, total steps and status
emit_event = function(event, step = NULL, current = NULL, total = NULL, status = NULL, message = NULL) {
  # Create the list of properties
  fields = list(step = step, current = current, total = total, status = status, message = message)
  # Include only populated fields
  valid_fields = fields[!vapply(fields, is.null, logical(1))]
  # Add the version and timestamp
  record = c(list(v = 1L, at = event_timestamp(), event = event), valid_fields)

  # Add prefix,  jsonify the record and output
  cat(HS_PREFIX,
      jsonlite::toJSON(record, auto_unbox = TRUE, na = "null", digits = 6),
      "\n", sep = "", file = stderr())
}

# Human readable progress logging #

# Returns a human readable label for regular progress logging, or the raw id if it has none
step_label = function(id) {
  id = as.character(id)
  if (id %in% names(HS_STEPS)) unname(HS_STEPS[[id]]) else id
}

# Starts a progress bar for a step
bar_start = function(step) {
  # Do not delay showing progress bars
  shown = options(cli.progress_show_after = getOption("cli.progress_show_after", 0))
  # When the function exits, restore defaults
  on.exit(options(shown), add = TRUE)

  # Create the cli bar
  bar = try_cli(cli::cli_progress_bar(
    step_label(step),
    total          = NA,                                       # unknown total draws a spinner, so uncountable work still shows as running
    clear          = getOption("cli.progress_clear", FALSE),   # leave the finished bar on screen
    auto_terminate = FALSE,                                    # would close the bar the moment the count reaches the total
    .auto_close    = FALSE,                                    # would close the bar when this function returns
    .envir         = globalenv()
  ))
  
  # Associate the progress bar id with the step
  if (!is.null(bar)) .progress_state$bars[[step]] = bar
}

# Attempt to draw a CLI bar, catching the error in case something goes wrong
# A failed progress bar is better than a failed analysis
try_cli = function(expr) {
  tryCatch(expr, error = function(e) NULL)
}

# Update a progress bars state
bar_update = function(step, current, total) {
  # Get the progress bar id
  bar = .progress_state$bars[[step]]
  # If its null, return
  if (is.null(bar)) return(invisible(FALSE))

  # Update the progress bar
  try_cli(cli::cli_progress_update(
    id     = bar,
    set    = if (is.na(current)) NULL else current,
    total  = if (is.na(total)) NULL else total,
    .envir = globalenv()
  ))
}

# Mark a progress bar as complete
bar_done = function(step) {
  # Get the progress bar id
  bar = .progress_state$bars[[step]]
  # If its null, return
  if (is.null(bar)) return(invisible(FALSE))

  # Mark the progress bar as done
  try_cli(cli::cli_progress_done(id = bar))
  # Clear the progress bar id
  .progress_state$bars[[step]] = NULL
}

# Shared step reporting #

# Report a step to the progress logging
# Called once at the top of an analysis function
# The step automatically closes itself when that function returns
#
# envir is the calling functions name
report_step = function(id, envir = parent.frame()) {
  # Add the blank progress item to the stack
  id = as.character(id)
  .progress_state$stack = c(.progress_state$stack, id)

  # Log as either a machine readable event, or start a progress bar
  if (machine_logging()) {
    emit_event("step_started", step = id)
  } else {
    bar_start(id)
  }

  # Register a call to close the step when it finishes
  # Automatically runs when the function returns
  finish = bquote(finish_step(.(id), returnValue(default = quote(.step_failed))))
  thunk  = as.call(list(function() eval(finish, envir)))
  do.call(base::on.exit, list(thunk, TRUE, TRUE), envir = envir)
}

# Close a step once its finished given an id
finish_step = function(id, returned) {
  # Check if the step failed
  failed = identical(returned, quote(.step_failed))
  status = if (failed) "failed" else "succeeded"

  # Get the error message if the step failed
  message = if (failed) sub("^Error[^:]*: ", "", trimws(geterrmessage())) else NULL

  # Remove the step from the stack
  stack = .progress_state$stack
  hit = which(stack == id)
  if (length(hit) > 0L) .progress_state$stack = stack[-hit[[length(hit)]]]

  # Either log the step as finished in machine code (machine readable) or just close the bar (human readable)
  if (machine_logging()) {
    emit_event("step_finished", step = id, status = status, message = message)
  } else {
    bar_done(id)
  }
}

# Helper functions #
# Coerce a progressr count to a whole number, or NA if it is not usable 
as_count = function(value) {
  count = suppressWarnings(as.integer(value))
  if (length(count) != 1L || is.na(count) || count < 0L) NA_integer_ else count
}

# Report a tick or progress update
report_tick = function(config, state) {
  # Get the current step
  id = current_step()
  if (is.null(id)) return(invisible(FALSE))

  # Get the total and current step counters
  total = as_count(config$max_steps)
  current = as_count(state$step)
  
  # Check the total and valid counts are valid
  if (!is.na(total) && total <= 0L) total = NA_integer_
  if (!is.na(total) && !is.na(current)) current = min(current, total)

  # Report a tick either in machine readable format, or update the progress bar
  if (machine_logging()) {
    emit_event("progress", step = id, current = current, total = total)
  } else {
    bar_update(id, current, total)
  }
}

# Provides callbacks for progressr including initiation, updates and finish 
progress_reporter = function() {
  report = function(config, state, ...) report_tick(config, state)
  list(initiate = report, update = report, finish = report)
}

# progressr handler that reports HapSelect's named steps 
# Usses progress_reporteras the core function
handler_hapselect = function(...) {
  progressr::make_progression_handler(
    "hapselect", progress_reporter(),
    # Throttle progress updates so we arent flooded with events, 200ms
    interval = getOption("HapSelect.progress_interval", 0.2), ...
  )
}

# Install HapSelect's progress handler, unless the caller has chosen their own
register_handler = function() {
  # If the consumer already has a handler defined, let it be
  if (!is.null(getOption("progressr.handlers"))) return(invisible(FALSE))

  # Otherwise register the handler
  progressr::handlers(handler_hapselect())
  invisible(TRUE)
}
