# R/progress.R
#
# Progress reporting. A progressr handler only sees counts, never which analysis is counting,
# so report_step() names the step and both backends read the name from there.
#
# HAPSELECT_MACHINE_LOGGING picks the backend: unset draws a cli bar per step, "1" writes
# one JSON line per event to stderr. Steps nest, so a reader needs a stack, not a slot.
#
#   @@HS {"v":"2.0.0","at":"2026-09-10T04:21:09.114Z","event":"step_started","step":"ld_blocking","label":"LD Blocking"}
#
#   step_started   step, label
#   progress       step, current, total, and message when the tick carried one. A status sent
#                  during uncounted work repeats the step's last counts, and has none before
#                  the step has counted anything
#   step_finished  step, status, and message when it failed
#   metric         step, name, values: a named set of numbers a step produced, such as a GA
#                  generation's fitness or a model's variance components


# --- the steps this package reports ---------------------------------------------------------

##### Every step HapSelect can report, and the name to show a person #####
# The ids go on the wire and are keyed on by whatever is reading, so treat them as interface.
# The names are for people: say what the step does, with an external tool it runs in brackets
HS_STEPS = c(
  pairwise_ld                 = "Pairwise LD",
  plink_ld                    = "Pairwise LD (PLINK)",
  beagle_impute               = "Imputation (Beagle)",
  beagle_phase                = "Phasing (Beagle)",
  ld_blocking                 = "LD Blocking",
  window_blocking             = "Window Blocking",
  graph_blocking              = "Graph Blocking",
  marker_effects              = "Marker Effects",
  cross_validation            = "Cross Validation",
  n_fold_cross_validation     = "N-Fold Cross Validation",
  gebv_chunk_prep             = "Preparing Blocks",
  local_gebv                  = "Local GEBV",
  haplotype_effects           = "Haplotype Effects",
  block_var_test              = "Haploblock Variance Test",
  haplotype_simulation        = "Haplotype vs Truncation Simulation",
  local_gebv_simulation       = "Local GEBV vs Truncation Simulation",
  local_gebv_parent_selection = "Local GEBV Parent Selection",
  haplotype_parent_selection  = "Haplotype Parent Selection",
  read_vcf                    = "Reading VCF"
)

# Machine readable progress logging #

# Return true if machine reporting is enabled
# Machine reporting is used for consuming platforms
machine_logging = function() {
  Sys.getenv("HAPSELECT_MACHINE_LOGGING", unset = "") == "1"
}

# Return true if a person is watching a progress bar, so a line of text can ride on the bar
# rather than be printed. When it is false the text should be printed: with machine logging
# nothing is drawn, and in a batch run progressr is switched off, so the bar would drop it
progress_drawn = function() {
  !machine_logging() && isTRUE(getOption("progressr.enable", interactive()))
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
# The version of the machine readable logging format, which is the package's own version.
# Read when the first event is sent rather than here, where the package is still being built
protocol_version = function() {
  if (is.null(.progress_state$version)) {
    .progress_state$version = as.character(utils::packageVersion("HapSelect"))
  }
  .progress_state$version
}

# Create an event timestamp in ISO8601 format for parsing by consumers
event_timestamp = function(time = Sys.time()) {
  format(as.POSIXct(time), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

# Emit an event log with the event and its named fields, such as step, current, total and status
emit_event = function(event, ...) {
  # Include only populated fields
  fields = list(...)
  valid_fields = fields[!vapply(fields, is.null, logical(1))]
  # Add the version and timestamp
  record = c(list(v = protocol_version(), at = event_timestamp(), event = event), valid_fields)

  # Add prefix,  jsonify the record and output
  # digits = NA writes numbers in full, so a metric as small as a marker variance is not rounded away
  cat(HS_PREFIX,
      jsonlite::toJSON(record, auto_unbox = TRUE, na = "null", digits = NA),
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

  # cli keeps one live bar per environment and closes it when another is started there, so
  # each step's bar gets an environment of its own. Otherwise a nested step's bar would close
  # the bar of the step around it, which would then never draw again
  envir = new.env(parent = globalenv())

  # Create the cli bar
  bar = try_cli(cli::cli_progress_bar(
    step_label(step),
    total          = NA,                                       # unknown total draws a spinner, so uncountable work still shows as running
    clear          = getOption("cli.progress_clear", FALSE),   # leave the finished bar on screen
    auto_terminate = FALSE,                                    # would close the bar the moment the count reaches the total
    .auto_close    = FALSE,                                    # would close the bar when this function returns
    .envir         = envir
  ))

  # Associate the progress bar id, and the environment it lives in, with the step
  if (!is.null(bar)) .progress_state$bars[[step]] = list(id = bar, envir = envir)

  # cli only draws on an update, so a step that never ticks would otherwise never be seen.
  # Only when progress is shown: with progressr off, cli would print the bar as a message
  if (!is.null(bar) && progress_drawn()) {
    try_cli(cli::cli_progress_update(id = bar, set = 0, force = TRUE, .envir = envir))
  }
}

# Attempt to draw a CLI bar, catching the error in case something goes wrong
# A failed progress bar is better than a failed analysis
# cli settles a bar's format when it first draws it, or when its total changes. Until a step
# ticks it has no total, and cli's own format for that counts "1 done (3.1/s)", which means
# nothing for work that is not counted, so the spinner shows the time taken instead.
# A format the user has set is left alone
try_cli = function(expr) {
  spinner = options(cli.progress_format_iterator_nototal = getOption(
    "cli.progress_format_iterator_nototal",
    "{cli::pb_spin} {cli::pb_name}| {cli::pb_status}{cli::pb_elapsed}"
  ))
  on.exit(options(spinner), add = TRUE)
  tryCatch(expr, error = function(e) NULL)
}

# Update a progress bars state
# status is optional text shown beside the bar, taken from the message a progressor ticks with
bar_update = function(step, current, total, status = NULL) {
  # Get the progress bar
  bar = .progress_state$bars[[step]]
  # If its null, return
  if (is.null(bar)) return(invisible(FALSE))

  # Update the progress bar
  try_cli(cli::cli_progress_update(
    id     = bar$id,
    set    = if (is.na(current)) NULL else current,
    total  = if (is.na(total)) NULL else total,
    # cli's default format puts the status directly before the ETA, so it carries its own separator
    status = if (is.null(status)) NULL else paste0(status, " | "),
    # progressr already limits updates to one per interval, and cli's own limit on top would
    # drop a status set just before a long wait, such as starting parallel workers
    force  = TRUE,
    .envir = bar$envir
  ))
}

# Mark a progress bar as complete
# failed: TRUE when the step stopped with an error, so the bar shows it failed rather than finished
bar_done = function(step, failed = FALSE) {
  # Get the progress bar
  bar = .progress_state$bars[[step]]
  # If its null, return
  if (is.null(bar)) return(invisible(FALSE))

  # Mark the progress bar as done
  try_cli(cli::cli_progress_done(id = bar$id, .envir = bar$envir,
                                 result = if (failed) "failed" else "done"))
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
    emit_event("step_started", step = id, label = step_label(id))
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

  # Remove the step from the stack, and the counts it last reported
  stack = .progress_state$stack
  hit = which(stack == id)
  if (length(hit) > 0L) .progress_state$stack = stack[-hit[[length(hit)]]]
  .progress_state$counts[[id]] = NULL

  # Either log the step as finished in machine code (machine readable) or just close the bar (human readable)
  if (machine_logging()) {
    emit_event("step_finished", step = id, status = status, message = message)
  } else {
    bar_done(id, failed = failed)
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

  # Kept so that a status sent later, during work with no ticks, can repeat where the step got to
  .progress_state$counts[[id]] = list(current = current, total = total)

  # A tick's message, such as the fold or rep being run, goes beside the bar or on the event
  message = state$message
  status = if (is.character(message) && length(message) == 1L && nzchar(message)) message else NULL

  # Report a tick either in machine readable format, or update the progress bar
  if (machine_logging()) {
    emit_event("progress", step = id, current = current, total = total, message = status)
  } else {
    bar_update(id, current, total, status)
  }
}

# Say what the current step is doing during work that has no ticks of its own, such as
# writing a file or building a plot after the counted work is done
# Shown beside the bar, or sent as a progress event carrying the message
report_status = function(message) {
  id = current_step()
  if (is.null(id)) return(invisible(FALSE))

  if (machine_logging()) {
    # The step's last counts, so a reader does not take the status as the count starting over
    counts = .progress_state$counts[[id]]
    emit_event("progress", step = id, current = counts$current, total = counts$total,
               message = message)
  } else if (progress_drawn()) {
    # Only when a bar is drawn: with progressr off, cli would print the bar as a message
    bar_update(id, NA_integer_, NA_integer_, message)
  }
  invisible(TRUE)
}

# Report a named set of numbers a step produced, such as a GA generation's fitness
# Only machine logging carries these; a person reads the same numbers from the printed summary
#
# name: what the numbers are, e.g. "ga_generation"
# ... : the numbers, each named, e.g. iter = 3, best = 1.2
report_metric = function(name, ...) {
  if (!machine_logging()) return(invisible(FALSE))
  emit_event("metric", step = current_step(), name = name, values = list(...))
  invisible(TRUE)
}

# Provides callbacks for progressr including initiation, updates and finish 
progress_reporter = function() {
  report = function(config, state, ...) report_tick(config, state)
  list(initiate = report, update = report, finish = report)
}

# progressr handler that reports HapSelect's named steps 
# Uses progress_reporter as the core function
handler_hapselect = function(...) {
  progressr::make_progression_handler(
    "hapselect", progress_reporter(),
    # Throttle progress updates so we arent flooded with events, 200ms
    interval = getOption("HapSelect.progress_interval", 0.2),
    # A "terminal" handler makes with_progress() hold back stdout, messages and warnings until it
    # exits, so they cannot break up a bar. Machine logging draws no bar, so it lets them through
    target = if (machine_logging()) "file" else "terminal", ...
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
