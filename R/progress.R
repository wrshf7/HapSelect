# R/progress.R
#
# Progress reporting. A progressr handler only sees counts, never which analysis is counting,
# so report_step() names the step and both backends read the name from there.
#
# HAPSELECT_MACHINE_LOGGING picks the backend: unset draws a cli bar, "1" writes one JSON
# line per event to stderr. Steps nest, so a reader needs a stack, not a slot.
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
# This is held in an environment so it persists across package calls
#   stack: the open steps, innermost last. Each is list(id, current, total, status): the step's
#          last count and status, which both backends read. current and total are NULL until
#          the step counts something
#   bar  : the one cli bar, list(id, envir), drawn for the step at the bottom of the stack
.progress_state = new.env(parent = emptyenv())

# The innermost open step's entry on the stack, or NULL outside any step
top_frame = function() {
  stack = .progress_state$stack
  if (length(stack) == 0L) NULL else stack[[length(stack)]]
}

# Returns the id of the innermost open step, or NULL outside any step
current_step = function() top_frame()$id

# Change fields of the innermost open step's entry, e.g. set_top(status = "Writing VCF").
# A field given as NULL is cleared
set_top = function(...) {
  n = length(.progress_state$stack)
  if (n == 0L) return(invisible(FALSE))
  fields = list(...)
  for (field in names(fields)) .progress_state$stack[[n]][field] = list(fields[[field]])
  invisible(TRUE)
}

# TRUE when a count is known: given, and usable
known = function(count) !is.null(count) && !is.na(count)

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
#
# cli draws one bar at a time, on the current line, so only the step at the bottom of the stack
# gets one. A step running inside it, such as reading Beagle's output back in, shows in that
# bar's status instead. A bar of its own would draw over the first, and the finished bars would
# be left on screen in the order they finished, the inner step's above the step it is part of
#
#   Imputation (Beagle) ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | Reading VCF 60% | Converting genotypes | ETA: 0s

# Returns a human readable label for regular progress logging, or the raw id if it has none
step_label = function(id) {
  id = as.character(id)
  if (id %in% names(HS_STEPS)) unname(HS_STEPS[[id]]) else id
}

# What a bar is left showing when its step fails. cli only uses this when the bar is closed with
# result = "failed", and falls back to the bar's running format without it, so a failed bar used
# to be left frozen where it stopped: a spinner, or a part filled bar, as if it were still going.
# The bar and percentage only show when the step got as far as counting, and the status says
# what the step was doing when it failed. cli's status token carries a trailing space, so it is
# trimmed and given one separator here
#
#   ✖ Imputation (Beagle) ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ 100% | Reading VCF 60% | Failed after 3s
#   ✖ Local GEBV | Preparing Blocks 40% | Failed after 210ms
BAR_FORMAT_FAILED = paste0(
  "{cli::col_red(cli::symbol$cross)} {cli::pb_name}{cli::pb_bar}",
  "{if (is.na(cli::pb_total)) '' else paste0(' ', trimws(cli::pb_percent), ' ')}| ",
  "{if (nzchar(trimws(cli::pb_status))) paste0(trimws(cli::pb_status), ' ') else ''}",
  "{cli::col_red('Failed')} after {cli::pb_elapsed}"
)

# Start the bar, for the step at the bottom of the stack
bar_start = function(step) {
  # Do not delay showing progress bars
  shown = options(cli.progress_show_after = getOption("cli.progress_show_after", 0))
  # When the function exits, restore defaults
  on.exit(options(shown), add = TRUE)

  # cli ties a bar to an environment. This one is never left, so the bar stays open until
  # bar_done() closes it
  envir = new.env(parent = globalenv())

  # Create the cli bar
  bar = try_cli(cli::cli_progress_bar(
    step_label(step),
    total          = NA,                                       # unknown total draws a spinner, so uncountable work still shows as running
    clear          = getOption("cli.progress_clear", FALSE),   # leave the finished bar on screen
    format_failed  = BAR_FORMAT_FAILED,                        # without it, a failed bar is left looking like it is still running
    auto_terminate = FALSE,                                    # would close the bar the moment the count reaches the total
    .auto_close    = FALSE,                                    # would close the bar when this function returns
    .envir         = envir
  ))

  .progress_state$bar = if (is.null(bar)) NULL else list(id = bar, envir = envir)
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

# The status a step running inside the bar's own step shows beside the bar: its name, how far
# through it is, and its own status
#
#   Reading VCF 60%
#   Reading VCF 100% | Converting genotypes
nested_status = function(frame) {
  shown = step_label(frame$id)
  if (known(frame$current) && known(frame$total) && frame$total > 0L) {
    shown = sprintf("%s %d%%", shown, as.integer(floor(100 * frame$current / frame$total)))
  }
  if (!is.null(frame$status)) shown = paste(shown, frame$status, sep = " | ")
  shown
}

# Draw the bar as the stack stands. The count is the bottom step's, which owns the bar, and the
# status is the innermost step's
# Only when progress is shown: with progressr off, cli would print the bar as a message
bar_render = function() {
  bar = .progress_state$bar
  stack = .progress_state$stack
  if (is.null(bar) || length(stack) == 0L || !progress_drawn()) return(invisible(FALSE))

  owner = stack[[1]]
  top = stack[[length(stack)]]
  status = if (length(stack) == 1L) top$status else nested_status(top)
  counted = known(owner$current)

  # cli's default format puts the status directly before the ETA, so it carries its own separator
  shown = if (is.null(status)) NULL else paste0(status, " | ")
  # cli keeps the last status it was given when given none, so a status no longer wanted is
  # cleared with an empty one. Only then: cli pads any status it was given, an empty one
  # included, so a bar that never had one would otherwise gain a space before its ETA
  cleared = is.null(shown) && !is.null(bar$shown)
  .progress_state$bar$shown = shown

  try_cli(cli::cli_progress_update(
    id     = bar$id,
    set    = if (counted) owner$current else NULL,
    # Without a count this only redraws. cli adds one when given neither
    inc    = if (counted) NULL else 0L,
    total  = if (known(owner$total)) owner$total else NULL,
    status = if (cleared) "" else shown,
    # progressr already limits updates to one per interval, and cli's own limit on top would
    # drop a status set just before a long wait, such as starting parallel workers
    force  = TRUE,
    .envir = bar$envir
  ))
  invisible(TRUE)
}

# Close the bar, once the step that owns it has finished
# failed: TRUE when the step stopped with an error, so the bar shows it failed rather than finished
bar_done = function(failed = FALSE) {
  bar = .progress_state$bar
  if (is.null(bar)) return(invisible(FALSE))
  .progress_state$bar = NULL
  try_cli(cli::cli_progress_done(id = bar$id, .envir = bar$envir,
                                 result = if (failed) "failed" else "done"))
  invisible(TRUE)
}

# Clear the bar off its line so text printed next starts on a clean line
# progressr calls this before printing output it held back during a with_progress() block,
# such as cat() output or a message, then calls bar_unhide(). Without it the held back text
# is printed from the start of the line the bar is drawn on, over the top of the bar
bar_hide = function() {
  if (!progress_drawn() || is.null(.progress_state$bar)) return(invisible(FALSE))

  # The connection cli draws its bars on, and the same way cli clears its own status line
  out = cli::cli_output_connection()
  if (cli::is_ansi_tty(out)) {
    cat("\r\033[K", file = out)
  } else if (cli::is_dynamic_tty(out)) {
    cat("\r", strrep(" ", cli::console_width()), "\r", sep = "", file = out)
  }
  invisible(TRUE)
}

# Draw the bar again, below the output progressr has just printed
bar_unhide = function() {
  bar = .progress_state$bar
  if (is.null(bar) || !progress_drawn()) return(invisible(FALSE))

  # inc = 0 redraws the bar as it stands, without counting anything
  try_cli(cli::cli_progress_update(id = bar$id, inc = 0L, force = TRUE, .envir = bar$envir))
  invisible(TRUE)
}

# Shared step reporting #

# Report a step to the progress logging
# Called once at the top of an analysis function
# The step automatically closes itself when that function returns
#
# envir is the calling functions name
report_step = function(id, envir = parent.frame()) {
  # Add the step to the stack, with nothing counted yet
  id = as.character(id)
  .progress_state$stack = c(.progress_state$stack,
                            list(list(id = id, current = NULL, total = NULL, status = NULL)))

  # Log as either a machine readable event, or draw it: the bottom step starts the bar, and a
  # step inside it shows in the bar's status
  if (machine_logging()) {
    emit_event("step_started", step = id, label = step_label(id))
  } else {
    if (length(.progress_state$stack) == 1L) bar_start(id)
    # cli only draws on an update, so a step that never ticks would otherwise never be seen
    bar_render()
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

  # Remove the step from the stack. Steps unwind innermost first, so this is the top entry,
  # but it is found by id so a step left open by mistake cannot close the wrong one
  stack = .progress_state$stack
  hit = which(vapply(stack, function(frame) frame$id, character(1)) == id)
  if (length(hit) > 0L) stack = stack[-hit[[length(hit)]]]
  .progress_state$stack = stack

  if (machine_logging()) {
    emit_event("step_finished", step = id, status = status, message = message)
  } else if (length(stack) == 0L) {
    # The step that owns the bar has finished
    bar_done(failed = failed)
  } else if (!failed) {
    # A step inside the bar's own step has finished. The status the step around it had when
    # this one began is not put back, as that work is over by the time this one finishes.
    # A failed step's status is left on the bar, which fails with it, to show where it went wrong
    set_top(status = NULL)
    bar_render()
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
  if (is.null(current_step())) return(invisible(FALSE))

  # Get the total and current step counters
  total = as_count(config$max_steps)
  current = as_count(state$step)

  # Check the total and valid counts are valid
  if (!is.na(total) && total <= 0L) total = NA_integer_
  if (!is.na(total) && !is.na(current)) current = min(current, total)

  # A tick's message, such as the fold or rep being run, goes beside the bar or on the event
  message = state$message
  status = if (is.character(message) && length(message) == 1L && nzchar(message)) message else NULL

  # Kept so that a status sent later, during work with no ticks, can repeat where the step got
  # to. A tick without a message keeps the status the step already had
  set_top(current = current, total = total)
  if (!is.null(status)) set_top(status = status)

  # Report a tick either in machine readable format, or update the progress bar
  if (machine_logging()) {
    emit_event("progress", step = current_step(), current = current, total = total,
               message = status)
  } else {
    bar_render()
  }
}

# Say what the current step is doing during work that has no ticks of its own, such as
# writing a file or building a plot after the counted work is done
# Shown beside the bar, or sent as a progress event carrying the message
report_status = function(message) {
  frame = top_frame()
  if (is.null(frame)) return(invisible(FALSE))

  if (machine_logging()) {
    # The step's last counts, so a reader does not take the status as the count starting over
    emit_event("progress", step = frame$id, current = frame$current, total = frame$total,
               message = message)
  } else if (progress_drawn()) {
    set_top(status = message)
    bar_render()
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

# --- what each chromosome is doing ----------------------------------------------------------
#
# Blocking works through chromosomes one per worker, so the status beside the bar names the
# chromosomes being worked on and what each is doing, and a line is printed above the bar as
# each one finishes:
#
#   ✔ Chromosome 3: 42 blocks (12 refilled, 0 dropped)
#   Graph Blocking ■■■■■■     28% | Chromosome 7 refilling, Chromosome 12 choosing | 2/16 done
#
# A worker says where a chromosome is up to with signal_chromosome_stage(). That is a
# condition of its own rather than a progress tick: future passes an immediateCondition back
# from a worker as it happens, and catching it around the map sees every one. The progress
# handler would not, since it is throttled and only sees the latest tick.

# Columns the bar takes besides its step name and status: the bar, percentage, ETA and the
# separators either side of the status, with room for a longer ETA
BAR_COLUMNS = 56L

# Say where a chromosome is up to, from wherever it is being worked on, including a worker
# chromosome: the chromosome's label
# stage     : what it is doing, e.g. "refilling", or "done" once it is finished
# ...       : numbers for the line printed when it is done, e.g. n_blocks = 42
signal_chromosome_stage = function(chromosome, stage, ...) {
  cond = structure(
    class = c("hs_chromosome_stage", "immediateCondition", "condition"),
    list(message = "", call = NULL, chromosome = as.character(chromosome), stage = stage,
         details = list(...))
  )
  signalCondition(cond)
  invisible(NULL)
}

# The status for a set of chromosomes, cut back until it fits in width characters
# stages: the stage of each chromosome, named by chromosome, NA for those not started yet
#
#   Chromosome 7 refilling, Chromosome 12 choosing | 5/16 done
#   Chr 7 refilling, Chr 12 choosing | 5/16 done
#   Chr 7 refilling, +1 more | 5/16 done
#   2 active | 5/16 done
#   5/16 done
chromosome_status = function(stages, width = Inf) {
  tally = sprintf("%d/%d done", sum(stages %in% "done"), length(stages))
  active = stages[!is.na(stages) & stages != "done"]
  if (length(active) == 0L) return(tally)

  listed = function(word, n) paste(word, names(active)[seq_len(n)], active[seq_len(n)], collapse = ", ")
  candidates = c(
    paste0(listed("Chromosome", length(active)), " | ", tally),
    paste0(listed("Chr", length(active)), " | ", tally),
    # As many as fit, then a count of the rest
    vapply(rev(seq_len(length(active) - 1L)), function(n) {
      sprintf("%s, +%d more | %s", listed("Chr", n), length(active) - n, tally)
    }, character(1)),
    sprintf("%d active | %s", length(active), tally)
  )
  fits = candidates[nchar(candidates) <= width]
  if (length(fits) > 0L) fits[[1]] else tally
}

# Room for the status beside the bar. Unlimited when no bar is drawn, since machine logging
# sends the status as an event rather than drawing it
status_width = function() {
  if (!progress_drawn()) return(Inf)
  stack = .progress_state$stack
  used = if (length(stack) == 0L) 0L else nchar(step_label(stack[[1]]$id))
  # A step inside the bar's own step puts its name and percentage before its status
  if (length(stack) > 1L) {
    top = stack[[length(stack)]]
    used = used + nchar(nested_status(list(id = top$id, current = 1L, total = 1L))) + nchar(" | ")
  }
  max(0L, cli::console_width() - used - BAR_COLUMNS)
}

# Run expr, which works through chromosomes calling signal_chromosome_stage(), keeping the
# status up to date and printing a line as each chromosome finishes
#
# chromosomes  : every chromosome expr works through, in the order to list them
# describe_done: function(chromosome, details) returning the line to print for a finished
#                chromosome, from the numbers it signalled with "done"
# p            : the progressor expr ticks, which carries the status to the bar
#
# The status goes out on the progressor, as a tick of nothing, rather than through
# report_status(). Stages change far more often than is worth drawing, with many workers, and
# progressr throttles ticks. It keeps the latest status for the next one it shows, so the
# status can lag by the throttle interval but always catches up. A tick is also drawn at once:
# cli draws with messages, which hs_with_progress() holds back when drawn outside a tick
track_chromosome_stages = function(expr, chromosomes, describe_done, p) {
  chromosomes = as.character(chromosomes)
  stages = stats::setNames(rep(NA_character_, length(chromosomes)), chromosomes)
  show = function() p(amount = 0, message = chromosome_status(stages, status_width()))

  # Replaces whatever was shown before, such as the workers starting
  show()
  withCallingHandlers(expr, hs_chromosome_stage = function(cond) {
    stages[[cond$chromosome]] <<- cond$stage
    if (identical(cond$stage, "done")) {
      # A message, so suppressMessages() quietens it. While a bar is up with_progress() holds
      # it back, and the tick show() sends straight after prints it above the bar
      message(describe_done(cond$chromosome, cond$details))
      do.call(report_metric, c(list("chromosome_done", chromosome = cond$chromosome),
                               cond$details))
    }
    show()
  })
}

# progressr's with_progress(), holding back only messages and warnings while the bar is up
# Use it in place of with_progress() wherever a step reports progress.
#
# By default with_progress() holds back every condition, its own progress ticks included.
# Holding anything back makes it clear the bar, print what it held, and draw the bar again, on
# every tick, throttled or not. With nothing printed in between, that is the bar blanking and
# redrawing up to a hundred times a second. Messages and warnings are what can be printed, so
# they are all it needs to hold back
hs_with_progress = function(expr) {
  progressr::with_progress(expr, delay_conditions = c("message", "warning"))
}

# Provides callbacks for progressr including initiation, updates and finish
progress_reporter = function() {
  report = function(config, state, ...) report_tick(config, state)
  # Having hide and unhide also tells progressr it can print held back output while the bar
  # is running, rather than holding all of it until with_progress() exits
  list(initiate = report, update = report, finish = report,
       hide = function(...) bar_hide(), unhide = function(...) bar_unhide())
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
