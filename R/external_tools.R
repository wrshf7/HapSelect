# R/external_tools.R
#
# Running the external tools HapSelect drives (PLINK, Beagle) as child processes.
#
# system2() hands a tool's output back only once the tool has exited, so nothing can report
# progress while it runs. Here the output is read as it arrives, one line at a time, so a
# caller can turn the tool's own progress lines into progress ticks.

##### Run an external command, streaming its output line by line #####
# command : path to the executable
# args    : character vector of arguments, passed as is - there is no shell, so nothing needs
#           quoting
# on_line : optional function called with each complete line of output, as it arrives
#
# stderr is merged into stdout, so lines arrive in the order the tool wrote them.
# Returns a list of:
#   status - the exit status
#   output - every line of output, in order
# nocov start
run_tool_command = function(command, args, on_line = NULL) {
  proc = processx::process$new(command, args, stdout = "|", stderr = "2>&1")
  # An interrupt or error in on_line would otherwise leave the tool running on its own
  on.exit(if (proc$is_alive()) proc$kill(), add = TRUE)

  output  = character()
  partial = ""

  # Split what has arrived into lines, holding back a trailing line that is not finished yet
  take = function(chunk, final = FALSE) {
    text = paste0(partial, chunk)
    if (!nzchar(text)) return(invisible(NULL))

    lines = strsplit(text, "\r?\n")[[1]]
    if (!final && !endsWith(text, "\n")) {
      partial <<- lines[length(lines)]
      lines = lines[-length(lines)]
    } else {
      partial <<- ""
    }

    output <<- c(output, lines)
    if (!is.null(on_line)) for (line in lines) on_line(line)
  }

  # Keeps reading after the tool exits, until everything it wrote has been read
  while (proc$is_incomplete_output()) {
    proc$poll_io(200)
    take(proc$read_output())
  }
  take("", final = TRUE)

  proc$wait()
  list(status = proc$get_exit_status(), output = output)
}
# nocov end
