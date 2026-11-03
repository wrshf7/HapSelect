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
# on_line : optional function called with each piece of output, as it arrives. A piece is a
#           line, or a part of one that the tool redraws in place with a carriage return or
#           backspaces, as PLINK does with a running percentage: "... 0%", "1%", "2%"
#
# stderr is merged into stdout, so lines arrive in the order the tool wrote them.
# Returns a list of:
#   status - the exit status
#   output - every line of output, in order, as a terminal would have shown it
# nocov start
run_tool_command = function(command, args, on_line = NULL) {
  proc = processx::process$new(command, args, stdout = "|", stderr = "2>&1")
  # An interrupt or error in on_line would otherwise leave the tool running on its own
  on.exit(if (proc$is_alive()) proc$kill(), add = TRUE)

  output = character()
  line_rest = ""
  piece_rest = ""

  take = function(chunk, final = FALSE) {
    lines = complete_pieces(paste0(line_rest, chunk), "\r?\n", final)
    line_rest <<- lines$rest
    output <<- c(output, render_redraws(lines$done))

    pieces = complete_pieces(paste0(piece_rest, chunk), "[\r\n\b]+", final)
    piece_rest <<- pieces$rest
    if (!is.null(on_line)) for (piece in pieces$done[nzchar(pieces$done)]) on_line(piece)
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

##### Split text into the pieces that have finished arriving #####
# text  : everything read and not yet handed on
# sep   : regular expression a piece ends at
# final : true once the tool has exited, so a last piece without an ending is finished too
# Returns a list of the finished pieces (done) and the unfinished text left over (rest)
complete_pieces = function(text, sep, final = FALSE) {
  if (!nzchar(text)) return(list(done = character(), rest = ""))

  pieces = strsplit(text, sep, perl = TRUE)[[1]]
  ended = grepl(paste0("(", sep, ")$"), text, perl = TRUE)
  if (final || ended) return(list(done = pieces, rest = ""))
  list(done = pieces[-length(pieces)], rest = pieces[length(pieces)])
}

##### Show lines as a terminal would, after the tool redrew parts of them #####
# A backspace moves back one character and a carriage return to the start of the line, and
# what follows writes over what was there. PLINK's "... 0%\b\b1%\b\b2%" becomes "... 2%"
render_redraws = function(lines) {
  vapply(lines, function(line) {
    if (!grepl("[\r\b]", line)) return(line)
    shown = character()
    at = 0L
    for (ch in strsplit(line, "")[[1]]) {
      if (ch == "\b") {
        at = max(0L, at - 1L)
      } else if (ch == "\r") {
        at = 0L
      } else {
        at = at + 1L
        shown[at] = ch
      }
    }
    sub("\\s+$", "", paste(shown, collapse = ""))
  }, character(1), USE.NAMES = FALSE)
}
