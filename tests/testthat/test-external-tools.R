# Tests: running external tools -----------------------------------------------
#
# run_tool_command() is what PLINK and Beagle run through. These drive it with Rscript as
# the child process, which is available wherever the tests run, on any platform.

rscript <- function() file.path(R.home("bin"), "Rscript")

# Run R code in a child Rscript through run_tool_command()
run_r <- function(code, on_line = NULL) {
  run_tool_command(rscript(), c("--vanilla", "-e", code), on_line = on_line)
}

test_that("run_tool_command returns every line of output and the exit status", {
  result <- run_r('cat("one\\ntwo\\n"); quit(status = 3)')

  expect_equal(result$output, c("one", "two"))
  expect_equal(result$status, 3)
})


test_that("run_tool_command hands each line over as it arrives, in order", {
  seen <- character()
  # The pauses split the output across reads, including part way through a line
  code <- 'cat("one\\n"); flush(stdout()); Sys.sleep(0.3); cat("tw"); flush(stdout());
           Sys.sleep(0.3); cat("o\\nthree")'

  result <- run_r(code, on_line = function(line) seen <<- c(seen, line))

  # A line split across reads arrives whole, and a last line without a newline still arrives
  expect_equal(seen, c("one", "two", "three"))
  expect_equal(result$output, seen)
  expect_equal(result$status, 0)
})


test_that("run_tool_command merges the tool's errors into its output", {
  # message() writes to stderr, which a failed tool's error message has to include
  result <- run_r('cat("out\\n"); flush(stdout()); Sys.sleep(0.2); message("err")')

  expect_equal(result$output, c("out", "err"))
})


test_that("run_tool_command stops the tool when following its output fails", {
  marker <- tempfile("still_running_")
  # Writes the file only if it outlives the error by a second
  code <- sprintf('cat("start\\n"); flush(stdout()); Sys.sleep(1.5); file.create("%s")',
                  gsub("\\\\", "/", marker))

  expect_error(run_r(code, on_line = function(line) stop("parser failed")), "parser failed")

  Sys.sleep(2)
  expect_false(file.exists(marker))
})
