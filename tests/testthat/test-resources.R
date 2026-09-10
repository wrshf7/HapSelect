# Tests: cpu_cores() and set_cpu_cores() ---------------------------------------
#
# The contract these pin down:
#   - with nothing configured, HapSelect uses every core available to the session
#   - a limit is honoured, but can never oversubscribe the machine
#   - a bad limit is an error, never silently ignored

available_cores <- function() max(1L, as.integer(future::availableCores()))

# Run body with no core limit configured, then restore whatever was set before.
with_no_limit <- function(body) {
  old_option <- getOption("HapSelect.cpu_cores")
  on.exit(options(HapSelect.cpu_cores = old_option), add = TRUE)

  options(HapSelect.cpu_cores = NULL)
  force(body)
}


# Default ----------------------------------------------------------------------

test_that("with no limit configured every available core is used", {
  with_no_limit({
    expect_equal(cpu_cores(), available_cores())
  })
})


# Setting a limit --------------------------------------------------------------

test_that("set_cpu_cores() limits the cores used and returns the value set", {
  with_no_limit({
    expect_equal(set_cpu_cores(1), 1L)
    expect_equal(cpu_cores(), 1L)
  })
})


test_that("set_cpu_cores(NULL) clears the limit", {
  with_no_limit({
    set_cpu_cores(1)
    expect_null(set_cpu_cores(NULL))
    expect_equal(cpu_cores(), available_cores())
  })
})


test_that("a limit is a ceiling, so it cannot oversubscribe the machine", {
  with_no_limit({
    # Asking for far more cores than exist must not start a worker per requested core.
    set_cpu_cores(available_cores() + 1000)
    expect_equal(cpu_cores(), available_cores())
  })
})


test_that("a fractional limit is floored rather than rejected", {
  with_no_limit({
    expect_equal(set_cpu_cores(2.9), 2L)
  })
})


# Bad input --------------------------------------------------------------------

test_that("an invalid limit is an error rather than being ignored", {
  with_no_limit({
    # Silently ignoring these would let a job quietly take over a shared machine.
    expect_error(set_cpu_cores(0), "greater than or equal to 1")
    expect_error(set_cpu_cores(-4), "greater than or equal to 1")
    expect_error(set_cpu_cores("four"), "single number")
    expect_error(set_cpu_cores(c(2, 4)), "single number")
    expect_error(set_cpu_cores(NA), "single number")
    expect_error(set_cpu_cores(Inf), "single number")

    # A rejected value must leave the previous setting alone.
    expect_equal(cpu_cores(), available_cores())
  })
})
