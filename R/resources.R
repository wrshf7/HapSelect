# R/resources.R
#
# Central resolution of the compute resources HapSelect is allowed to use.
#
# Every parallel entry point in the package asks cpu_cores() for its worker count instead of
# calling detectCores() itself.

##### Set the maximum number of CPU cores HapSelect may use #####
# cores: a positive whole number, or NULL to clear the limit and go back to using every core
#        available to the session.
#
# The value is a ceiling, not a target. HapSelect never uses more cores than are actually
# available, so setting a limit higher than the machine has does nothing.
set_cpu_cores = function(cores = NULL) {
  # Clear the setting if NULL is passed
  if (is.null(cores)) {
    options(HapSelect.cpu_cores = NULL)
    return(invisible(NULL))
  }

  # Validate the value, or stop with an error explaining the limit was not updated.
  resolved = check_cpu_cores(
    cores,
    "The CPU core limit was not changed. Pass NULL to clear it and use every available core."
  )

  # Set the option for the rest of the session
  options(HapSelect.cpu_cores = resolved)
  invisible(resolved)
}

##### Validate a CPU core limit #####
# Accepts a single positive number and returns it as an integer. Anything else is an error,
# since silently ignoring a bad limit would let a job quietly take over a shared machine.
check_cpu_cores = function(cores, consequence) {
  if (length(cores) != 1L || !is.numeric(cores) || is.na(cores) || !is.finite(cores) || cores < 1) {
    stop("A CPU core limit must be a single number greater than or equal to 1.\n", consequence)
  }

  as.integer(floor(cores))
}

##### Number of CPU cores HapSelect will use for parallel work #####
# Resolution order, first hit wins:
#   1. The limit given to set_cpu_cores()
#   2. Every core available to the session
#
# A limit never exceeds the cores actually available, so it cannot be used to oversubscribe a
# machine.
#.
cpu_cores = function() {
  available = max(1L, as.integer(future::availableCores()))

  configured = getOption("HapSelect.cpu_cores")
  if (!is.null(configured)) {
    limit = check_cpu_cores(
      configured,
      "This is the CPU core limit configured for HapSelect. Set a valid one with set_cpu_cores(), or clear it with set_cpu_cores(NULL)."
    )
    return(max(1L, min(available, limit)))
  }

  available
}
