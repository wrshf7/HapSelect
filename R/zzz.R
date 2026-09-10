# R/zzz.R
#
# Load and attach hooks. .onLoad registers the progress handler and runs whenever the
# namespace loads; .onAttach prints the startup banner and runs only when the package is
# attached, by library(HapSelect) or requireNamespace(..., attach = TRUE).
#

.onLoad <- function(libname, pkgname) {
  # Only runs if progressr's own environment variable is NOT set, as the consumer has already decided
  if (!nzchar(Sys.getenv("R_PROGRESSR_ENABLE", unset = ""))) {
    # Set the progressr.enable to true if machine logging is enabled, the session is interactive or a terminal that is unpiped
    options(progressr.enable = machine_logging() || interactive() || isatty(stderr()))
  }

  # Register the progress handler, either outputting CLI based progress bars, or machine readable progress JSON if enabled
  register_handler()

  # Return nothing by convention
  invisible()
}
 
.onAttach <- function(libname, pkgname) {
  version  <- as.character(utils::packageVersion("HapSelect"))
  docs_url <- sprintf(
    "https://wrshf7.github.io/HapSelect-Docs/%s/",
    version
  )
 
  version_line <- sprintf(
    "                     %s     ",
    cli::format_inline("HapSelect {.dim v{version}}")
  )

  banner <- c(
    "",
    "                               o",
    "                            o-----o",
    "                           o-------o",
    "                            o-----o",
    "                               o",
    "                            o-----o",
    "                           o-------o",
    "                            o-----o",
    "                               o",
    "",
    "            __  __           _____      __          __",
    "           / / / /___ _____ / ___/___  / /__  _____/ /_",
    "          / /_/ / __ `/ __ \\\\__ \\/ _ \\/ / _ \\/ ___/ __/",
    "         / __  / /_/ / /_/ /__/ /  __/ /  __/ /__/ /_",
    "        /_/ /_/\\__,_/ .___/____/\\___/_/\\___/\\___/\\__/",
    "                   /_/",
    "",
    version_line,
    "           Will Shaffer, Zane Carter & Victor Papin"
  )
  packageStartupMessage(paste(banner, collapse = "\n"))
 
  cli::cat_line()
  old_width <- options(cli.width = 64)
  cli::cli_h1("")
  options(old_width)
  cli::cli_text(
    "{.emph Haplotype-based selection and genomic prediction}"
  )
  cli::cat_line()
  cli::cli_alert_success("Welcome to HapSelect!")
  cli::cli_text(
    "Documentation: {.url {docs_url}}"
  )
  cli::cli_alert_warning(paste(
    "Please use the documentation for your installed version.",
    "Features and functions may differ between versions."
  ))
  cli::cat_line()
  cli::cli_text(
    "{.strong Questions & discussion:} ",
    "{.url https://github.com/wrshf7/HapSelect/discussions}"
  )
  cli::cli_text(
    "{.strong Report bugs:} ",
    "{.url https://github.com/wrshf7/HapSelect/issues}"
  )
  cli::cat_line()
}