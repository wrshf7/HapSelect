########################################
### Beagle Installation Verification ###
########################################

# nocov start

##### Name of the java executable on this platform #####
java_exe_name = function() {
  if (.Platform$OS.type == "windows") {
    return("java.exe")
  } else {
    return("java")
  }
}

##### Find a Java runtime to run Beagle with #####
# Beagle is distributed as a .jar file rather than a native executable, so it's always run as
# `java -jar beagle.jar ...`. That means a Java runtime (JRE 8 or later) has to be on the PATH
# before Beagle can be run at all, separately from finding the Beagle .jar itself.
find_java = function() {
  java_path = Sys.which(java_exe_name())

  if (!nzchar(java_path)) {
    stop(
      "No Java runtime found on the PATH.\n",
      "Beagle requires a Java runtime (JRE 8 or later) to run. Install one and ensure `",
      java_exe_name(), "` is on the PATH."
    )
  }

  return(unname(java_path))
}

##### Directory the install scripts place Beagle in #####
# Mirrors plink_default_install_dir(): install_linux.sh / install_windows.ps1 place PLINK in
# $HOME/bin. Beagle has no installer of its own (it's just a .jar you download), so it's expected
# to live in that same directory, alongside PLINK.
beagle_default_install_dir = function() {
  home = if (.Platform$OS.type == "windows") {
    Sys.getenv("USERPROFILE", unset = Sys.getenv("HOME", unset = ""))
  } else {
    Sys.getenv("HOME", unset = "")
  }

  # If the home directory is not set, return an empty string so the caller can skip this candidate.
  if (!nzchar(home)) {
    return("")
  }

  file.path(home, "bin")
}

##### Puts a path into a canonical form for comparison and display #####
# The same file can be reached through different spellings, especially on Windows: mixed
# / and \ separators, differing case, and short 8.3 components such as PROGRA~1. Paths are
# normalised before being compared or shown to the user.
normalize_beagle_path = function(path) {
  # winslash is ignored on non-Windows platforms
  normalizePath(path, winslash = "\\", mustWork = FALSE)
}

##### Normalise a user-supplied Beagle location to a .jar path #####
# Accepts either the .jar itself or the directory containing it. If a directory is given and it
# contains exactly one beagle*.jar, that file is used; otherwise "beagle.jar" is assumed, so the
# caller gets a clear "file not found" error rather than silently picking the wrong jar.
resolve_beagle_candidate = function(path) {
  path = path.expand(path)

  if (dir.exists(path)) {
    jars = list.files(path, pattern = "^beagle.*\\.jar$", full.names = TRUE, ignore.case = TRUE)

    if (length(jars) == 1) {
      path = jars[1]
    } else {
      path = file.path(path, "beagle.jar")
    }
  }

  return(normalize_beagle_path(path))
}

##### Identify a .jar as Beagle 5.5 by running it and reading its banner #####
# Running a Beagle .jar with no arguments prints a version banner and usage message, then exits
# with a non-zero status (it hasn't been told what to do). For example:
#
#   beagle.27Feb25.75f.jar (version 5.5)
#   Copyright (C) 2014-2024 Brian L. Browning
#   Usage: java -jar beagle.27Feb25.75f.jar [arguments]
#
# HapSelect is built against Beagle 5.5, so a jar is only accepted if its banner contains
# "(version 5.5)". Any other beagle*.jar is rejected as the wrong version rather than silently
# used, since Beagle's argument names and output format can change between major versions.
identify_beagle = function(path) {
  version_output = tryCatch(
    suppressWarnings(
      system2(find_java(), args = c("-jar", shQuote(path)), stdout = TRUE, stderr = TRUE, timeout = 30)
    ),
    error = function(e) character()
  )

  # Prefer a line that names Beagle for the banner shown to the user, since Beagle's usage text
  # has multiple lines and the version line is the useful one to show.
  banner = version_output[nzchar(trimws(version_output))]
  named = banner[grepl("beagle", banner, ignore.case = TRUE)]
  banner = if (length(named) > 0) trimws(named[1]) else if (length(banner) > 0) trimws(banner[1]) else "no output"

  # HapSelect requires Beagle 5.5, whose banners all contain "(version 5.5)"
  if (any(grepl("(version 5.5)", version_output, fixed = TRUE))) {
    return(list(status = "ok", banner = banner))
  }

  # A banner that names Beagle but not version 5.5 is a real Beagle jar, just the wrong release
  if (any(grepl("beagle", version_output, ignore.case = TRUE))) {
    return(list(status = "wrong_version", banner = banner))
  }

  # Otherwise, this doesn't look like a Beagle jar at all
  list(status = "not_beagle", banner = banner)
}

##### Describe a rejected candidate for use in an error or warning #####
describe_beagle_rejection = function(path, identified) {
  reason = if (identified$status == "wrong_version") {
    "Beagle version is wrong, HapSelect requires Beagle 5.5"
  } else {
    "This does not look like a Beagle .jar file"
  }

  paste0(path, "\n    ", identified$banner, "\n    ", reason)
}

##### Resolve an explicitly given Beagle location, or stop trying #####
# Used for any path the user supplied, set through set_beagle_path().
# path: the Beagle .jar, or the directory containing it.
# consequence: completes the error message, describing what will not happen as a result.
verify_beagle_path = function(path, consequence) {
  if (!is.character(path) || length(path) != 1 || !nzchar(path)) {
    stop("The Beagle path must be a single non-empty character string.\n", consequence)
  }

  resolved = resolve_beagle_candidate(path)

  if (!file.exists(resolved)) {
    stop("No Beagle .jar found at: ", resolved, "\n", consequence)
  }

  identified = identify_beagle(resolved)

  if (identified$status != "ok") {
    stop("The given file is not Beagle 5.5:\n  ", describe_beagle_rejection(resolved, identified), "\n", consequence)
  }

  return(resolved)
}

##### Point HapSelect at a specific Beagle .jar for the rest of the session #####
# path: the Beagle .jar, or the directory containing it. Pass NULL to clear the setting and go
# back to detecting Beagle automatically.
set_beagle_path = function(path = NULL) {
  if (is.null(path)) {
    options(HapSelect.beagle_path = NULL)
    return(invisible(NULL))
  }

  resolved = verify_beagle_path(
    path,
    "The Beagle path was not changed. Pass NULL to clear it and detect Beagle automatically."
  )

  options(HapSelect.beagle_path = resolved)
  invisible(resolved)
}

##### Work out which Beagle .jar to use #####
# Checked in order, first hit wins:
#   1. Whatever was set with set_beagle_path()
#   2. The install scripts' shared tools directory (~/bin, alongside PLINK)
#   3. Every directory on the PATH, checked by actually running each beagle*.jar found there
find_beagle = function() {
  configured = getOption("HapSelect.beagle_path")

  # Re-verified here too, in case the installation changed or was upgraded after the path was set.
  if (!is.null(configured)) {
    return(verify_beagle_path(
      configured,
      paste0("This is the Beagle path configured for HapSelect. Point HapSelect at a Beagle 5.5 .jar\n",
             "with set_beagle_path(), or clear the setting with set_beagle_path(NULL)\n",
             "to search for one automatically.")
    ))
  }

  # Search the install location from the install scripts, then the PATH, for any beagle*.jar file.
  install_dir = beagle_default_install_dir()
  path_dirs = strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  search_dirs = c(install_dir, path_dirs)
  search_dirs = search_dirs[nzchar(search_dirs) & dir.exists(search_dirs)]

  candidates = unlist(lapply(search_dirs, function(dir) {
    list.files(dir, pattern = "^beagle.*\\.jar$", full.names = TRUE, ignore.case = TRUE)
  }))

  # Normalize the candidate paths and remove duplicates, ignoring case on Windows.
  candidates = normalize_beagle_path(candidates)
  keys = if (.Platform$OS.type == "windows") tolower(candidates) else candidates
  candidates = candidates[!duplicated(keys)]

  # Take the first candidate that identifies itself as Beagle 5.5, recording why the others were
  # rejected so the error can explain what was found instead.
  rejected = character()
  for (candidate in candidates) {
    identified = identify_beagle(candidate)

    if (identified$status == "ok") {
      return(candidate)
    }

    rejected = c(rejected, describe_beagle_rejection(candidate, identified))
  }

  # Nothing usable was found - explain what (if anything) was found instead, and how to fix it.
  stop(
    if (length(rejected) > 0) {
      paste0(
        "No usable Beagle installation was found. Files matching beagle*.jar were found,\n",
        "but none of them is Beagle 5.5:\n  ",
        paste(rejected, collapse = "\n  "), "\n"
      )
    } else {
      paste0(
        "Beagle .jar not found.\n",
        "Expected it in the install location (", file.path(install_dir, "beagle.jar"), ") or on the PATH.\n"
      )
    },
    "Download Beagle 5.5 from https://faculty.washington.edu/browning/beagle/beagle.html and place\n",
    "it in ", install_dir, ", or point HapSelect at it with set_beagle_path(\"/path/to/beagle.jar\")."
  )
}

########################################
######## Beagle-Based Functions ########
########################################

##### Run Beagle, automatically finding Java and the Beagle jar #####
# args: character vector of Beagle's key=value arguments, e.g. c("gt=in.vcf.gz", "out=imputed")
call_beagle = function(args, stdout = TRUE, stderr = TRUE) {
  system2(find_java(), args = c("-jar", shQuote(find_beagle()), args), stdout = stdout, stderr = stderr)
}

##### Run a Beagle command and stop with Beagle's own output if it fails #####
run_beagle_command = function(args){
  # Run the beagle command
  beagle_output = call_beagle(args, stdout = TRUE, stderr = TRUE)
  # Get the output code
  beagle_status = attr(beagle_output, "status")

  # If the command failed or outputted a non-zero status code, something went wrong
  if(!is.null(beagle_status) && beagle_status != 0){
    # Stop execution, beagle failed
    stop(
      "Beagle imputation/phasing failed.\n",
      paste(beagle_output, collapse = "\n")
    )
  }

  invisible(beagle_output)
}
# nocov end
