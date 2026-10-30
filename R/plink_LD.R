#######################################
### PLINK Installation Verification ###
#######################################

# nocov start

##### Name of the PLINK executable on this platform #####
plink_exe_name = function() {
  # On Windows, the executable is plink.exe; on other platforms it is just plink.
  if (.Platform$OS.type == "windows") {
    return("plink.exe")
  } else {
    return("plink")
  }
}

##### Directory the install scripts place PLINK in #####
# install_linux.sh / install_mac.sh uses $HOME/bin
# install_windows.ps1 uses %USERPROFILE%\bin.
plink_default_install_dir = function() {
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
normalize_plink_path = function(path) {
  # winslash is ignored on non-Windows platforms
  normalizePath(path, winslash = "\\", mustWork = FALSE)
}

##### Normalise a user-supplied PLINK location to an executable path #####
# Accepts either the executable itself or the directory containing it.
resolve_plink_candidate = function(path) {
  path = path.expand(path)

  # Add the executable name if the user supplied a directory. 
  if (dir.exists(path)) {
    path = file.path(path, plink_exe_name())
  }

  return(normalize_plink_path(path))
}

##### Identify an executable by running "plink --version" #####
# PuTTY's SSH link tool is also called plink, so an executable of the right name is not
# necessarily PLINK. Checking the version banner distinguishes them and rules out
# unsupported PLINK versions at the same time:
#   PLINK 1.9 -> "PLINK v1.90b7.2..."
#   PLINK 2.0 -> "PLINK v2.00a5.10LM..."
#   PuTTY     -> "plink: Release 0.83..."
# Returns a list with the outcome ("ok", "wrong_version" or "not_plink") and the banner
identify_plink = function(path) {
  # Get the version output, suppressing warnings and errors
  version_output = tryCatch(
    suppressWarnings(
      system2(path, "--version", stdout = TRUE, stderr = TRUE, timeout = 30)
    ),
    error = function(e) character()
  )

  # Describe whatever was run using its version line, preferring a line that names PLINK
  # because some versions print a decorative header first
  banner = version_output[nzchar(trimws(version_output))]
  named = banner[grepl("PLINK", banner, fixed = TRUE)]
  # Use the first non-empty line of output as the banner, preferring a line that names PLINK
  banner = if (length(named) > 0) trimws(named[1]) else if (length(banner) > 0) trimws(banner[1]) else "no output"

  # HapSelect requires PLINK 1.9, whose banners all begin "PLINK v1.9"
  if (any(grepl("PLINK v1.9", version_output, fixed = TRUE))) {
    return(list(status = "ok", banner = banner))
  }

  # Any other PLINK, such as 2.0 or the much older 1.07, is the wrong version
  if (any(grepl("PLINK", version_output, fixed = TRUE))) {
    return(list(status = "wrong_version", banner = banner))
  }

  # Otherwise, the executable is not PLINK at all (e.g. PuTTY's plink.exe)s
  list(status = "not_plink", banner = banner)
}

##### Describe a rejected candidate for use in an error or warning #####
describe_plink_rejection = function(path, identified) {
  reason = if (identified$status == "wrong_version") {
    "PLINK version is wrong, HapSelect requires PLINK 1.9"
  } else {
    "This may not be PLINK, do you have PuTTY installed?"
  }

  paste0(path, "\n    ", identified$banner, "\n    ", reason)
}

##### Resolve an explicitly given PLINK location, or stop trying #####
# Used for any path the user supplied, set through set_plink_path().
# path: the plink executable, or the directory containing it.
# consequence: completes the error message, describing what will not happen as a result.
verify_plink_path = function(path, consequence) {
  # Check the path is a single non-empty string
  if (!is.character(path) || length(path) != 1 || !nzchar(path)) {
    stop("The PLINK path must be a single non-empty character string.\n", consequence)
  }

  # Attempt to resolve the path to an executable
  resolved = resolve_plink_candidate(path)

  # Check the resolved path exists and is a file
  if (!file.exists(resolved)) {
    stop("No PLINK executable found at: ", resolved, "\n", consequence)
  }

  # Attempt to identify the executable
  identified = identify_plink(resolved)

  # If the identification is not ok, stop with a message describing what was found instead
  if (identified$status != "ok") {
    stop("The given PLINK executable is not PLINK 1.9:\n  ",
         describe_plink_rejection(resolved, identified), "\n", consequence)
  }

  return(resolved)
}

##### Set the PLINK executable used by HapSelect for the rest of the session #####
# path: the plink executable, or the directory containing it. Pass NULL to clear
# the setting and fall back to automatic detection.
set_plink_path = function(path = NULL) {
  # Clear the setting if NULL is passed
  if (is.null(path)) {
    options(HapSelect.plink_path = NULL)
    return(invisible(NULL))
  }

  # Resolve the path and verify it is a PLINK 1.9 executable, or stop with an error explaining the path was not updated.
  resolved = verify_plink_path(
    path,
    "The PLINK path was not changed. Pass NULL to clear it and detect PLINK automatically."
  )

  # Set the option for the rest of the session
  options(HapSelect.plink_path = resolved)
  invisible(resolved)
}

##### Resolve the PLINK executable path #####
# Resolution order, first hit wins:
#   1. The path given to set_plink_path()
#   2. The location used by the install scripts (~/bin)
#   3. Automatic detection across the PATH, confirmed by running "plink --version"
find_plink = function() {
  configured = getOption("HapSelect.plink_path")

  # Verified here as well for good measure, in case the installation was changed or upgraded after the path was set.
  if (!is.null(configured)) {
    return(verify_plink_path(
      configured,
      paste0("This is the PLINK path configured for HapSelect. Point HapSelect at a PLINK 1.9\n",
             "installation with set_plink_path(), or clear the setting with set_plink_path(NULL)\n",
             "to search for one automatically.")
    ))
  }

  # Search the install location from the install scripts, then the PATH, for any executable named plink or plink.exe.
  install_dir = plink_default_install_dir()
  path_dirs = strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  search_dirs = c(install_dir, path_dirs)
  candidates = file.path(search_dirs[nzchar(search_dirs)], plink_exe_name())
  candidates = candidates[file.exists(candidates)]

  # Normalize the candidate paths and remove duplicates, ignoring case on Windows.
  candidates = normalize_plink_path(candidates)
  keys = if (.Platform$OS.type == "windows") tolower(candidates) else candidates
  candidates = candidates[!duplicated(keys)]

  # Take the first candidate that identifies itself as PLINK 1.9, recording why the others
  # were rejected so the error can explain what was found instead.
  rejected = character()
  for (candidate in candidates) {
    identified = identify_plink(candidate)

    # If a valid PLINK 1.9 executable was found, return it immediately.
    if (identified$status == "ok") {
      return(candidate)
    }

    # Add the candidate to the list of rejected executables, along with a description of why it was rejected.
    rejected = c(rejected, describe_plink_rejection(candidate, identified))
  }

  # If no valid PLINK 1.9 executable was found, stop with an error message describing what was found instead.
  stop(
    if (length(rejected) > 0) {
      paste0(
        "No PLINK 1.9 installation was found. Executables named ", plink_exe_name(),
        " were found,\nbut none of them is PLINK 1.9:\n  ",
        paste(rejected, collapse = "\n  "), "\n"
      )
    } else {
      paste0(
        "PLINK executable not found.\n",
        "Expected it in the install location (", file.path(install_dir, plink_exe_name()), ") or on the PATH.\n"
      )
    },
    "Either run the installer for your platform in inst/scripts/install/, or point HapSelect at an\n",
    "existing installation with set_plink_path(\"/path/to/", plink_exe_name(), "\")."
  )
}
# nocov end

#######################################
######## PLINK-Based Functions ########
#######################################

##### Add PLINK's thread limit to its arguments #####
# PLINK uses every core unless told otherwise, so it is held to cpu_cores() like the rest of
# HapSelect. A --threads the caller passed in extra_args is left as it is.
plink_thread_args = function(args) {
  if ("--threads" %in% args) return(args)
  c(args, "--threads", as.character(cpu_cores()))
}

##### Run PLINK with platform-aware executable resolution #####
# nocov start
# on_line: optional function called with each piece of PLINK's output as it arrives
call_plink = function(args, on_line = NULL) {
  run_tool_command(find_plink(), plink_thread_args(args), on_line = on_line)
}
# nocov end

##### Read a PLINK .bim file and assign per-chromosome locus indices #####
read_plink_bim = function(path){
  # Check the bim file exists and is not empty
  if(!file.exists(path)){
    stop("PLINK .bim file not found: ", path, "\nExpected format: tab-delimited, no header, columns: chromosome, SNP ID, cM position, base-pair position (and optionally allele columns).")
  }
  if(file.info(path)$size == 0){
    stop("PLINK .bim file is empty: ", path)
  }

  # Read file. Marker names are read as text, as format_plink_ld() reads them from the .ld
  # file, so a name that looks like a number (e.g. 007) still matches between the two
  bim = utils::read.table(path, header = FALSE, stringsAsFactors = FALSE,
                          colClasses = c(V2 = "character"))

  # Check bim has the required columns
  if(ncol(bim) < 4){
    stop("PLINK .bim file must contain at least 4 columns: chromosome, SNP, cM, position.")
  }

  # Only the first 4 columns are needed, drop the rest
  bim = bim[, 1:4]

  # Rename columns
  colnames(bim) = c("Chrom", "SNP", "CM", "Position")

  # Convert row values to numeric
  bim$Chrom = type.convert(bim$Chrom, as.is = TRUE)
  bim$Position = as.numeric(bim$Position)

  # Assign marker indices that restart at 1 within each chromosome.
  bim$Locus = ave(seq_len(nrow(bim)), bim$Chrom, FUN = seq_along)

  return(bim)
}

##### Convert a PLINK .ld output file into the HapSelect LD format #####
# bim can be either a path to a .bim file or the parsed data frame from read_plink_bim().
# If you wish to run PLINK manually, this function is needed to convert the PLINK LD output into fast stack ld format #
format_plink_ld = function(ld_path, bim){
  # Check ld output file exists
  if(!file.exists(ld_path)){
    stop("PLINK did not produce an .ld output file.")
  }

  # If bim is a path, read it
  if(is.character(bim) && length(bim) == 1){
    bim = read_plink_bim(bim)
  }

  if(file.info(ld_path)$size == 0) {
    # No LD pairs were reported, so return an empty HapSelect table.
    return(data.frame(
      Chrom = bim$Chrom[FALSE],
      Locus1 = integer(),
      Locus2 = integer(),
      Name1 = character(),
      Name2 = character(),
      LD = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  # Ensure columns are valid, from the header alone
  required_cols = c("SNP_A", "SNP_B", "R2")
  if(!all(required_cols %in% names(data.table::fread(ld_path, nrows = 0)))){
    stop("PLINK .ld output must contain the columns SNP_A, SNP_B, and R2.")
  }

  # Read the ld output file. fread rather than read.table, since a genome-wide table runs to
  # gigabytes, and only the three columns used below are read. Marker names are always read
  # as text, so a name that looks like a number keeps its exact spelling.
  ld = data.table::fread(ld_path, select = required_cols,
                         colClasses = list(character = c("SNP_A", "SNP_B")),
                         data.table = FALSE, showProgress = FALSE)

  # Find the position of each SNP in the bim table
  idx1 = match(ld$SNP_A, bim$SNP)
  idx2 = match(ld$SNP_B, bim$SNP)

  # Check all LD SNPs are found in the bim file
  if(any(is.na(idx1)) || any(is.na(idx2))){
    stop("One or more SNPs in the PLINK .ld output were not found in the .bim file.")
  }

  # Look up the chromosome of each LD SNP pair
  chrom1 = bim$Chrom[idx1]
  chrom2 = bim$Chrom[idx2]

  # Filter cross-chromosome pairs, which PLINK should not produce but occasionally does
  if(any(chrom1 != chrom2)){
    warning("PLINK .ld output contained cross-chromosome marker pairs; they will be dropped.")
    keep   = chrom1 == chrom2
    ld     = ld[keep, ]
    idx1   = idx1[keep]
    idx2   = idx2[keep]
    chrom1 = chrom1[keep]
  }

  # Rebuild the same long-form LD structure used by HapSelect's internal LD path.
  ld_df = data.frame(
    Chrom = chrom1,
    Locus1 = bim$Locus[idx1],
    Locus2 = bim$Locus[idx2],
    Name1 = ld$SNP_A,
    Name2 = ld$SNP_B,
    LD = ld$R2,
    stringsAsFactors = FALSE
  )

  # Sort by chromosome, then loci
  ld_df = ld_df[order(ld_df$Chrom, ld_df$Locus1, ld_df$Locus2), ]
  row.names(ld_df) = NULL

  return(ld_df)
}

##### Runs a plink command given a set of arguments #####
# nocov start
# on_line is passed through to call_plink(), to follow PLINK's progress as it runs
run_plink_command = function(args, on_line = NULL){
  # Run the plink command
  result = call_plink(args, on_line = on_line)

  # If the command outputted a non-zero status code, something went wrong
  if(result$status != 0){
    # Stop execution, plink failed
    stop(
      "PLINK LD calculation failed.\n",
      paste(result$output, collapse = "\n")
    )
  }

  invisible(result$output)
}
# nocov end

##### Progress through the stages of a PLINK run, on one bar #####
# Each stage, such as writing the input or running --r2, gets an equal share of the bar
# whatever it takes, since how long the stages take relative to one another depends on the
# data. A stage moves through its share as it reports, and is filled to the end of its share
# when it finishes. The bar only moves forward.
#
# p: progressor with PLINK_STAGE_STEPS steps per stage
# Returns a list of:
#   within(fraction, message) - the current stage is this far through, from 0 to 1
#   finish()                  - the current stage is done, so the next one starts
PLINK_STAGE_STEPS = 100L

plink_stages = function(p) {
  finished = 0L
  ticked = 0L

  move_to = function(to, message = NULL) {
    to = as.integer(floor(to))
    if (to > ticked) {
      p(amount = to - ticked, message = message)
      ticked <<- to
    } else if (!is.null(message)) {
      p(amount = 0, message = message)
    }
  }

  list(
    within = function(fraction, message = NULL) {
      move_to(PLINK_STAGE_STEPS * (finished + min(max(fraction, 0), 1)), message)
    },
    finish = function() {
      finished <<- finished + 1L
      move_to(PLINK_STAGE_STEPS * finished)
    }
  )
}

##### Follow PLINK's running percentages through one command #####
# PLINK works through a command in passes, drawing a percentage for each in place:
#   Scanning .ped file... 0%  1%  ...  100%.ped scan complete (for binary autoconversion).
#   Performing single-pass .bed write (6070 variants, 444 people).
#   0%  1%  ...  94%--file: x-temporary.bed + x-temporary.bim + x-temporary.fam written.
#   --r2 to x.ld ... 0% [processing]  writing]  done.
# run_tool_command() hands these over a percentage at a time. A pass starts when a
# percentage follows a "..." or falls back below the last one, and the stage is split
# equally between the passes the command is expected to make. A command that makes more
# than expected fills its stage early and stays there, rather than overrunning it.
#
# stages : from plink_stages()
# passes : the passes the command is expected to make
# Returns a function taking one piece of output, for run_plink_command()'s on_line
#
# The passes HapSelect's commands make, each matched by how PLINK starts its line, and what
# to call it
PLINK_PASS_NAMES = list(
  c(plink = "^Scanning \\.ped file",               name = "Reading genotypes"),
  c(plink = "^Performing .*\\.bed write",           name = "Converting genotypes"),
  c(plink = "^Calculating allele frequencies",      name = "Calculating allele frequencies"),
  c(plink = "^--make-bed",                          name = "Writing PLINK files"),
  c(plink = "^--r2",                                name = "Calculating LD")
)

plink_progress_parser = function(stages, passes) {
  pass = -1L
  last = -1
  label = NULL
  text = NULL

  # What a pass does, said for a person rather than in PLINK's words. A pass PLINK names in
  # some other way, as a new version might, is just PLINK running
  pass_name = function(x) {
    x = trimws(x)
    for (pass in PLINK_PASS_NAMES) {
      if (grepl(pass[["plink"]], x)) return(pass[["name"]])
    }
    "Running PLINK"
  }

  function(piece) {
    labelled = regmatches(piece, regexec("^(.*\\S)\\s*\\.\\.\\.\\s*([0-9]{1,3})%", piece))[[1]]
    bare = regmatches(piece, regexec("^\\s*([0-9]{1,3})%", piece))[[1]]

    if (length(labelled) == 3) {
      pct = as.numeric(labelled[[3]])
      pass <<- pass + 1L
      label <<- pass_name(labelled[[2]])
    } else if (length(bare) == 2) {
      pct = as.numeric(bare[[2]])
      if (pct < last || pass < 0L) {
        pass <<- pass + 1L
        label <<- if (is.null(text)) NULL else pass_name(text)
      }
    } else {
      # Text between passes names the pass that follows it when that pass has no name of its own
      if (nzchar(trimws(piece))) text <<- piece
      return(invisible(NULL))
    }

    last <<- pct
    message = paste0(if (is.null(label)) "Running PLINK" else label, " ", pct, "%")
    stages$within((pass + pct / 100) / passes, message)
    invisible(NULL)
  }
}


##### Write PLINK text input files (.ped / .map) from a genotype data frame #####
# geno: data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#       cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
# progress: optional function called with the fraction of individuals written so far
write_plink_ped_map = function(geno, prefix, progress = NULL) {

  # PED stores two discrete alleles per genotype, so a dosage has to be a whole number of
  # ALT copies. A fractional dosage has no representation here, and the conversion below
  # would read every one of them as a heterozygote whatever its value, so a dosage just
  # above zero and one just below two would come out identical.
  dosages = as.matrix(geno[, -(1:3), drop = FALSE])
  observed = dosages[!is.na(dosages)]
  if (length(observed) > 0 && any(observed != round(observed))) {
    stop("geno holds fractional dosages, which PLINK's text format cannot represent: each ",
         "one would be written as a heterozygote. read_vcf_geno(prefer_ds = TRUE) returns ",
         "fractional dosages for imputed data - re-read the VCF with prefer_ds = FALSE to ",
         "get whole-number calls before passing it to PLINK.")
  }

  utils::write.table(
    # Positions in full, as write_vcf_geno() writes them, rather than e.g. 5e+05
    data.frame(CHR = geno[[2]], SNP = geno[[1]], CM = 0,
               BP = format(geno[[3]], scientific = FALSE, trim = TRUE),
               stringsAsFactors = FALSE),
    file      = paste0(prefix, ".map"),
    quote     = FALSE,
    sep       = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

  dosage_to_calls = function(x) {
    a1 = ifelse(is.na(x), "0", ifelse(x == 2L, "G", "A"))
    a2 = ifelse(is.na(x), "0", ifelse(x == 0L, "A", "G"))
    c(rbind(a1, a2))
  }

  geno_matrix = as.matrix(geno[, -(1:3)])

  # Written a few individuals at a time, which is where the time goes on a large set, so
  # progress can be reported between them. The file is the same as writing it in one go
  ped_file = file(paste0(prefix, ".ped"), open = "w")
  on.exit(close(ped_file), add = TRUE)
  individuals = seq_len(ncol(geno_matrix))
  chunks = split(individuals, ceiling(individuals / max(1L, ceiling(length(individuals) / 50))))

  for (chunk in chunks) {
    ped = do.call(rbind, lapply(chunk, function(i) {
      c(paste0("F", i), paste0("I", i), "0", "0", "0", "-9",
        dosage_to_calls(geno_matrix[, i]))
    }))

    utils::write.table(
      ped,
      file      = ped_file,
      quote     = FALSE,
      sep       = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
    if (!is.null(progress)) progress(chunk[length(chunk)] / length(individuals))
  }
}

##### Run PLINK pairwise LD on a binary fileset and return it in HapSelect format #####
# Shared by plink_pairwise_ld() and plink_pairwise_ld_geno(), which each report the
# plink_ld step themselves, so a call reports one step whichever input it starts from.
# stages: the caller's plink_stages(), moved through the two stages run here: running --r2,
#         then reading its output
plink_pairwise_ld_run = function(prefix, ld_window, ld_window_kb, ld_window_r2, extra_args,
                                 stages){
  required_files = paste0(prefix, c(".bed", ".bim", ".fam"))
  missing_files  = required_files[!file.exists(required_files)]

  if(length(missing_files) > 0){
    stop(
      "Missing required PLINK input files: ",
      paste(basename(missing_files), collapse = ", ")
    )
  }

  out_prefix = tempfile("hapselect_plink_out_")
  on.exit(
    unlink(paste0(out_prefix, c(".ld", ".log", ".nosex")), force = TRUE),
    add = TRUE
  )

  args = c(
    "--bfile", prefix,
    "--r2",
    "--ld-window",    as.character(ld_window),
    "--ld-window-kb", as.character(ld_window_kb),
    "--ld-window-r2", as.character(ld_window_r2),
    extra_args,
    "--out", out_prefix
  )

  # Allele frequencies, then the --r2 pass itself
  stages$within(0, "Calculating LD")
  run_plink_command(args, on_line = plink_progress_parser(stages, passes = 2L))
  stages$finish()

  # Nothing to count while the table is read, which on a large one takes longer than PLINK
  stages$within(0, "Reading LD results")
  ld = format_plink_ld(paste0(out_prefix, ".ld"), paste0(prefix, ".bim"))
  stages$finish()

  ld
}

##### Run PLINK pairwise LD and return the result in HapSelect format #####
# prefix should point to a PLINK binary fileset without extension (.bed/.bim/.fam)
plink_pairwise_ld = function(prefix, ld_window = 999999, ld_window_kb = 1000000,
                             ld_window_r2 = 0, extra_args = character()){
  report_step("plink_ld")

  # Two stages on one bar: running --r2, then reading its output
  with_progress({
    stages = plink_stages(progressor(steps = 2L * PLINK_STAGE_STEPS))
    ld = plink_pairwise_ld_run(prefix, ld_window = ld_window, ld_window_kb = ld_window_kb,
                               ld_window_r2 = ld_window_r2, extra_args = extra_args,
                               stages = stages)
  })

  ld
}

##### Run PLINK pairwise LD from a genotype data frame #####
# Writes temporary PLINK text files, delegates to plink_pairwise_ld_run, then cleans up.
# geno: data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#       cols 4+ = dosage values (0 / 1 / 2 / NA) per individual. PLINK reads hard calls, so
#       these must be whole numbers - read a VCF with prefer_ds = FALSE rather than TRUE,
#       whose fractional imputed dosages PLINK's text format cannot carry.
plink_pairwise_ld_geno = function(geno, ld_window = 999999, ld_window_kb = 1000000,
                                  ld_window_r2 = 0, extra_args = character()){
  # The same step as plink_pairwise_ld(), opened here so it also covers writing and
  # converting the input
  report_step("plink_ld")

  if(!is.data.frame(geno) || ncol(geno) < 4){
    stop("geno must be a data frame with columns: marker, chromosome, position, and at least one genotype column.")
  }

  #extract autosome number
  chr_num = max(geno[,2])

  in_prefix = tempfile("hapselect_plink_in_")
  on.exit(unlink(paste0(in_prefix, c(".ped", ".map", ".bed", ".bim", ".fam", ".log", ".nosex")),
                 force = TRUE), add = TRUE)

  # Four stages on one bar: writing the input and converting it, ahead of the two stages
  # shared with plink_pairwise_ld()
  with_progress({
    stages = plink_stages(progressor(steps = 4L * PLINK_STAGE_STEPS))

    stages$within(0, "Writing genotypes for PLINK")
    write_plink_ped_map(geno, in_prefix,
                        progress = function(fraction) stages$within(fraction, "Writing genotypes for PLINK"))
    stages$finish()

    # Convert text files to binary so plink_pairwise_ld_run can consume them. PLINK scans the
    # .ped, writes the .bed, counts allele frequencies, then writes the final fileset
    stages$within(0, "Converting to PLINK format")
    run_plink_command(c("--file", in_prefix, "--make-bed", "--out", in_prefix, "--chr-set", as.character(chr_num)),
                      on_line = plink_progress_parser(stages, passes = 4L))
    stages$finish()

    extra_args = c("--chr-set", as.character(chr_num), extra_args)

    ld = plink_pairwise_ld_run(in_prefix, ld_window = ld_window, ld_window_kb = ld_window_kb,
                               ld_window_r2 = ld_window_r2, extra_args = extra_args,
                               stages = stages)
  })

  ld
}
