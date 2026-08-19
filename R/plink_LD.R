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

##### Run PLINK with platform-aware executable resolution #####
# nocov start
call_plink = function(args, stdout = TRUE, stderr = TRUE) {
  system2(find_plink(), args = args, stdout = stdout, stderr = stderr)
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

  # Read file
  bim = utils::read.table(path, header = FALSE, stringsAsFactors = FALSE)

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

  # Read the ld output file
  ld = utils::read.table(ld_path, header = TRUE, stringsAsFactors = FALSE)

  # Ensure columns are valid
  required_cols = c("SNP_A", "SNP_B", "R2")
  if(!all(required_cols %in% names(ld))){
    stop("PLINK .ld output must contain the columns SNP_A, SNP_B, and R2.")
  }

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
run_plink_command = function(args){
  # Run the plink command
  plink_output = call_plink(args, stdout = TRUE, stderr = TRUE)
  # Get the output code
  plink_status = attr(plink_output, "status")

  # If the command failed or outputted a non-zero status code, something went wrong
  if(!is.null(plink_status) && plink_status != 0){
    # Stop execution, plink failed
    stop(
      "PLINK LD calculation failed.\n",
      paste(plink_output, collapse = "\n")
    )
  }

  invisible(plink_output)
}
# nocov end


##### Write PLINK text input files (.ped / .map) from a genotype data frame #####
# geno: data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#       cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
write_plink_ped_map = function(geno, prefix) {
  utils::write.table(
    data.frame(CHR = geno[[2]], SNP = geno[[1]], CM = 0, BP = geno[[3]],
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
  ped = do.call(rbind, lapply(seq_len(ncol(geno_matrix)), function(i) {
    c(paste0("F", i), paste0("I", i), "0", "0", "0", "-9",
      dosage_to_calls(geno_matrix[, i]))
  }))

  utils::write.table(
    ped,
    file      = paste0(prefix, ".ped"),
    quote     = FALSE,
    sep       = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}

##### Run PLINK pairwise LD and return the result in HapSelect format #####
# prefix should point to a PLINK binary fileset without extension (.bed/.bim/.fam)
plink_pairwise_ld = function(prefix, ld_window = 999999, ld_window_kb = 1000000,
                             ld_window_r2 = 0, extra_args = character()){
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

  run_plink_command(args)

  format_plink_ld(paste0(out_prefix, ".ld"), paste0(prefix, ".bim"))
}

##### Run PLINK pairwise LD from a genotype data frame #####
# Writes temporary PLINK text files, delegates to plink_pairwise_ld, then cleans up.
# geno: data frame with col 1 = marker name, col 2 = chromosome, col 3 = position,
#       cols 4+ = dosage values (0 / 1 / 2 / NA) per individual
plink_pairwise_ld_geno = function(geno, ld_window = 999999, ld_window_kb = 1000000,
                                  ld_window_r2 = 0, extra_args = character()){
  if(!is.data.frame(geno) || ncol(geno) < 4){
    stop("geno must be a data frame with columns: marker, chromosome, position, and at least one genotype column.")
  }

  #extract autosome number
  chr_num = max(geno[,2])

  in_prefix = tempfile("hapselect_plink_in_")
  on.exit(unlink(paste0(in_prefix, c(".ped", ".map", ".bed", ".bim", ".fam", ".log", ".nosex")),
                 force = TRUE), add = TRUE)

  write_plink_ped_map(geno, in_prefix)

  # Convert text files to binary so plink_pairwise_ld can consume them
  run_plink_command(c("--file", in_prefix, "--make-bed", "--out", in_prefix, "--chr-set", as.character(chr_num)))

  extra_args = c("--chr-set", as.character(chr_num), extra_args)

  plink_pairwise_ld(in_prefix, ld_window = ld_window, ld_window_kb = ld_window_kb,
                    ld_window_r2 = ld_window_r2, extra_args = extra_args)
}
