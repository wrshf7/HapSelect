#######################################
##### PLINK-backed LD calculation #####
#######################################

##### Resolve the PLINK executable path #####
find_plink = function() {
  if (.Platform$OS.type == "windows") {
    install_dir = file.path(Sys.getenv("USERPROFILE", unset = ""), "bin")
    plink_path = file.path(install_dir, "plink.exe")

    if (nzchar(plink_path) && file.exists(plink_path)) {
      return(plink_path)
    }

    stop(
      "PLINK executable not found at the expected Windows location: ", plink_path, "\n",
      "Run inst/scripts/install/install_windows.ps1, or ensure PLINK is installed in that directory."
    )
  }

  plink = Sys.which("plink")
  if (nzchar(plink)) {
    return(plink)
  }

  stop("PLINK executable not found.")
}

##### Run PLINK with Windows-aware executable resolution #####
call_plink = function(args, stdout = TRUE, stderr = TRUE) {
  system2(find_plink(), args = args, stdout = stdout, stderr = stderr)
}

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

  # Look up the chromosone of each LD SNP pair
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

  # Sort by chromosone, then loci
  ld_df = ld_df[order(ld_df$Chrom, ld_df$Locus1, ld_df$Locus2), ]
  row.names(ld_df) = NULL

  return(ld_df)
}

##### Runs a plink command given a set of arguments #####
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
