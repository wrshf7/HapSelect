########################################
######### Native LD Calculation ########
########################################

# ld_func ----------------------------------------------------------------------
# Calculates r^2 between marker pairs on one chromosome, in R. The compiled
# equivalent is ld_func_c(); pairwise_ld() and pairwise_ld_c() call these once per
# chromosome.
#
# genotypes : one chromosome's markers, with row names set to the marker names:
#             SNP ID, chromosome and position in columns 1 to 3, one dosage column
#             per individual after that
# window,
# min_r2,
# min_obs   : as documented on pairwise_ld()
ld_func = function(genotypes, window = NULL, min_r2 = NULL, min_obs = 2L){
  #the columns returned when no pair survives
  empty = data.frame(Chrom = character(), Locus1 = integer(), Locus2 = integer(),
                     Name1 = character(), Name2 = character(), LD = numeric(),
                     stringsAsFactors = FALSE)

  #pull the chromosome number
  chromo = genotypes[1,2]

  #pull the marker names of the chromosome
  marker_names = row.names(genotypes)

  #remove chromo and snp name columns, transpose the df so that markers are columns
  genotypes = genotypes[,-(1:3)]
  genotypes = as.data.frame(t(as.matrix(genotypes)))

  n_marker = length(marker_names)
  if (is.null(window)) window = n_marker

  #iterate over marker 1 through marker n-1 (all but the last marker) and row bind the output
  ld_df = map_dfr(1:(n_marker-1), function(i){

    #iterate over the remaining markers from i+1 to n, or as far as the window reaches
    j_end = min(n_marker, i + window)
    if (j_end <= i) return(NULL)

    marker_ld = map_dfr((i+1):j_end, function(j){
      #extract marker names of the current 2 markers being compared
      snp1 = marker_names[i]
      snp2 = marker_names[j]

      #compute r^2 of the current two markers
      snp_cor2 = suppressWarnings(
        cor(genotypes[,i], genotypes[,j], use = "pairwise.complete.obs")^2
      )

      # cor() is already NA below two shared individuals, so counting them is only
      # needed when a stricter minimum was asked for
      if (min_obs > 2L && !is.na(snp_cor2)) {
        n_obs = sum(!is.na(genotypes[,i]) & !is.na(genotypes[,j]))
        if (n_obs < min_obs) snp_cor2 = NA_real_
      }

      # an undefined r^2 says nothing about the pair, so it is never reported
      if (is.na(snp_cor2)) return(NULL)
      if (!is.null(min_r2) && snp_cor2 < min_r2) return(NULL)

      #give the relevant info to the markers being compared in a df
      return_df = data.frame(
        Chrom = chromo,
        Locus1 = i,
        Locus2 = j,
        Name1 = snp1,
        Name2 = snp2,
        LD = snp_cor2
      )

      #return the relevant info and rbind it for the ith iteration
      return(return_df)
    })

    #return the data frame from the inner loop and rbind it
    return(marker_ld)
  })

  if (nrow(ld_df) == 0) return(empty)

  ld_df
}

# ld_func_c --------------------------------------------------------------------
# Wrapper around the C++ implementation of ld_func. Accepts the same arguments and
# returns the same result, but delegates to compiled code, which is much faster
# than building one data frame per marker pair in R.
#
# genotypes : as for ld_func(), though the marker names are read from column 1
#             rather than the row names
# window,
# min_r2,
# min_obs   : as documented on pairwise_ld()
ld_func_c = function(genotypes, window = NULL, min_r2 = NULL, min_obs = 2L){
  chromo = genotypes[1,2]
  marker_names = as.character(genotypes[,1])

  #markers as columns, individuals as rows, which is what the C++ side reads
  dosages = t(as.matrix(genotypes[,-(1:3)]))
  storage.mode(dosages) = "double"

  empty = data.frame(Chrom = character(), Locus1 = integer(), Locus2 = integer(),
                     Name1 = character(), Name2 = character(), LD = numeric(),
                     stringsAsFactors = FALSE)

  n_marker = length(marker_names)
  if (n_marker < 2) return(empty)

  pairs = pairwise_ld_cpp(
    geno    = dosages,
    window  = as.integer(if (is.null(window)) n_marker else window),
    #a negative floor keeps everything, since r^2 is never below zero
    min_r2  = as.numeric(if (is.null(min_r2)) -1 else min_r2),
    min_obs = as.integer(min_obs)
  )

  #a chromosome where nothing cleared the filters still has to return the columns
  if (nrow(pairs) == 0) return(empty)

  data.frame(
    Chrom  = chromo,
    Locus1 = pairs$Locus1,
    Locus2 = pairs$Locus2,
    Name1  = marker_names[pairs$Locus1],
    Name2  = marker_names[pairs$Locus2],
    LD     = pairs$LD,
    stringsAsFactors = FALSE
  )
}

# pairwise_ld_run --------------------------------------------------------------
# The chromosome splitting, parallelisation and progress reporting shared by
# pairwise_ld() and pairwise_ld_c(). Only the per-chromosome function differs.
#
# chromosome_fn : ld_func() for the R implementation, ld_func_c() for the compiled
#                 one
# the remaining arguments are documented on pairwise_ld()
pairwise_ld_run = function(genotype_matrix, parallelize, window, min_r2, min_obs,
                           chromosome_fn){
  # Validate the input genotype matrix structure and content before proceeding with LD calculations
  check_ld_matrix(genotype_matrix)

  #set row names to the marker names for the internal loop
  row.names(genotype_matrix) = genotype_matrix[,1]

  #split the genotype matrix by chromosome for parallelization
  genotype_matrix = split(genotype_matrix, genotype_matrix[,2])

  #setup parallelization using future and parallel package and utilize all but 1 core
  if(parallelize){
    future::plan(multisession, workers = parallel::detectCores() - 1)
    on.exit(future::plan(sequential), add = TRUE)
  }

  #setup progress bar
  handlers("txtprogressbar")

  chromosome_ld = function(genotypes){
    chromo_ld = chromosome_fn(genotypes, window = window, min_r2 = min_r2,
                              min_obs = min_obs)
    p()
    chromo_ld
  }

  #call progress bar and perform main function
  with_progress({

    #define the progress bar - it has as many iterations (along) as the list provided
    p = progressor(along = genotype_matrix)

    #parallelize the different chromosomes with furrr, provide their genotype data frames, and row bind all chromosomes back together
    #seed = TRUE gives each worker a parallel-safe RNG stream; LD draws no random
    #numbers, but loading packages on a worker can, and future warns when it does
    all_ld = if(parallelize){
      furrr::future_map_dfr(genotype_matrix, chromosome_ld,
                            .options = furrr::furrr_options(seed = TRUE))
    } else {
      purrr::map_dfr(genotype_matrix, chromosome_ld)
    }

  })

  #return the entire dataframe with pairwise marker LDs within chromosome
  return(all_ld)
}

# pairwise_ld ------------------------------------------------------------------
# Calculates r^2 between pairs of markers, within each chromosome. Markers on
# different chromosomes are never compared. pairwise_ld_c() does the same thing in
# compiled code and is much faster on anything but a small dataset.
#
# genotype_matrix : marker map and dosages in one data frame, markers as rows:
#                   SNP ID (character) in column 1, chromosome (numeric) in column
#                   2, position (numeric) in column 3, and one column per
#                   individual after that holding dosages of 0, 1, 2 or NA
# parallelize     : if TRUE, process chromosomes in parallel using all available
#                   cores minus one
# window          : only compare markers at most this many positions apart within
#                   the chromosome. NULL compares every pair, which is quadratic
#                   in the marker count and so the expensive choice on dense data
# min_r2          : drop pairs whose r^2 falls below this. NULL keeps every pair
# min_obs         : individuals that must be observed at both markers before an
#                   r^2 is defined. Two individuals always lie on a line, so a
#                   pair sharing only two gives r^2 = 1 whatever the genotypes
#                   were; raising this to 3 refuses those, at the cost of losing
#                   pairs on sparsely observed markers
#
# Returns a data frame with one row per marker pair and columns:
#   Chrom          chromosome the pair sits on
#   Locus1, Locus2 the two markers' positions within that chromosome, counting
#                  from 1 in the order they appear
#   Name1, Name2   the two markers' SNP IDs
#   LD             r^2, over the individuals observed at both markers
#
# A pair whose r^2 is undefined - one of its markers monomorphic, or too few
# individuals observed at both - is left out rather than returned with LD = NA.
pairwise_ld = function(genotype_matrix, parallelize = TRUE, window = NULL,
                       min_r2 = NULL, min_obs = 2L){
  pairwise_ld_run(genotype_matrix, parallelize = parallelize, window = window,
                  min_r2 = min_r2, min_obs = min_obs,
                  chromosome_fn = ld_func)
}

# Marker count per chromosome above which an unwindowed run is usually worth
# parallelising. Starting workers takes time that the work has to earn back, and
# a window keeps the pair count linear in the marker count while leaving it out
# makes it quadratic - so the two kinds of run sit on opposite sides of that.
# A rough threshold rather than a precise one, since machines differ.
LD_PARALLEL_MARKER_HINT = 1000

# advise_parallel_ld -----------------------------------------------------------
# Messages when a serial run looks large enough that parallelising would pay, so
# that a long job does not run serially just because that is the default.
#
# genotype_matrix : as passed to pairwise_ld_c()
# window          : forward marker window, or NULL for every pair
advise_parallel_ld = function(genotype_matrix, window){
  # Leave a malformed input to check_ld_matrix(), which reports it properly
  if (!is.data.frame(genotype_matrix) || ncol(genotype_matrix) < 4) return(invisible(NULL))

  # A window is the case parallelising does not pay for, whatever the marker count
  if (!is.null(window)) return(invisible(NULL))

  per_chromosome = as.vector(table(genotype_matrix[[2]]))

  # Work is split by chromosome, so a single one cannot be spread over workers
  if (length(per_chromosome) < 2) return(invisible(NULL))
  if (max(per_chromosome) <= LD_PARALLEL_MARKER_HINT) return(invisible(NULL))

  message("pairwise_ld_c() is running serially over every within-chromosome pair, ",
          "with up to\n  ", format(max(per_chromosome), big.mark = ","),
          " markers on a chromosome. Without a window the pair count grows with ",
          "the\n  square of the markers, and parallelising is usually worth it ",
          "from somewhere\n  around ", format(LD_PARALLEL_MARKER_HINT, big.mark = ","),
          " markers per chromosome upwards. Pass parallelize = TRUE to try\n  it, ",
          "or parallelize = FALSE to keep it serial without this message.")
  invisible(NULL)
}

# pairwise_ld_c ----------------------------------------------------------------
# Calculates r^2 between pairs of markers within each chromosome, in compiled
# code. Takes the same arguments as pairwise_ld() and returns the same result,
# but is much faster, so it is the one to reach for on anything but a small
# dataset.
#
# genotype_matrix,
# window,
# min_r2,
# min_obs         : as documented on pairwise_ld()
# parallelize     : if TRUE, process chromosomes in parallel using all available
#                   cores minus one. Defaults to FALSE, unlike pairwise_ld():
#                   starting workers costs time that the R implementation earns
#                   back easily but the compiled one often does not, being fast
#                   enough already. A run large enough to benefit says so
#
# Returns the same columns as pairwise_ld(): Chrom, Locus1, Locus2, Name1, Name2
# and LD.
pairwise_ld_c = function(genotype_matrix, parallelize = FALSE, window = NULL,
                         min_r2 = NULL, min_obs = 2L){
  # Only advise when the default was left in place, not when serial was chosen
  if (missing(parallelize) && !parallelize) {
    advise_parallel_ld(genotype_matrix, window)
  }

  pairwise_ld_run(genotype_matrix, parallelize = parallelize, window = window,
                  min_r2 = min_r2, min_obs = min_obs,
                  chromosome_fn = ld_func_c)
}

# check_ld_matrix --------------------------------------------------------------
# Stops unless the genotype matrix has the structure pairwise_ld() documents: SNP
# IDs, chromosome and position in columns 1 to 3, numeric dosages after that.
#
# genotype_matrix : the data frame to check
check_ld_matrix = function(genotype_matrix) {
  # Check column structure and types of the genotype matrix
  if(!is.data.frame(genotype_matrix) || ncol(genotype_matrix) < 4){
    stop("genotype_matrix must be a data frame with at least 4 columns: SNP ID (character), Chromosome (numeric), Position (numeric), and one or more individual dosage columns (numeric integer, values 0/1/2/3+/NA).")
  }
  if(!is.character(genotype_matrix[,1])){
    stop("Column 1 of genotype_matrix must be character SNP IDs.")
  }
  if(!is.numeric(genotype_matrix[,2])){
    stop("Column 2 of genotype_matrix must be numeric chromosome identifiers.")
  }
  if(!is.numeric(genotype_matrix[,3])){
    stop("Column 3 of genotype_matrix must be numeric marker positions (physical or genetic map position).")
  }
  if(!is.numeric(genotype_matrix[,4])){
    stop("Columns 4 and onward of genotype_matrix must be numeric dosage values (0, 1, 2, 3+, or NA).")
  }
}
