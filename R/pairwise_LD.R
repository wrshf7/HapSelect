########################################
######### Native LD Calculation ########
########################################

#####ld_func() - function to calculate LD within a chromosome utilizing the subsetted marker genotype file####
#file str is row names are marker names, first column is the marker name, second column is the chromosome, and the rest are genotyped individuals
ld_func = function(genotypes){
  #pull the chromosome number
  chromo = genotypes[1,2]

  #pull the marker names of the chromosome
  marker_names = row.names(genotypes)

  #remove chromo and snp name columns, transpose the df so that markers are columns
  genotypes = genotypes[,-(1:3)]
  genotypes = as.data.frame(t(as.matrix(genotypes)))

  #iterate over marker 1 through marker n-1 (all but the last marker) and row bind the output
  ld_df = map_dfr(1:(length(marker_names)-1), function(i){

    #iterate over the remaining markers from i+1 to n
    marker_ld = map_dfr((i+1):length(marker_names), function(j){
      #j = as.numeric(j)

      #extract marker names of the current 2 markers being compared
      snp1 = marker_names[i]
      snp2 = marker_names[j]

      #compute r^2 of the current two markers
      snp_cor2 = cor(genotypes[,i],genotypes[,j], use = "pairwise.complete.obs")^2

      #give the relevant info to the markers being compared in a df
      return_df = data.frame(
        Chrom = chromo,
        Locus1 = i,
        Locus2 = j,
        Name1 = snp1,
        Name2 = snp2,
        LD = snp_cor2
      )

      #return_df = c(chromo, i, j, snp1, snp2, snp_cor2)

      #return the relevant info and rbind it for the ith iteration
      return(return_df)
    })

    #return the data frame from the inner loop and rbind it
    return(marker_ld)
  })

}

#####pairwise_ld() - function to calculate LD across all chromos and parallelize chromos - calls ld_func()####
#requires a genotype matrix with SNP name as the first column, chromosome identifier as the second column, markers as rows, and individuals genotyped as columns 3 onwards
#genotypes should be dosages: 0,1,2
pairwise_ld = function(genotype_matrix, parallelize = TRUE){
  # Validate the input genotype matrix structure and content before proceeding with LD calculations
  check_ld_matrix(genotype_matrix)

  #set row names to the marker names for the internal loop
  row.names(genotype_matrix) = genotype_matrix[,1]

  #split the genotype matrix by chromosome for parallelization
  genotype_matrix = split(genotype_matrix, genotype_matrix[,2])

  #setup parallelization using future and parallel package and utilize all but 1 core
  if(parallelize){
    future::plan(multisession, workers = parallel::detectCores() - 1)
  }

  #setup progress bar
  handlers("txtprogressbar")

  #call progress bar and perform main function
  with_progress({

    #define the progress bar - it has as many iterations (along) as the list provided
    p = progressor(along = genotype_matrix)

    #parallelize the different chromosomes with furrr, provide their genotype data frames, and row bind all chromosomes back together
    if(parallelize){
      all_ld = furrr::future_map_dfr(genotype_matrix, function(genotypes){

        #call the ld_func() on the chromosome
        chromo_ld = ld_func(genotypes)

        #once it's done, progress the progress bar and return the chromo data frame
        p()
        return(chromo_ld)
      })
    } else {
      all_ld = purrr::map_dfr(genotype_matrix, function(genotypes){

        #call the ld_func() on the chromosome
        chromo_ld = ld_func(genotypes)

        #once it's done, progress the progress bar and return the chromo data frame
        p()
        return(chromo_ld)
      })
    }

  })

  #release all of the cores and undo parallelization
  if(parallelize){
    future::plan(sequential)
  }


  #bind all of the chromosomes tog
  #all_ld = bind_rows(all_ld)

  #return the entire dataframe with pairwise marker LDs within chromosome
  return(all_ld)
}

# Checks the ld matrix and stops if it does not meet the expected structure and content for the pairwise_ld function. 
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
