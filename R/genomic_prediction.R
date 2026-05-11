#####################################
##### Basic Genomic Prediction ######
#####################################


#function to check that all phenotypes have a corresponding genotype
#I did not export this function as it is an internal call. But do I need to export it
#for the package to "see" it?
check_BLUE = function(BLUE, geno){
  check = all(BLUE[,1] %in% colnames(geno)[4:ncol(geno)])
  return(check)
}

#solve the marker effects
solve_marker_effects = function(geno, BLUE, h2_method, ploidy){

  if(!is.data.frame(BLUE) || ncol(BLUE) < 2){
    stop("BLUE must be a data frame with at least 2 columns: individual ID (column 1, character) and a single adjusted phenotype, BLUE, or de-regressed BLUP (column 2, numeric).")
  }
  if(!is.numeric(BLUE[,2])){
    stop("Column 2 of BLUE must be numeric (phenotype, BLUE, or de-regressed BLUP values).")
  }

  if(length(ploidy) != 1 || is.na(ploidy) || is.null(ploidy) || !is.integer(ploidy)){
    stop("Please provide a singular integer value for the ploidy. This can be specified using the #L format, where # is the ploidy number (2L is the default).")
  }

  if(!check_BLUE(BLUE, geno)){
    stop("Ensure the first column is the exact genotype name in the genotype file and the second column is a singular adjusted phenotype, singular phenotype, BLUE, or de-regressed BLUP (no additional fixed or random effects can be fit).")
  }

  #checking to ensure the first genotype column is numeric and not the position/chromosome column.
  #checking to ensure the first column are snp marker names
  #I don't check all columns, but it's easy to update it if we should.
  if(!is.numeric(geno[,4]) | max(geno[,4], na.rm = TRUE) > 20 | min(geno[,4], na.rm = TRUE) < 0 | !is.character(geno[,1])){
    stop("Ensure that columns 1:3 of the genotype file correspond to the map (SNP, Chromosome, Position), the SNP are characters, and columns 4:ncol(geno) are numeric values (only ploidy dosage values between 0 and 20 are accepted.")
  }

  #only select columns that appear in the phenotype file and tranpose the matrix. Rows are individuals
  #and columns are markers (Z matrix) to connect phenos to marker effects.
  geno_mat = t(geno[, BLUE[,1]])

  #allele frequencies
  means = colMeans(geno_mat, na.rm = TRUE)

  ploidy = max(geno_mat)

  p = means / ploidy

  # scaling for h2
  if(length(h2_method) == 1 && h2_method == "VanRaden"){

    h2_scale = ploidy * sum(p * (1-p))

  } else if(length(h2_method) == 1 && h2_method == "marker_num"){

    h2_scale = length(p)

  } else stop("Please specify a valid method to compute the total additive genetic variance: either 'VanRaden' or 'marker_num'.")

  # center
  Z = sweep(geno_mat, 2, means, "-")

  #solve marker effects
  marker_solutions = rrBLUP::mixed.solve(y = BLUE[,2], Z = Z)

  marker_variance = marker_solutions$Vu
  residual_variance = marker_solutions$Ve
  additive_genetic_variance = marker_variance * h2_scale
  h2 = additive_genetic_variance / (additive_genetic_variance + residual_variance)

  cat(sprintf("Marker Variance:           %.4e \nAdditive Genetic Variance: %.4e \nResidual Variance:         %.4e \nh2:                        %.2f \n",
              marker_variance, additive_genetic_variance, residual_variance, h2))

  return(marker_solutions)
}


#compute model prediction accuracy
compute_prediction_accuracy = function(geno, marker_effects, BLUE){


  #only select columns that appear in the phenotype file and tranpose the matrix. Rows are individuals
  #and columns are markers (Z matrix) to connect phenos to marker effects.
  geno_mat = t(geno[, BLUE[,1]])

  #create a matrix of n individuals and columns equal to nubmer of markers
  #each column is a repeated vector of a given marker's mean
  geno_means = matrix(data = 1, nrow = nrow(geno_mat), ncol = 1) %*% t(colMeans(geno_mat))

  #center the marker matrix and compute GEBV - centering not strictly needed for accuracy
  #but this is low-cost and adds prediction capabilities in the future
  GEBV = (geno_mat - geno_means) %*% marker_effects$u

  #accuracy = correlation between pheno and predicted GEBV
  predict_acc = cor(GEBV, BLUE[,2])

  return(predict_acc)
}

#obtain data frame needed for next steps
create_marker_effects_file = function(geno, BLUE, h2_method = c("VanRaden", "marker_num"), ploidy = 2L){

  h2_method = match.arg(h2_method)

  marker_solutions = solve_marker_effects(geno, BLUE, h2_method, ploidy)


  prediction_accuracy = compute_prediction_accuracy(geno, marker_solutions, BLUE)
  cat(paste0("\nPrediction Accuracy (cor(GEBV, BLUE or deregressed BLUP or adj pheno)) is: ", round(prediction_accuracy,2)))

  #output to return
  output_df = data.frame(
    SNP = geno[,1],
    Effect = marker_solutions$u
  )
}

#perform cross-validation to assess prediction accuracy
n_fold_cross_validation = function(geno, BLUE, nfold = 5L, h2_method = c("VanRaden", "marker_num"), ploidy = 2L){
  if(nrow(BLUE) <= 200){
    warning("Trying to do CV on a small population size may not work well!\n")
  }

  if(length(nfold) > 1 || is.na(nfold) || !is.integer(nfold)){
    stop("Provide a valid integer designation for the number of folds, nfold, in #L format (e.g., 5L).")
  }

  h2_method = match.arg(h2_method)

  fold_id = sample(rep(1:nfold, length.out = nrow(BLUE)))

  CV = sapply(1:nfold, function(fold){
    cat(paste0("\nPerforming cross-validation for fold: ", fold, "\n"))
    train_BLUE = BLUE[fold_id != fold, ]
    valid_BLUE = BLUE[fold_id == fold, ]

    marker_solutions = solve_marker_effects(geno, train_BLUE, h2_method, ploidy)
    prediction_accuracy = compute_prediction_accuracy(geno, marker_solutions, valid_BLUE)
    return(prediction_accuracy)
  })

  mean_CV = mean(CV)

  cat(paste0("\n",nfold, " fold CV Mean Prediction Accuracy is: ", round(mean_CV, 2), "\n"))

  return(CV)
}

cross_validation = function(geno, BLUE, train_prop = 0.9, fold = 30, h2_method = c("VanRaden", "marker_num"), ploidy = 2L){
  if(nrow(BLUE) <= 200){
    warning("Trying to do CV on a small population size may not work well!\n")
  }

  if(length(fold) > 1 || is.na(fold) || !is.integer(fold)){
    stop("Provide a valid integer designation for the number of folds, nfold, in #L format (e.g., 5L).")
  }

  h2_method = match.arg(h2_method)

  CV = sapply(1:fold, function(fold){
    cat(paste0("\nPerforming cross-validation for fold: ", fold, "\n"))
    id = sample(1:nrow(BLUE), ceiling(nrow(BLUE) * train_prop))
    train_BLUE = BLUE[id, ]
    valid_BLUE = BLUE[!(1:nrow(BLUE) %in% id), ]

    marker_solutions = solve_marker_effects(geno, train_BLUE, h2_method, ploidy)
    prediction_accuracy = compute_prediction_accuracy(geno, marker_solutions, valid_BLUE)
    return(prediction_accuracy)
  })

  mean_CV = mean(CV)

  cat(paste0("\n",fold, " fold CV Mean Prediction Accuracy is: ", round(mean_CV, 2), "\n"))

  return(CV)
}
