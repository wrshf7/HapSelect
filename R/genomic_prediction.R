#function to check that all phenotypes have a corresponding genotype
#I did not export this function as it is an internal call. But do I need to export it
#for the package to "see" it?
check_BLUE = function(BLUE, geno){
  check = all(BLUE[,1] %in% colnames(geno)[4:ncol(geno)])
  return(check)
}

#solve the marker effects
solve_marker_effects = function(geno, BLUE){

  if(!check_BLUE(BLUE, geno)){
    stop("Ensure the first column is the exact genotype name in the genotype file and the second column is a singular adjusted phenotype, singular phenotype, BLUE, or de-regressed BLUP (no additional fixed or random effects can be fit).")
  }

  #checking to ensure the first genotype column is numeric and not the position/chromosome column.
  #checking to ensure the first column are snp marker names
  #I don't check all columns, but it's easy to update it if we should.
  if(!is.numeric(geno[,4]) | max(geno[,4], na.rm = TRUE) > 6 | min(geno[,4], na.rm = TRUE) < -6 | !is.character(geno[,1])){
    stop("Ensure that columns 1:3 of the genotype file correspond to the map (rs, chrom, pos), the SNP are characters, and columns 4:ncol(geno) are numeric values (only values between -6 and 6 are accepted.")
  }

  #only select columns that appear in the phenotype file and tranpose the matrix. Rows are individuals
  #and columns are markers (Z matrix) to connect phenos to marker effects.
  geno_mat = t(geno[, BLUE[,1]])

  #solve marker effects
  marker_solutions = rrBLUP::mixed.solve(y = BLUE[,2], Z = geno_mat)

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
create_marker_effects_file = function(geno, BLUE){
  marker_solutions = solve_marker_effects(geno, BLUE)


  prediction_accuracy = compute_prediction_accuracy(geno, marker_solutions, BLUE)
  cat(paste0("Prediction Accuracy (cor(GEBV, BLUE or deregressed BLUP or adj pheno)) is: ", round(prediction_accuracy,2)))

  #output to return
  output_df = data.frame(
    SNP = geno[,1],
    Effect = marker_solutions$u
  )
}

#perform cross-validation to assess prediction accuracy
n_fold_cross_validation = function(geno, BLUE, nfold = 5){
  if(nrow(BLUE) <= 200){
    warning("Trying to do CV on a small population size may not work well!\n")
  }
  fold_id = sample(rep(1:nfold, length.out = nrow(BLUE)))

  CV = sapply(1:nfold, function(fold){
    cat(paste0("Performing cross-validation for fold: ", fold, "\n"))
    train_BLUE = BLUE[fold_id != fold, ]
    valid_BLUE = BLUE[fold_id == fold, ]

    marker_solutions = solve_marker_effects(geno, train_BLUE)
    prediction_accuracy = compute_prediction_accuracy(geno, marker_solutions, valid_BLUE)
    return(prediction_accuracy)
  })

  mean_CV = mean(CV)

  cat(paste0(nfold, " fold CV Mean Prediction Accuracy is: ", round(mean_CV, 2), "\n"))

  return(CV)
}

cross_validation = function(geno, BLUE, train_prop = 0.9, fold = 30){
  if(nrow(BLUE) <= 200){
    warning("Trying to do CV on a small population size may not work well!\n")
  }

  CV = sapply(1:fold, function(fold){
    cat(paste0("Performing cross-validation for fold: ", fold, "\n"))
    id = sample(1:nrow(BLUE), ceiling(nrow(BLUE) * train_prop))
    train_BLUE = BLUE[id, ]
    valid_BLUE = BLUE[!(1:nrow(BLUE) %in% id), ]

    marker_solutions = solve_marker_effects(geno, train_BLUE)
    prediction_accuracy = compute_prediction_accuracy(geno, marker_solutions, valid_BLUE)
    return(prediction_accuracy)
  })

  mean_CV = mean(CV)

  cat(paste0(fold, " fold CV Mean Prediction Accuracy is: ", round(mean_CV, 2), "\n"))

  return(CV)
}
