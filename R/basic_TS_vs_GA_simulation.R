#########################################
####### Basic GA vs TS Simulation #######
#########################################

# main entry function to do a basic GA vs TS simulation
# GA_output: the output from the genetic algorithm in the package
# geno: the input used throughout the pipeline
# marker_effects: the input used in the localGEBV step
# map: the input used for haploblocking and LD calculations
# genetic_map_position: an optional argument. The user should supply a vector of the map distances for each marker.
#   note: this should be in the same order as the map file!
# num_sim_reps: the number of times to repeat the simulation to quantify Monte Carlo variability/error
# num_TS_parents: the number of parents selected via truncation selection to compare against.
#   note: the default, NULL, uses the same number as the GA parents supplied
# mean_adjust: the same parameter as the localGEBV step. It should be kept as TRUE in almost all circumstances.
# max_cM_chr: default value of 100 used per chromosome. This argument is only used if genetic_map_position is not supplied.
# colors: the colors of the GA selected, TS selected, GA/TS selected overlapping parents, and all other individuals

# function to check compatability
check_geno_marker_compatibility = function(geno, marker_effects, map){

  # Extract marker IDs
  geno_ids   = geno[[1]]
  effect_ids = marker_effects$SNP
  map_ids    = map$SNP

  # Type coercion
  if(!is.character(geno_ids)){
    geno_ids = as.character(geno_ids)
    warning("Coerced first column of geno to character.")
  }

  if(!is.character(effect_ids)){
    effect_ids = as.character(effect_ids)
    warning("Coerced marker_effects$SNP to character.")
  }

  if(!is.character(map_ids)){
    map_ids = as.character(map_ids)
    warning("Coerced map$SNP to character.")
  }

  # Missing values
  if(anyNA(geno_ids)){
    stop("Missing values found in geno marker IDs.")
  }

  if(anyNA(effect_ids)){
    stop("Missing values found in marker_effects marker IDs.")
  }

  if(anyNA(map_ids)){
    stop("Missing values found in map marker IDs.")
  }

  # Duplicates
  if(anyDuplicated(geno_ids)){
    stop("Duplicate marker IDs found in geno.")
  }

  if(anyDuplicated(effect_ids)){
    stop("Duplicate marker IDs found in marker_effects.")
  }

  if(anyDuplicated(map_ids)){
    stop("Duplicate marker IDs found in map.")
  }

  # Set equality
  if(!setequal(geno_ids, effect_ids) || !setequal(geno_ids, map_ids)){
    stop("Mismatch between geno, marker_effects, and map marker sets.")
  }
}


#write out the simulation files
write_sim_files = function(GA_geno, TS_geno, map, marker_effects, temp_string){
  file_list = c(
    map_file_name = paste0("Genomic_Simulation_genetic_map_temp_",temp_string,".txt"),
    geno_GA_file_name = paste0("Genomic_Simulation_GA_genotypes_temp_",temp_string,".txt"),
    geno_TS_file_name = paste0("Genomic_Simulation_TS_genotypes_temp_",temp_string,".txt"),
    marker_eff_file_name = paste0("Genomic_Simulation_marker_effects_temp_",temp_string,".txt")
  )

  write.table(map, file = file_list["map_file_name"], sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(GA_geno, file = file_list["geno_GA_file_name"], sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(TS_geno, file = file_list["geno_TS_file_name"], sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(marker_effects, file = file_list["marker_eff_file_name"], sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

  return(file_list)
}

#clean up the simulation files
delete_sim_files = function(file_list){
  invisible(file.remove(file_list))
}

# simulation function
run_basic_simulation = function(temp_files, num_gen, num_sim_reps,
                                num_cross_per_gen, num_GA_parents, num_TS_parents,
                                maximize){

  #tracking
  mean_BV_GA = list()
  mean_BV_TS = list()


  #run it for the GA parents
  cat("\nRunning simulation for GA selected parents:\n\n")

  #initiate populations
  init = genomicSimulation::load.data(allele.file = temp_files[[2]],
                                      map.file = temp_files[[1]],
                                      effect.file = temp_files[[4]])
  g0 = init$groupNum

  #simulate across reps
  for(rep in seq_len(num_sim_reps)){
    current_pop = genomicSimulation::make.random.crosses(g0, n.crosses = num_cross_per_gen, give.names = FALSE, name.prefix = "F1.")
    mean_BV_GA[[rep]] = sim_gens(num_gen = num_gen, init = init,
                                 num_parents = num_GA_parents, num_cross_per_gen = num_cross_per_gen,
                                 current_pop = current_pop, maximize = maximize)
  }

  groups = genomicSimulation::see.existing.groups()$Group
  for(group in groups){
    invisible(utils::capture.output(genomicSimulation::delete.group(group)))

  }

  #run it for the TS parents
  cat("\nRunning simulation for TS selected parents:\n\n")

  #initiate populations
  init = genomicSimulation::load.data(allele.file = temp_files[[3]],
                                      map.file = temp_files[[1]],
                                      effect.file = temp_files[[4]])
  g0 = init$groupNum

  for(rep in seq_len(num_sim_reps)){
    current_pop = genomicSimulation::make.random.crosses(g0, n.crosses = num_cross_per_gen, give.names = FALSE, name.prefix = "F1.")
    mean_BV_TS[[rep]] = sim_gens(num_gen = num_gen, init = init,
                                 num_parents = num_TS_parents, num_cross_per_gen = num_cross_per_gen,
                                 current_pop = current_pop, maximize = maximize)
  }

  groups = genomicSimulation::see.existing.groups()$Group
  for(group in groups){
    invisible(utils::capture.output(genomicSimulation::delete.group(group)))

  }

  # GA
  df_GA = do.call(rbind, lapply(seq_along(mean_BV_GA), function(i) {
    data.frame(
      rep = i,
      gen = seq_along(mean_BV_GA[[i]]),
      value = mean_BV_GA[[i]],
      method = "GA"
    )
  }))

  # TS
  df_TS = do.call(rbind, lapply(seq_along(mean_BV_TS), function(i) {
    data.frame(
      rep = i,
      gen = seq_along(mean_BV_TS[[i]]),
      value = mean_BV_TS[[i]],
      method = "TS"
    )
  }))

  df_all = rbind(df_GA, df_TS)

  df_summary = df_all %>%
    dplyr::group_by(method, gen) %>%
    summarise(
      mean = mean(value),
      se = sd(value) / sqrt(n()),
      .groups = "drop"
    )
  return(df_summary)
}

# main function to run the simulation for specified number of generations
sim_gens = function(num_gen, init, num_parents, num_cross_per_gen, current_pop, maximize){
  mean_BV = c()
  for(gen in seq_len(num_gen)){
    # Record mean BV of current population
    mean_BV[gen] = mean(genomicSimulation::see.GEBVs(current_pop))

    # Select top parents by GEBV
    selected = genomicSimulation::break.group.by.GEBV(
      current_pop,
      number = num_parents,
      #if maximizing is TRUE, then selecting by low should be FALSE
      #if maximize is FALSE, then select by low should be TRUE
      low.score.best = !maximize
    )

    # Generate next generation (1000 offspring)
    current_pop = genomicSimulation::make.random.crosses(
      group       = selected,
      n.crosses   = num_cross_per_gen,
      give.names  = TRUE,
      name.prefix = paste0("F", gen + 1, ".")
    )

    invisible(utils::capture.output(genomicSimulation::delete.group(selected)))
  }

  invisible(utils::capture.output(genomicSimulation::delete.group(current_pop)))

  return(mean_BV)
}

# Simulation plotting function
generate_sim_traject_plot = function(summary_df, colors){
  traject_plot = ggplot2::ggplot(summary_df, aes(x = gen, y = mean, color = method, fill = method)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se, group = method), alpha = 0.2, show.legend = FALSE) +
    labs(x = "Generation", y = "Mean Breeding Value") + theme_cowplot() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors, guide = "none")

  return(traject_plot)
}


# PCA parent selection plot
generate_selection_pca_plot = function(genotype_matrix_centered,
                                       GA_parents_indices,
                                       TS_indices,
                                       colors, alpha){

  # PCA
  pca = stats::prcomp(
    genotype_matrix_centered,
    center = TRUE,
    scale. = TRUE
  )

  # Extract first 10 PCs
  pc_df = as.data.frame(pca$x[, 1:10])

  # Rename columns explicitly
  colnames(pc_df) = paste0("PC", 1:10)

  # Combine with accession names and grouping
  pca_df = data.frame(
    accession = rownames(genotype_matrix_centered),
    pc_df,
    group = "Not Selected"
  )

  # Assign groups
  pca_df$group[GA_parents_indices] = "GA"
  pca_df$group[TS_indices] = "TS"

  # Overlap
  overlap = dplyr::intersect(GA_parents_indices, TS_indices)

  if(length(overlap) > 0){
    pca_df$group[overlap] = "Both"
  }

  # Factor ordering
  pca_df$group = factor(
    pca_df$group,
    levels = c("Not Selected", "TS", "GA", "Both")
  )

  # Colors
  plot_colors = c(
    `Not Selected` = colors[4],
    TS = colors[2],
    GA = colors[1],
    Both = colors[3]
  )

  # Alpha
  plot_alphas = c(
    `Not Selected` = alpha[4],
    TS = alpha[2],
    GA = alpha[1],
    Both = alpha[3]
  )

  # Size
  plot_sizes = c(
    `Not Selected` = 1,
    TS = 2.5,
    GA = 2.5,
    Both = 3
  )

  # Variance explained
  var_exp = summary(pca)$importance[2,1:2] * 100

  # Plot
  pca_plot = ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, y = PC2, color = group, alpha = group, size = group)
  ) + theme_cowplot() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
      color = "Selection", alpha = "Selection", size = "Selection"
    ) +
    ggplot2::scale_color_manual(values = plot_colors) +
    ggplot2::scale_alpha_manual(values = plot_alphas) +
    ggplot2::scale_size_manual( values = plot_sizes)

  return_list = list(PCA_df = pca_df, PCA_plot = pca_plot)

  return(return_list)
}

encode_phased_haplotypes = function(hap1, hap2){

  if(length(hap1) != length(hap2)){
    stop("Haplotype columns must have same length.")
  }

  # ---- handle missing phase BEFORE encoding ----
  # If either hap is NA, randomly assign a valid phased state
  missing <- is.na(hap1) | is.na(hap2)

  if(any(missing)){
    n_missing <- sum(missing)

    # random 0/1 pairs but enforcing heterozygous dosage = 1
    # (since missing genotype still assumed known dosage)
    rand_phase <- sample(c("AT", "TA"), n_missing, replace = TRUE)

    hap1[missing] <- ifelse(rand_phase == "AT", 1, 0)
    hap2[missing] <- ifelse(rand_phase == "AT", 0, 1)
  }

  out <- ifelse(
    hap1 == 0 & hap2 == 0, "TT",
    ifelse(
      hap1 == 0 & hap2 == 1, "TA",
      ifelse(
        hap1 == 1 & hap2 == 0, "AT",
        ifelse(
          hap1 == 1 & hap2 == 1, "AA", NA
        )
      )
    )
  )

  return(out)
}

prepare_phased_inputs = function(geno_phased){

  snp_info <- geno_phased[,1:3]

  parent_names <- unique(sub("_[12]$", "", setdiff(colnames(geno_phased), names(snp_info))))

  phased_matrix <- lapply(parent_names, function(p){

    h1 <- geno_phased[[paste0(p, "_1")]]
    h2 <- geno_phased[[paste0(p, "_2")]]

    encode_phased_haplotypes(h1, h2)
  })

  phased_matrix <- do.call(cbind, phased_matrix)
  colnames(phased_matrix) <- parent_names

  return(phased_matrix)
}


prepare_dosage_inputs = function(geno){

  individual_names <- geno[,1]

  genotype_matrix <- t(geno[,4:ncol(geno)])
  colnames(genotype_matrix) <- individual_names

  # missing -> heterozygous placeholder
  genotype_matrix[is.na(genotype_matrix)] <- 1

  return(genotype_matrix)
}

prepare_phased_dosage_inputs = function(geno_phased){

  snp_info <- geno_phased[,1:3]

  parent_names <- unique(
    sub("_[12]$", "",
        setdiff(colnames(geno_phased), names(snp_info)))
  )

  dosage_matrix <- lapply(parent_names, function(p){

    h1 <- geno_phased[[paste0(p, "_1")]]
    h2 <- geno_phased[[paste0(p, "_2")]]

    dosage <- h1 + h2

    # preserve missing values
    dosage[is.na(h1) | is.na(h2)] <- NA

    dosage
  })

  dosage_matrix <- do.call(cbind, dosage_matrix)

  colnames(dosage_matrix) <- parent_names
  rownames(dosage_matrix) <- geno_phased[[1]]

  return(dosage_matrix)
}

localGEBV_vs_TS_simulation = function(
    GA_output,
    geno,
    marker_effects,
    map,
    genetic_map_position = NULL,
    num_gen = 50,
    num_sim_reps = 30,
    num_cross_per_gen = 1000,
    num_TS_parents = NULL,
    mean_adjust = TRUE,
    maximize = TRUE,
    max_cM_chr = 100,
    PCA = TRUE,
    colors = c("green", "#d95f02", "#A01FF0", "gray80"),
    alpha = c(1,1,1,0.5)
){

  check_geno_marker_compatibility(geno, marker_effects, map)

  #check the number of crosses per generation
  if(length(num_cross_per_gen) >1  || anyNA(as.numeric(num_cross_per_gen))) stop("The num_cross_per_gen value must either be NULL or a single numeric and cannot be NA!")

  #check the number of crosses per generation
  if(length(num_gen) >1  || anyNA(as.numeric(num_gen))) stop("The num_gen value must be a single numeric and cannot be NA!")

  #check maximize is specified correctly
  if(length(maximize) != 1 || any(is.na(maximize)) || any(is.null(maximize)) || any(!is.logical(maximize))){
    stop("Please ensure the maximize argument is only specified as TRUE or FALSE.")
  }

  #check number of colors
  if(
    !.check_color(x = colors, n = 4, allow_na = FALSE)
  ){
    stop("Please provide 4 valid R colors for plotting purposes. The first two are for GA and TS selected parents, respectively, and the latter two are for the overlap and non-selected parents on the PCA plot.")
  }

  #check alphas
  if(
    length(alpha) != 4 ||
    anyNA(alpha) ||
    !is.numeric(alpha) ||
    any(alpha < 0 | alpha > 1)
  ){
    stop("alpha must contain 4 numeric values between 0 and 1.")
  }

  GA_parents_names <- GA_output$selected_founders$individuals
  GA_parents_indices <- GA_output$selected_founders$indices

  num_GA_parents <- length(GA_parents_indices)

  if(is.null(num_TS_parents)){
    num_TS_parents <- num_GA_parents
  }

  snp_names <- geno[,1]

  genotype_matrix <- prepare_dosage_inputs(geno)

  if(mean_adjust){
    genotype_matrix_centered <- center_genotypes(geno)
    genotype_matrix_centered <- t(genotype_matrix_centered[,4:ncol(genotype_matrix_centered)])
    genotype_matrix_centered[is.na(genotype_matrix_centered)] = 0
  } else {
    genotype_matrix_centered <- t(genotype_matrix[,4:ncol(genotype_matrix)])
    marker_means                          = colMeans(genotype_matrix_centered, na.rm = TRUE)
    col_indices                           = which(is.na(genotype_matrix_centered), arr.ind = TRUE)
    genotype_matrix_centered[col_indices] = marker_means[genotype_matrix_centered[, 2]]
  }

  colnames(genotype_matrix_centered) <- snp_names

  marker_effects <- marker_effects[match(colnames(genotype_matrix_centered), marker_effects$SNP), ]

  GEBV <- genotype_matrix_centered %*% marker_effects[,2]
  GEBV <- cbind(GEBV, 1:nrow(GEBV))

  #sort based on selecting for min or max
  if(maximize){
    GEBV <- GEBV[order(GEBV[,1], decreasing=TRUE), ]
  } else{
    GEBV <- GEBV[order(GEBV[,1], decreasing=FALSE), ]
  }

  TS_indices <- GEBV[1:num_TS_parents,2]

  genotype_matrix[is.na(genotype_matrix)] <- 1

  GA_geno <- genotype_matrix[GA_parents_indices, ]
  TS_geno <- genotype_matrix[TS_indices, ]
  GA_geno <- data.frame(ID = rownames(GA_geno), GA_geno)
  TS_geno <- data.frame(ID = rownames(TS_geno), TS_geno)

  # map handling unchanged
  if(is.null(genetic_map_position)){
    map <- map %>%
      dplyr::group_by(Chromosome) %>%
      dplyr::mutate(cM = max_cM_chr*(Position-min(Position))/(max(Position)-min(Position))) %>%
      dplyr::ungroup()
  } else {
    map$cM <- genetic_map_position
  }

  map <- map[,c("SNP","Chromosome","cM")]
  map <- data.frame(chr=map$Chromosome, marker=map$SNP, pos=map$cM)

  colnames(marker_effects) <- c("marker","eff")
  marker_effects$allele <- "A"
  marker_effects <- marker_effects[,c(1,3,2)]

  temp_string <- paste0(sample(c(letters,LETTERS,0:9),12,replace=TRUE), collapse="")

  temp_files <- write_sim_files(
    GA_geno = GA_geno,
    TS_geno = TS_geno,
    map = map,
    marker_effects = marker_effects,
    temp_string = temp_string
  )

  summary_df <- run_basic_simulation(
    temp_files,
    num_gen,
    num_sim_reps,
    num_cross_per_gen,
    num_GA_parents,
    num_TS_parents,
    maximize = maximize
  )

  delete_sim_files(temp_files)

  gc()

  sim_plot <- generate_sim_traject_plot(summary_df, colors[1:2])

  if(PCA){
    PCA_list <- generate_selection_pca_plot(
      genotype_matrix_centered,
      GA_parents_indices,
      TS_indices,
      colors,
      alpha
    )
  }


  if(PCA){
    return(list(
      Simulation_Plot = sim_plot,
      PCA_Plot = PCA_list$PCA_plot,
      Simulation_Summary = summary_df,
      PCA_df = PCA_list$PCA_df
    ))
  } else {
    return(list(Simulation_Plot = sim_plot,
           Simulation_Summary = summary_df))
  }
}


Haplotype_vs_TS_simulation = function(
    GA_output,
    geno_phased,
    marker_effects,
    map,
    genetic_map_position = NULL,
    num_gen = 50,
    num_sim_reps = 30,
    num_cross_per_gen = 1000,
    num_TS_parents = NULL,
    mean_adjust = TRUE,
    maximize = TRUE,
    max_cM_chr = 100,
    PCA = TRUE,
    colors = c("green", "#d95f02", "#A01FF0", "gray80"),
    alpha = c(1,1,1,0.5)
){

  check_geno_marker_compatibility(geno_phased, marker_effects, map)

  #check the number of crosses per generation
  if(length(num_cross_per_gen) >1  || anyNA(as.numeric(num_cross_per_gen))) stop("The num_cross_per_gen value must either be NULL or a single numeric and cannot be NA!")

  #check the number of crosses per generation
  if(length(num_gen) >1  || anyNA(as.numeric(num_gen))) stop("The num_gen value must be a single numeric and cannot be NA!")

  #check maximize is specified correctly
  if(length(maximize) != 1 || any(is.na(maximize)) || any(is.null(maximize)) || any(!is.logical(maximize))){
    stop("Please ensure the maximize argument is only specified as TRUE or FALSE.")
  }

  #check number of colors
  if(
    !.check_color(x = colors, n = 4, allow_na = FALSE)
  ){
    stop("Please provide 4 valid R colors for plotting purposes. The first two are for GA and TS selected parents, respectively, and the latter two are for the overlap and non-selected parents on the PCA plot.")
  }

  #check alphas
  if(
    length(alpha) != 4 ||
    anyNA(alpha) ||
    !is.numeric(alpha) ||
    any(alpha < 0 | alpha > 1)
  ){
    stop("alpha must contain 4 numeric values between 0 and 1.")
  }


  # ---- GA parents now by NAME ----
  GA_parent_names <- GA_output$selected_founders$individuals

  if(is.null(num_TS_parents)){
    num_TS_parents <- length(GA_parent_names)
  }

  num_GA_parents = length(GA_parent_names)

  phased_matrix <- prepare_phased_inputs(geno_phased)

  dosage_matrix <- prepare_phased_dosage_inputs(geno_phased)

  # # ensure alignment safety
  # phased_matrix <- phased_matrix[, GA_parent_names %in% colnames(phased_matrix)]

  individual_names <- colnames(phased_matrix)


  colnames(dosage_matrix) <- colnames(phased_matrix)
  rownames(dosage_matrix) <- rownames(phased_matrix)

  dosage_matrix = cbind(geno_phased[,1:3], dosage_matrix)

  if(mean_adjust){
    genotype_matrix_centered <- center_genotypes(dosage_matrix)
    genotype_matrix_centered <- t(genotype_matrix_centered[,4:ncol(genotype_matrix_centered)])
    genotype_matrix_centered[is.na(genotype_matrix_centered)] = 0
  } else {
    genotype_matrix_centered <- t(genotype_matrix[,4:ncol(dosage_matrix)])
    marker_means                          = colMeans(genotype_matrix_centered, na.rm = TRUE)
    col_indices                           = which(is.na(genotype_matrix_centered), arr.ind = TRUE)
    genotype_matrix_centered[col_indices] = marker_means[genotype_matrix_centered[, 2]]
  }

  GEBV <- genotype_matrix_centered %*% marker_effects[,2]
  GEBV <- cbind(GEBV, 1:nrow(GEBV))

  #sort based on selecting for min or max
  if(maximize){
    GEBV <- GEBV[order(GEBV[,1], decreasing=TRUE), ]
  } else{
    GEBV <- GEBV[order(GEBV[,1], decreasing=FALSE), ]
  }

  TS_parent_names <- rownames(GEBV)[1:num_TS_parents]

  GA_geno <- cbind(marker_effects$SNP, phased_matrix[, GA_parent_names, drop=FALSE])
  colnames(GA_geno)[1] <- "name"
  TS_geno <- cbind(marker_effects$SNP, phased_matrix[, TS_parent_names, drop=FALSE])
  colnames(TS_geno)[1] <- "name"
  # ---- map ----
  if(is.null(genetic_map_position)){
    map  = map %>% dplyr::group_by(Chromosome) %>% dplyr::mutate(
      cM = max_cM_chr * (Position - min(Position)) / (max(Position) - min(Position))
    ) %>%
      dplyr::ungroup()
    map = map[, c("SNP", "Chromosome", "cM")]

    #if provided, just use the provided values
  } else{
    if(length(genetic_map_position) != nrow(map)  || anyNA(genetic_map_position) ) stop("The genetic_map_position vector must be the same length as the map file, in the same order, and not contain any NA values!")

    map$cM = genetic_map_position
    map    = map[, c("SNP", "Chromosome", "cM")]
  }

  #header names required for GenomicSimulation
  map                         <- map[ , c("Chromosome", "SNP", "cM")]
  colnames(map)               <- c("chr", "marker", "pos")

  colnames(marker_effects) <- c("marker","eff")
  marker_effects$allele <- "A"
  marker_effects <- marker_effects[,c(1,3,2)]

  temp_string <- paste0(sample(c(letters,LETTERS,0:9),12,replace=TRUE), collapse="")

  temp_files <- write_sim_files(
    GA_geno, TS_geno, map, marker_effects, temp_string
  )

  summary_df <- run_basic_simulation(
    temp_files,
    num_gen,
    num_sim_reps,
    num_cross_per_gen,
    ncol(GA_geno),
    ncol(TS_geno),
    maximize = maximize
  )

  delete_sim_files(temp_files)

  gc()

  sim_plot <- generate_sim_traject_plot(summary_df, colors[1:2])

  if(PCA){

    if(mean_adjust){
      genotype_matrix_centered <- center_genotypes(dosage_matrix)
      genotype_matrix_centered <- t(genotype_matrix_centered[,4:ncol(genotype_matrix_centered)])
      genotype_matrix_centered[is.na(genotype_matrix_centered)] = 0
    } else {
      genotype_matrix_centered <- t(genotype_matrix[,4:ncol(dosage_matrix)])
      marker_means                          = colMeans(genotype_matrix_centered, na.rm = TRUE)
      col_indices                           = which(is.na(genotype_matrix_centered), arr.ind = TRUE)
      genotype_matrix_centered[col_indices] = marker_means[genotype_matrix_centered[, 2]]
    }

    PCA_list <- generate_selection_pca_plot(
      genotype_matrix_centered,
      match(GA_parent_names, rownames(genotype_matrix_centered)),
      match(TS_parent_names, rownames(genotype_matrix_centered)),
      colors,
      alpha
    )
  }



  if(PCA){
    return(list(
      Simulation_Plot = sim_plot,
      PCA_Plot = PCA_list$PCA_plot,
      Simulation_Summary = summary_df,
      PCA_df = PCA_list$PCA_df
    ))
  } else {
    return(list(
      Simulation_Plot = sim_plot,
      Simulation_Summary = summary_df
    ))
  }
}
