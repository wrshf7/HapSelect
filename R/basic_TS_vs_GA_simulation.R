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


GA_vs_TS_simulation = function(GA_output, geno, marker_effects, map, genetic_map_position = NULL, num_gen = 50, num_sim_reps = 30,
                               num_cross_per_gen = 1000, num_TS_parents = NULL, mean_adjust = TRUE, max_cM_chr = 100,
                               colors = c("#A01FF0", "#A7A8AA")){

  #check compatability
  check_geno_marker_compatibility(geno = geno, marker_effects = marker_effects, map = map)

  #extract the parent names as well as the parent indices in the genotype file (after removing map info)
  GA_parents_indices = GA_output$One_Solution$Indices
  num_GA_parents = length(GA_parents_indices)

  #default TS parent nubmer to GA parent number
  if(missing(num_TS_parents) || is.na(num_TS_parents) || !is.numeric(num_TS_parents)){
    num_TS_parents = num_GA_parents
  }

  #check the number of crosses per generation
  if(length(num_cross_per_gen) >1  || anyNA(as.numeric(num_cross_per_gen))) stop("The num_cross_per_gen value must either be NULL or a single numeric and cannot be NA!")

  #check the number of crosses per generation
  if(length(num_gen) >1  || anyNA(as.numeric(num_gen))) stop("The num_gen value must either be NULL or a single numeric and cannot be NA!")



  #extract the genotype matrix and center it using centering function from localGEBV_calculation.R
  genotype_matrix = t(geno[,4:ncol(geno)])

  #mean adjust and fill in missing values with the mean to compute GEBV
  if(mean_adjust){
    genotype_matrix_centered   = center_genotypes(geno)
    genotype_matrix_centered   = t(genotype_matrix_centered[,4:ncol(genotype_matrix_centered)])
    genotype_matrix_centered[is.na(genotype_matrix_centered)] = 0

  } else{
    genotype_matrix_centered              = genotype_matrix
    marker_means                          = colMeans(genotype_matrix_centered, na.rm = TRUE)
    col_indices                           = which(is.na(genotype_matrix_centered), arr.ind = TRUE)
    genotype_matrix_centered[col_indices] = marker_means[col_indices[, 2]]
  }


  marker_effects = marker_effects[match(marker_effects$SNP, colnames(genotype_matrix_centered)), ]

  #compute GEBV, order them, and extract indices of num_TS_parents
  GEBV       = genotype_matrix_centered %*% marker_effects$Effect
  GEBV       = cbind(GEBV, 1:nrow(GEBV))
  GEBV       = GEBV[order(GEBV[,1], decreasing = TRUE), ]
  TS_indices = GEBV[1:num_TS_parents, 2]

  #unfortunately GenomicSimulation doesn't handle missing values well right now, so we set it to heterozygous (must take hard dosage calls)
  #for true haplotypes there may be another option from what I recall in GenomicSimulation
  genotype_matrix[is.na(genotype_matrix)] = 1

  #extract genotypes of TS parents and GA parents
  genotype_matrix_TS = genotype_matrix[TS_indices, ]
  genotype_matrix_GA = genotype_matrix[GA_parents_indices, ]

  GA_parent_names    = row.names(genotype_matrix_GA)
  TS_parent_names    = row.names(genotype_matrix_TS)

  #sets the genetic distance between first and last marker to be 100 cM, spaced proportionally to physical distance
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
  map                         = map[ , c("Chromosome", "SNP", "cM")]
  colnames(map)               = c("chr", "marker", "pos")
  colnames(marker_effects)    = c("marker", "eff")
  marker_effects$allele       = "A"
  marker_effects              = marker_effects[,c(1,3,2)]

  #I put a random string to allow simultaneous runs... but the package doesn't support that last I checked (package made own temp files with same name)
  #Simultaneous runs would need to be in different directories

  temp_string = paste0(
    sample(c(letters, LETTERS, 0:9), 12, replace = TRUE),
    collapse = ""
  )



  #write out the temp files needed for GenomicSimulation
  temp_files = write_sim_files(GA_geno = genotype_matrix_GA, TS_geno = genotype_matrix_TS, map = map, marker_effects = marker_effects, temp_string = temp_string)


  #run simulation
  summary_df = run_basic_simulation(temp_files = temp_files, num_gen = num_gen, num_sim_reps = num_sim_reps, num_cross_per_gen = num_cross_per_gen,
                                    num_GA_parents = num_GA_parents, num_TS_parents = num_TS_parents)
  #delete temp files
  delete_sim_files(temp_files)

  #make plot
  sim_plot = generate_sim_traject_plot(summary_df = summary_df, colors = colors)

  return(sim_plot)
}

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
run_basic_simulation = function(temp_files, num_gen, num_sim_reps, num_cross_per_gen, num_GA_parents, num_TS_parents){

  #tracking
  mean_BV_GA = list()
  mean_BV_TS = list()


  #run it for the GA parents
  #initiate populations
  init = genomicSimulation::load.data(allele.file = temp_files[[2]],
                   map.file = temp_files[[1]],
                   effect.file = temp_files[[4]])
  g0 = init$groupNum
  current_pop = genomicSimulation::make.random.crosses(g0, n.crosses = num_cross_per_gen, give.names = FALSE, name.prefix = "F1.")

  #simulate across reps
  for(rep in seq_len(num_sim_reps)){
    mean_BV_GA[[rep]] = sim_gens(num_gen = num_gen, init = init,
                                 num_parents = num_GA_parents, num_cross_per_gen = num_cross_per_gen,
                                 current_pop = current_pop)
  }

  #run it for the TS parents
  #initiate populations
  init = genomicSimulation::load.data(allele.file = temp_files[[3]],
                                      map.file = temp_files[[1]],
                                      effect.file = temp_files[[4]])
  g0 = init$groupNum
  current_pop = genomicSimulation::make.random.crosses(g0, n.crosses = num_cross_per_gen, give.names = FALSE, name.prefix = "F1.")

  for(rep in seq_len(num_sim_reps)){
    mean_BV_TS[[rep]] = sim_gens(num_gen = num_gen, init = init,
                                 num_parents = num_TS_parents, num_cross_per_gen = num_cross_per_gen,
                                 current_pop = current_pop)
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
sim_gens = function(num_gen, init, num_parents, num_cross_per_gen, current_pop){
  mean_BV = c()
  for(gen in seq_len(num_gen)){
    # Record mean BV of current population
    mean_BV[gen] = mean(genomicSimulation::see.GEBVs(current_pop))

    # Select top parents by GEBV
    selected = genomicSimulation::break.group.by.GEBV(
      current_pop,
      number = num_parents
    )

    # Generate next generation (1000 offspring)
    current_pop = genomicSimulation::make.random.crosses(
      group       = selected,
      n.crosses   = num_cross_per_gen,
      give.names  = TRUE,
      name.prefix = paste0("F", gen + 1, ".")
    )
  }

  return(mean_BV)
}

# plotting function
generate_sim_traject_plot = function(summary_df, colors){
  traject_plot = ggplot2::ggplot(summary_df, aes(x = gen, y = mean, color = method, fill = method)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se, group = method), alpha = 0.2, show.legend = FALSE) +
    labs(x = "Generation", y = "Mean Breeding Value") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors, guide = "none")

  return(traject_plot)
}

