library(tidyverse)
library(GA)

load("Example_Files/gapit_haploblock_obj.R")

haploblock_effects = gapit_haploblock_obj$Haploblocks
haploblock_effects = haploblock_effects[order(haploblock_effects$Block_Effect_Var, decreasing = TRUE), ]

haploblock_top_blocks = haploblock_effects[1:15, ]

localGEBV = gapit_haploblock_obj$Haplotype_Effect_Matrix
localGEBV = localGEBV[row.names(localGEBV) %in% haploblock_effects$Block_ID, ]

localGEBV = as.data.frame(t(as.matrix(localGEBV)))



GA_output = genetic_algorithm(localGEBV = localGEBV, n_founders = 20, popSize = 10, maxiter = 300, run = 150, selfing = FALSE, pmutation = 0.2, pcrossover = 0.8, pelite = 0.5)

####main function####
genetic_algorithm = function(localGEBV, n_founders = 20, popSize = 100, maxiter = 500, run = 50, selfing = FALSE, pmutation = 0.1, pcrossover = 0.8, pelite = 0.5){
  # Define number of individuals and number of haplotype blocks based on input matrix
  n_individuals = nrow(localGEBV)   # <-- total number of individuals to choose from
  n_blocks = ncol(localGEBV)        # <-- number of genomic regions (haploblocks)
  
  #### fitness function ####
  # Calculates total score by summing the best expected offspring GEBV per block
  ohs_fitness_custom = function(localGEBV, selected_ind, same_ind_ok){
    # Subset localGEBV matrix to only selected individuals
    sel_localGEBV = localGEBV[selected_ind, , drop = FALSE]
    
    # Create all pairwise combinations of selected individuals (no selfing)
    combos = combn(1:nrow(sel_localGEBV), 2)

    if(selfing == TRUE){ # If selfing is true allow for selfing pairs
      combos2 = rbind(1:nrow(sel_localGEBV), 1:nrow(sel_localGEBV)) # Introduce selfing pairs
      combos = cbind(combos, combos2)
    }
    
    # Initialize running total for combined fitness
    total_score = 0
    
    # For each haplotype block
    for(block in 1:n_blocks){
      # Extract local GEBV for that block
      block_values = sel_localGEBV[ , block]
      
      # Compute average offspring GEBV for each pair
      avg_scores = colMeans(matrix(block_values[combos], nrow = 2))
      
      # Take the highest-scoring pair for the block
      if (all(is.na(avg_scores))) {
        block_score = -1e6   # large negative penalty
      } else {
        block_score = max(avg_scores, na.rm = TRUE)
      }
      total_score = total_score + block_score
    }
    
    return(total_score)
  }
  
  #### fitness wrapper function ####
  # Wrapper that limits the solution to n_founders and calls the fitness function - GA function doesn't allow for multiple arguments
  fitness_wrapper = function(selected_ind){
    selected_ind = selected_ind[1:n_founders]  # Ensure the input is the right length
    total_score = ohs_fitness_custom(localGEBV = localGEBV, selected_ind = selected_ind, same_ind_ok = FALSE)
    return(total_score)
  }
  
  #### define founder pops function ####
  # Generates a random initial population of solutions (each is a founder set)
  # popSize number of founder solutions with n_founder number of founders per solution set (rows are populations, columns are the founders)
  custom_population = function(object){
    pops = t(replicate(object@popSize, sample(1:n_individuals, n_founders, replace = FALSE)))  # Each row is a candidate solution
    return(pops)
  }
  
  #### custom mutation function #####
  # Randomly mutate one element of a founder set by replacing it with an unused individual
  # Selects the replacement from the total base pool minus those already in the solution set (population)
  custom_mutation = function(object, parents){
    mutated = matrix(NA, nrow = length(parents), ncol = n_founders)
    for (i in seq_along(parents)) {
      ind = object@population[parents[i], ]                         # Get parent solution
      pos = sample(n_founders, 1)                                   # Random position to mutate
      new_val = sample(setdiff(1:n_individuals, ind), 1)            # Sample a new individual not already selected
      ind[pos] = new_val                                            # Replace
      mutated[i, ] = ind                                            # Save mutated solution
    }
    return(mutated)
  }
  
  #### custom crossover function ####
  # Combines two parent solutions into two new children, preserving some structure and enforcing uniqueness
  # If not enough individuals, pulls randomly from the base population sans individuals already in the unique set
  # If too many unique individuals, randomly samples them
  # Future direction: pull from the top percentage of population solutions?
  # custom_crossover = function(object, parents){
  #   parent1 = object@population[parents[1], ]
  #   parent2 = object@population[parents[2], ]
  #   
  #   half = floor(n_founders / 2)
  #   
  #   # Keep first half from one parent and fill in from the other
  #   base1 = parent1[1:half]
  #   base2 = parent2[1:half]
  #   
  #   add1 = setdiff(parent2, base1)
  #   add2 = setdiff(parent1, base2)
  #   
  #   child1 = unique(c(base1, add1))
  #   child2 = unique(c(base2, add2))
  #   
  #   # Ensure exactly n_founders in each child - pull randomly from base population to fill in
  #   # Sample from existing if too many
  #   if (length(child1) < n_founders) {
  #     filler1 = setdiff(1:n_individuals, child1)
  #     child1 = c(child1, sample(filler1, n_founders - length(child1)))
  #   } else if (length(child1) > n_founders) {
  #     child1 = sample(child1, n_founders)
  #   }
  #   
  #   if (length(child2) < n_founders) {
  #     filler2 = setdiff(1:n_individuals, child2)
  #     child2 = c(child2, sample(filler2, n_founders - length(child2)))
  #   } else if (length(child2) > n_founders) {
  #     child2 = sample(child2, n_founders)
  #   }
  #   
  #   # Combine into matrix and compute fitness for both children
  #   children = rbind(as.integer(child1), as.integer(child2))
  #   fitness_values = apply(children, 1, fitness_wrapper)
  #   
  #   return(list(children = children, fitness = fitness_values))
  # }
  
  custom_crossover = function(object, parents){
    # Extract the two parent solutions (each is a vector of founder IDs)
    parent1 = object@population[parents[1], ]
    parent2 = object@population[parents[2], ]
    
    # Take the first half of each parent to start building the children
    half = floor(n_founders / 2)
    base1 = parent1[1:half]
    base2 = parent2[1:half]
    
    # Add founders from the other parent that are not already in the base
    add1 = setdiff(parent2, base1)
    add2 = setdiff(parent1, base2)
    
    # Merge bases and additions to get initial child solutions
    child1 = unique(c(base1, add1))
    child2 = unique(c(base2, add2))
    
    ## ============================================================
    ## NEW: Limit filler pool to top X% of population by fitness
    ## ============================================================
    
    # Get the fitness values for the whole current population
    pop_fitness = object@fitness
    
    # Calculate number of elite individuals to keep (at least 1)
    elite_size = max(1, floor(length(pop_fitness) * pelite))
    
    # Get indices of the top-performing individuals (highest fitness first)
    elite_indices = order(pop_fitness, decreasing = TRUE)[1:elite_size]
    
    # Collect all founders used by these elite individuals (unique set)
    elite_founders = unique(as.vector(object@population[elite_indices, ]))
    
    ## ------------------------------------------------------------
    ## Fill missing founders in each child
    ## ------------------------------------------------------------
    
    # Fill missing positions in Child 1 using elite founders
    if (length(child1) < n_founders) {
      filler1 = setdiff(elite_founders, child1)  # avoid duplicates
      # If not enough elite founders available, fill from full population
      if (length(filler1) < (n_founders - length(child1))) {
        filler1 = c(filler1, setdiff(1:n_individuals, child1))
      }
      # Randomly sample from filler pool to complete the child
      child1 = c(child1, sample(filler1, n_founders - length(child1)))
    } 
    # If child has too many founders, trim down to exactly n_founders
    else if (length(child1) > n_founders) {
      child1 = sample(child1, n_founders)
    }
    
    # Fill missing positions in Child 2 using elite founders
    if (length(child2) < n_founders) {
      filler2 = setdiff(elite_founders, child2)
      if (length(filler2) < (n_founders - length(child2))) {
        filler2 = c(filler2, setdiff(1:n_individuals, child2))
      }
      child2 = c(child2, sample(filler2, n_founders - length(child2)))
    } 
    else if (length(child2) > n_founders) {
      child2 = sample(child2, n_founders)
    }
    
    ## ------------------------------------------------------------
    ## Evaluate the children
    ## ------------------------------------------------------------
    
    # Combine children into a matrix (rows = individuals, cols = founders)
    children = rbind(as.integer(child1), as.integer(child2))
    
    # Compute fitness for each child using your wrapper function
    fitness_values = apply(children, 1, fitness_wrapper)
    
    # Return list containing children and their fitness values
    return(list(children = children, fitness = fitness_values))
  }
  
  
  #### Run GA ####
  ga_result = ga(
    type = "real-valued",                            # Each solution is a numeric vector (indices)
    fitness = fitness_wrapper,                       # Function to evaluate solution fitness
    population = custom_population,                  # Custom initial population generator
    mutation = custom_mutation,                      # Custom mutation operator
    crossover = custom_crossover,                    # Custom crossover operator
    lower = rep(1, n_founders),                      # Lower bound (index 1)
    upper = rep(n_individuals, n_founders),          # Upper bound (index = n_individuals)
    popSize = popSize,                               # Number of founder sets/solutions
    maxiter = maxiter,                               # Max number of iterations to run the GA
    run = run,                                       # Stop if best solution doesn't improve for this many generations
    monitor = TRUE                                   # Show progress
  )
  
  # Get one best solution (first row) and return both its indices and names
  one_solution = sort(ga_result@solution[1,])
  one_solution = data.frame(
    Indices = one_solution,
    Names = rownames(localGEBV)[one_solution]        # Look up individual names by row index
  )
  
  return(list(GA = ga_result, One_Solution = one_solution))  # Return GA object and best solution
}
