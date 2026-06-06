#############################################
#### Unified Genetic Algorithm Framework ####
#############################################

############################
# 0. Block Selection       #
############################
select_top_blocks = function(haploblock_obj, n = NULL, perc_total = NULL, perc_of_total_var = NULL){

  #check only one value was provided
  if(sum(c(!is.null(n), !is.null(perc_total), !is.null(perc_of_total_var))) != 1 | any(c(is.character(n), is.character(perc_total), is.character(perc_of_total_var)))){
    stop("Please ensure only one of the following is non-null and numerical: number of blocks, percentage of total blocks, or blocks comprising percentage of total block variance is specified.")
  }

  #if percentages were specified, check they're within range
  if(!is.null(perc_total) && (perc_total <= 0 || perc_total > 1)){
    stop("Please ensure percentages are specified between 0 and 1.")
  }

  if(!is.null(perc_of_total_var) && (perc_of_total_var < 0 || perc_of_total_var > 1)){
    stop("Please ensure percentages are specified between 0 and 1.")
  }

  #order the blocks by variance
  haploblocks_df = haploblock_obj$Haploblocks
  haploblocks_df = haploblocks_df[order(haploblocks_df$Block_Var, decreasing = TRUE), ]

  #select blocks based on which parameter was not null - a set number
  if(!is.null(n)){
    haploblocks_df = haploblocks_df[1:n, ]
    block_id = haploblocks_df$Block_ID
    haploblock_obj$Haploblocks_GA = haploblocks_df
    haploblock_obj$Haplotype_Effect_Matrix_GA = as.data.frame(t(haploblock_obj$Haplotype_Effect_Matrix[block_id, ]))

    #at least a certain percentage of blocks (ceiling)
  } else if(!is.null(perc_total)){
    haploblocks_df = haploblocks_df[1:ceiling(nrow(haploblocks_df) * perc_total), ]
    block_id = haploblocks_df$Block_ID
    haploblock_obj$Haploblocks_GA = haploblocks_df
    haploblock_obj$Haplotype_Effect_Matrix_GA = as.data.frame(t(haploblock_obj$Haplotype_Effect_Matrix[block_id, ]))

    #percentage of blocks explaining at least a certain percentage of variance (ceiling)
  } else if(!is.null(perc_of_total_var)){
    cumulative_var = cumsum(haploblocks_df$Block_Var) / sum(haploblocks_df$Block_Var)
    haploblocks_df = haploblocks_df[1:which(cumulative_var >= perc_of_total_var)[1], ]
    block_id = haploblocks_df$Block_ID
    haploblock_obj$Haploblocks_GA = haploblocks_df
    haploblock_obj$Haplotype_Effect_Matrix_GA = as.data.frame(t(haploblock_obj$Haplotype_Effect_Matrix[block_id, ]))
  } else{
    stop("Please ensure only one of the selection methods is specified as non-NULL and is an appropriate numeric value.")
  }

  return(haploblock_obj)
}

############################
# 1. Strategy definitions  #
############################

VALID_STRATEGIES_LOCAL <- c(
  "selfing",
  "no_selfing"
)

VALID_STRATEGIES_OHS <- c(
  "self_allow_duplicate_chromosomes",
  "self_no_duplicate_chromosomes",
  "no_self_unique_individuals"
)

validate_strategy <- function(strategy, type = c("localGEBV", "OHS")) {
  type <- match.arg(type)

  if (type == "localGEBV" && !(strategy %in% VALID_STRATEGIES_LOCAL)) {
    stop("Invalid strategy for localGEBV")
  }

  if (type == "OHS" && !(strategy %in% VALID_STRATEGIES_OHS)) {
    stop("Invalid strategy for OHS")
  }
}


#############################################################################
###### 2. Helper: extract active matrices and build haplotype metadata ######
############################################################################

get_effect_matrix <- function(obj) {
  if (!is.null(obj$Haplotype_Effect_Matrix_GA)) {
    obj$Haplotype_Effect_Matrix_GA
  } else {
    stop("Be sure to use `select_top_blocks()` and use the object output here!")
  }
}

get_block_matrix <- function(obj) {
  if (!is.null(obj$Haploblocks_GA)) {
    obj$Haploblocks_GA
  } else {
    stop("Be sure to use `select_top_blocks()` and use the object output here!")
  }
}


build_row_metadata <- function(effect_matrix){

  rn <- rownames(effect_matrix)

  # OHS format
  if(any(grepl("_[0-9]+$", rn))){

    individual <- sub("_[0-9]+$", "", rn)
    chromosome  <- sub("^.*_([0-9]+)$", "\\1", rn)

    data.frame(
      row = seq_len(nrow(effect_matrix)),
      individual = individual,
      chromosome = chromosome,
      stringsAsFactors = FALSE
    )

  } else {

    # localGEBV format
    data.frame(
      row = seq_len(nrow(effect_matrix)),
      individual = rn,
      chromosome = NA,
      stringsAsFactors = FALSE
    )

  }
}

#############################################################
#### 3. CORE GA ENGINE (individual-based, fully shared) #####
#############################################################

.genetic_algorithm_core <- function(
    fitness_fn,
    n_individuals,
    n_founders,
    popSize = 100,
    maxiter = 500,
    run = 50,
    pmutation = 0.1,
    pcrossover = 0.1,
    monitor = TRUE
){


  # population initializer
  # create the initial populations by random sampling
  custom_population <- function(object){
    t(replicate(object@popSize,
                sample(1:n_individuals, n_founders, replace = FALSE)))
  }

  # Custom mutation operator for subset-selection GA.
  #
  # INPUTS
  #   object:
  #     GA object passed automatically by GA::ga().
  #     Contains the current population in:
  #       object@population
  #
  #   parents:
  #     Integer vector of row indices identifying the candidate
  #     solutions selected by GA for mutation.
  #
  #     IMPORTANT:
  #     These are NOT the founder IDs themselves.
  #     They are row numbers of solutions within the current
  #     GA population.
  #
  #     Example:
  #       parents = c(4, 17, 32)
  #
  #     means mutate population members:
  #       object@population[4, ]
  #       object@population[17, ]
  #       object@population[32, ]
  #
  #
  # MUTATION LOGIC
  #
  #   Each candidate solution is a vector of founder IDs:
  #
  #       c(12, 44, 7, 98, 31)
  #
  #   representing a subset of selected individuals.
  #
  #   For each solution selected for mutation:
  #
  #     1. Randomly choose one founder position.
  #
  #     2. Identify all individuals not already present in
  #        the current subset.
  #
  #     3. Replace the chosen founder with one randomly
  #        sampled unused individual.
  #
  #   This preserves:
  #
  #     - fixed subset size (n_founders)
  #     - uniqueness of founder IDs
  #     - valid founder IDs between 1:n_individuals
  #
  #
  # OUTPUT
  #
  #   Returns a matrix with:
  #
  #     rows = mutated solutions
  #     cols = founder positions
  #
  #   Dimensions:
  #
  #     length(parents) x n_founders
  #
  #   Example:
  #
  #     [,1] [,2] [,3] [,4] [,5]
  #   [1,]  12   44   55   98   31
  #   [2,]   7   21   83   16   99
  #
  #   GA expects the returned matrix to contain the mutated
  #   versions of the population members specified by
  #   'parents'. These rows replace the corresponding
  #   individuals in the next generation.
  #
  #
  # NOTE
  #
  #   Mutation probability is controlled by GA::ga()
  #   through the pmutation argument.
  #
  #   This function only defines HOW a selected solution
  #   is mutated, not WHETHER mutation occurs.
  #
  custom_mutation <- function(object, parents){

    #matrix of popSize (number of pops) x # of parents (founders) in a pop
    out <- matrix(
      NA,
      nrow = length(parents), #this is providing the pops that need to be mutated!
      ncol = n_founders
    )


    #parents is a list of the popSize containing n_founders if I recall correctly
    for(i in seq_along(parents)){


      #extract theset of founders from the current pop to make a mutation
      sol <- object@population[parents[i], ]

      #choose a random individual in that population
      pos <- sample(seq_along(sol), 1)

      #find new individuals that could replace it
      available <- setdiff(
        seq_len(n_individuals),
        sol
      )

      #sample from those new individuals to replace a population member
      if(length(available) > 0){
        sol[pos] <- sample(available, 1)
      }

      #return the new set of parents to the current population
      out[i, ] <- sol
    }

    out
  }

  # crossover
  # Custom crossover operator for subset-selection GA.
  #
  # INPUTS
  #   object:
  #     GA object passed automatically by GA::ga().
  #
  #   parents:
  #     Integer vector of length 2 containing the row
  #     indices of the two candidate solutions selected
  #     for crossover.
  #
  #     Example:
  #
  #       parents = c(8, 21)
  #
  #     means crossover:
  #
  #       object@population[8, ]
  #       object@population[21, ]
  #
  #
  # REPRESENTATION
  #
  #   Each solution is a vector of founder IDs:
  #
  #       c(12, 44, 7, 98, 31)
  #
  #   representing a subset of selected individuals.
  #
  #   Founder order has no biological meaning and the
  #   solution should be interpreted as a set.
  #
  #
  # CROSSOVER LOGIC
  #
  #   1. Randomly sample approximately half of the founders
  #      from parent 1.
  #
  #   2. Combine those founders with all founders from
  #      parent 2.
  #
  #   3. Remove duplicates.
  #
  #   4. Repeat the process in the opposite direction to
  #      create a second child.
  #
  #   This allows each child to inherit founders from both
  #   parental subsets.
  #
  #
  # REPAIR STEP
  #
  #   After combining founders, the resulting child may
  #   contain too many or too few founders.
  #
  #   Too many founders:
  #     Randomly sample down to n_founders.
  #
  #   Too few founders:
  #     Randomly add individuals not already present
  #     until n_founders founders are present.
  #
  #   This guarantees:
  #
  #     - fixed subset size
  #     - unique founder IDs
  #     - valid founder IDs
  #
  #
  # OUTPUT
  #
  #   Returns a list containing:
  #
  #     children:
  #       Matrix of offspring solutions.
  #
  #       Dimensions:
  #         2 x n_founders
  #
  #     fitness:
  #       Fitness values of the offspring.
  #
  #   Example:
  #
  #     $children
  #          [,1] [,2] [,3] [,4] [,5]
  #     child1   5   12   44   61   88
  #     child2   3   17   44   72   91
  #
  #     $fitness
  #     [1] 152.3 149.7
  #
  #   This is the structure expected by GA::ga().
  #
  #
  # NOTE
  #
  #   This is a set-based crossover rather than a classical
  #   positional crossover.
  #
  #   Founder order is not meaningful, so the operator
  #   effectively mixes parental subsets and then repairs
  #   the resulting offspring to satisfy the subset-size
  #   constraint.
  #
  custom_crossover <- function(object, parents){

    #select the vectors of parents provided by the GA
    p1 <- object@population[parents[1], ]
    p2 <- object@population[parents[2], ]

    #select a random half from each population
    take1 <- sample(seq_along(p1), ceiling(length(p1)/2))
    take2 <- sample(seq_along(p2), ceiling(length(p2)/2))

    #create the unique sets of parents
    child1 <- unique(c(
      p1[take1],
      p2
    ))

    child2 <- unique(c(
      p2[take2],
      p1
    ))

    #if there are too few/too many parents, we need to bring in more (mutation)
    #or extract the needed number
    repair <- function(x){

      #unique individuals
      x <- unique(x)

      #sample up or down
      if(length(x) > n_founders){
        x <- sample(x, n_founders)
      }

      if(length(x) < n_founders){

        available <- setdiff(
          seq_len(n_individuals),
          x
        )

        x <- c(
          x,
          sample(
            available,
            n_founders - length(x)
          )
        )
      }

      #return the new population
      x
    }

    child1 <- repair(child1)
    child2 <- repair(child2)

    children <- rbind(child1, child2)

    #calculate fitness
    fitness_values = apply(children, 1, fitness_fn)

    #return the new pops and fitness values
    list(
      children = children,
      fitness = fitness_values
    )
  }

  # run GA
  GA::ga(
    type = "real-valued",
    fitness = fitness_fn,
    lower = rep(1, n_founders),
    upper = rep(n_individuals, n_founders),
    popSize = popSize,
    maxiter = maxiter,
    run = run,
    pcrossover = pcrossover,
    pmutation = pmutation,
    population = custom_population,
    mutation = custom_mutation,
    crossover = custom_crossover,
    monitor = monitor
  )
}

########################################################
#####          4. LOCAL GEBV FITNESS               #####
########################################################

fitness_localGEBV <- function(
    effect_matrix,
    strategy = c(
      "selfing",
      "no_selfing"
    )
){

  strategy <- match.arg(strategy)

  function(selected_ind){

    #pull localGEBV
    selected <- effect_matrix[selected_ind, , drop = FALSE]

    # individual pairings

    if(nrow(selected) < 2 && strategy == "no_selfing"){
      return(0)
    }

    #make the individual pairings
    combos <- combn(
      seq_len(nrow(selected)),
      2,
      simplify = TRUE
    )

    combos <- t(combos)

    #add self to self if valid
    if(strategy == "selfing"){

      self_pairs <- cbind(
        seq_len(nrow(selected)),
        seq_len(nrow(selected))
      )

      combos <- rbind(
        combos,
        self_pairs
      )
    }

    total_score <- 0

    #compute max progeny GEBV for each block
    #for each block, evaluate each pair's fitness by averaging localGEBV (averaging across the 4 chromosomes)
    for(block in seq_len(ncol(selected))){

      vals <- selected[, block]

      scores <- (vals[combos[,1]] + vals[combos[,2]]) / 2 # future spot to update for ploidy difference

      total_score <- total_score + max(scores)
    }

    total_score
  }
}

########################################################
##### 5. OHS FITNESS (full haplotype-aware engine) #####
########################################################

fitness_OHS <- function(
    effect_matrix,
    row_metadata,
    strategy = c(
      "self_allow_duplicate_chromosomes",
      "self_no_duplicate_chromosomes",
      "no_self_unique_individuals"
    )
){

  strategy <- match.arg(strategy)

  individual_map <- split(
    row_metadata$row,
    row_metadata$individual
  )

  function(selected_ind){

    # expand selected individuals into chromosomes

    selected_names <- names(individual_map)[selected_ind]

    selected_rows <- unlist(
      individual_map[selected_names]
    )

    meta <- row_metadata[selected_rows, ]

    # chromosome pairings

    if(length(selected_rows) < 2){
      return(0)
    }

    # unique chromosome pairs
    combos <- combn(
      seq_along(selected_rows),
      2,
      simplify = TRUE
    )

    combos <- t(combos)

    # optionally allow selfing

    if(strategy == "self_allow_duplicate_chromosomes"){

      self_pairs <- cbind(
        seq_along(selected_rows),
        seq_along(selected_rows)
      )

      combos <- rbind(
        combos,
        self_pairs
      )
    }

    # strategy filtering - depending on the strategy chosen

    r1 <- combos[,1]
    r2 <- combos[,2]

    same_row <- r1 == r2

    same_individual <-
      meta$individual[r1] ==
      meta$individual[r2]

    # strategy rules

    keep <- rep(TRUE, nrow(combos))

    if(strategy == "self_no_duplicate_chromosomes"){
      keep <- !same_row
    }

    if(strategy == "no_self_unique_individuals"){
      keep <- !same_individual
    }

    combos <- combos[keep, , drop = FALSE]

    # compute optimal score

    total_score <- 0

    #edge case for
    if(nrow(combos) == 0){
      return(total_score)
    }


    #compute pairwise chromosome localGEBV (addition of haplotype GEBV)
    for(block in seq_len(ncol(effect_matrix))){



      #pull out chromosomes of selected parents
      vals <- effect_matrix[selected_rows, block]



      #for each valid chromosome combo (based on filtering above), compute haplo sum
      #currently limited to diploid! This is a point of future expansion.
      scores <- vals[combos[,1]] + vals[combos[,2]] #eventually need to alter for polyploid

      total_score <- total_score + max(scores)
    }

    total_score
  }
}

########################################################
#####       6. localGEBV wrapper                   #####
########################################################

local_gebv_parent_selection <- function(
    haploblock_obj,
    strategy = c("selfing", "no_selfing"),
    n_founders = 20,
    popSize = 100,
    maxiter = 500,
    run = 50,
    pmutation = 0.1,
    pcrossover = 0.1,
    monitor = TRUE
){

  strategy <- match.arg(strategy)

  validate_strategy(
    strategy,
    type = "localGEBV"
  )

  effect_matrix <- as.matrix(
    get_effect_matrix(haploblock_obj)
  )

  fitness_fn <- fitness_localGEBV(
    effect_matrix = effect_matrix,
    strategy = strategy
  )

  ga_result <- .genetic_algorithm_core(
    fitness_fn = fitness_fn,
    n_individuals = nrow(effect_matrix),
    n_founders = n_founders,
    popSize = popSize,
    maxiter = maxiter,
    run = run,
    pmutation = pmutation,
    pcrossover = pcrossover,
    monitor = monitor
  )

  solution <- ga_result@solution

  if(is.matrix(solution)){
    solution <- solution[1, ]
  }

  solution <- sort(as.integer(solution))

  unique_individuals <- row.names(effect_matrix)

  One_Solution = data.frame(
    Index = solution,
    Individual = unique_individuals[solution]
  )


    return_obj = list(GA = ga_result, One_Solution = One_Solution)
    return(return_obj)
}

########################################################
#####             7. OHS wrapper                   #####
########################################################

ohs_parent_selection <- function(
    haploblock_obj,
    strategy = c(
      "self_allow_duplicate_chromosomes",
      "self_no_duplicate_chromosomes",
      "no_self_unique_individuals"
    ),
    n_founders = 20,
    popSize = 100,
    maxiter = 500,
    run = 50,
    pmutation = 0.1,
    pcrossover = 0.1,
    monitor = TRUE
){

  strategy <- match.arg(strategy)

  validate_strategy(
    strategy,
    type = "OHS"
  )

  effect_matrix <- as.matrix(
    get_effect_matrix(haploblock_obj)
  )

  row_metadata <- build_row_metadata(
    effect_matrix
  )

  n_individuals <- length(
    unique(row_metadata$individual)
  )

  fitness_fn <- fitness_OHS(
    effect_matrix = effect_matrix,
    row_metadata = row_metadata,
    strategy = strategy
  )

  ga_result <- .genetic_algorithm_core(
    fitness_fn = fitness_fn,
    n_individuals = n_individuals,
    n_founders = n_founders,
    popSize = popSize,
    maxiter = maxiter,
    run = run,
    pmutation = pmutation,
    pcrossover = pcrossover,
    monitor = monitor
  )

  solution <- ga_result@solution

  if(is.matrix(solution)){
    solution <- solution[1, ]
  }

  solution <- sort(as.integer(solution))

  unique_individuals <- unique(
    row_metadata$individual
  )

  One_Solution = data.frame(
    Index = solution,
    Individual = unique_individuals[solution]
  )

  return_obj = list(GA = ga_result, One_Solution = One_Solution)
  return(return_obj)
}
