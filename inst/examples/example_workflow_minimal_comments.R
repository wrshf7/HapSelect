##################################################
####  Example Workflow with Minimal Comments  ####
##################################################

library(HapSelect)

#############################################################
####  Pairwise LD with Native Function (not recommended) ####
#############################################################

data("map", package = "HapSelect")

#simulate unordered map file - original file is already ordered, so this is just an example
map2 = map[sample(1:nrow(map), nrow(map)), ]

#check that the map file is ordered - using simulated unordered map
map2 = order_map(map = map2)

#loading the genotype file from an R object
#first 3 columns of the genotype file should be the map file and the rest of the columns represent accessions/individuals with cells
#being the genotypes at each locus for each individual
data("geno", package = "HapSelect")


#compute ld (painfully slow)
#ld_pairs = pairwise_ld(geno, parallelize = FALSE)

#This is very slow and is just for demonstration purposes only!



#load the example file to see the structure if you do not want to run the function
data("pairwise_ld", package = "HapSelect")



####################################################
####  Pairwise LD with PLINK Call (recommended) ####
####################################################

data("map", package = "HapSelect")
data("geno", package = "HapSelect")


#compute ld with PLINK

#once LD File Script is updated run this
ld_pairs = plink_pairwise_ld(prefix = "example_plink", ld_window = 999999, ld_window_kb = 1000000,
                             ld_window_r2 = 0, extra_args = character())

###########################################
##### Load LD File (Skip Computation) #####
###########################################

#load the example file to see the structure if you do not want to run the function
data("pairwise_ld", package = "HapSelect")

#########################################
########  Compute Haploblocks  ##########
#########################################

#haploblocking - using a threshold of 0.4, a tolerance of 2 as an example,
#resetting the tolerance between successful marker additions,
#forming blocks around the highest LD pairs, and using LD between adjacent markers (flanking)
#(i.e., not average LD to the block)

haploblocks = def_blocks(ld = ld_pairs, map = map, method = "flanking",
                         threshold = 0.2, tolerance = 4, tol_reset = TRUE,
                         start = "LD", parallel = FALSE)

#turn the block object into a data frame
haploblocks = block_obj_to_df(haploblocks, map)

#summarize block information
block_summary(block_df = haploblocks)

#######################################
#### Compute localGEBV  #####
#######################################


#the first column in the marker effects file should be the SNP ID (same as the map file) and the second
#column should be the marker effect estimated from a model.
data("marker_effects", package = "HapSelect")

#compute localGEBV
haploblock_obj = compute_local_GEBV(geno = geno, marker_effects = marker_effects, haploblocks_df = haploblocks,
                                    set_missing_NA = TRUE, mean_adjust = TRUE)
##################################
####  Visualizations  #####
##################################

#must provide the column of the marker effects as well as the column indicating the chromosome and position on the chromosome
#if the map and marker effects plot are in the same order, you can do it as it is done here, else merge the two data frames together
marker_plot = marker_effects_plot(marker_effects = marker_effects$Effect, chr = map$Chromosome, pos = map$Position)
marker_plot

#all of these plots utilize the haploblock_obj created
haplo_eff_plot = unique_haplo_effects_plot(haplo_obj = haploblock_obj)
haplo_eff_plot

#mean_line = TRUE adds a line to the block variance (scaled)
funnel_plot = block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE)
funnel_plot

haploblock_plot = plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks)
haploblock_plot

#marker density plot
marker_density_plot = plot_marker_density(map_df = map, bin_size = 500000)
marker_density_plot

#LD decay plot
#plot the LD decay utilizing the map and LD objects
ld_decay_plot = plot_ld_decay(map = map, ld = ld_pairs, max_kb = 500, span = 0.3, k = 10, method = "gam_cr")
ld_decay_plot


##############################################
#### Parent Selection with the GA ####
##############################################

#There are three methods we have implemented to select haploblocks:
#1: select the top n blocks based on Block_Var
#2: select the top x% of blocks (round up)
#3: select the top blocks explaining x% of the total sum of block variance (round up)


#1: select top 15 haploblocks (arbitrary)

haploblock_obj = select_top_blocks(haploblock_obj = haploblock_obj, n = 15)

#Objects are added to the haploblock_obj for the GA
nrow(haploblock_obj$Haploblocks_GA)

#2 select top 50% of haploblocks
haploblock_obj = select_top_blocks(haploblock_obj = haploblock_obj, perc_total = 0.5)
nrow(haploblock_obj$Haploblocks_GA)

#3 select the top blocks explaining at least 90% of the total block variance
haploblock_obj = select_top_blocks(haploblock_obj = haploblock_obj, perc_of_total_var = 0.9)
nrow(haploblock_obj$Haploblocks_GA)


#run the GA
#n_founders  is the number of parents you want to select
#popSize     is the number of groups of 20 you want - more popSize takes longer per iteration, but leads to faster convergence
#maxiter     is the maximum number of iterations to reach an optimum solution
#run         is the maximum number of iterations without improvement before the GA stops
#selfing     determines whether a parent's localGEBV can be passed down to offspring twice (inbreeding with selfing) in the fitness function
#pmutation   is the probability a random individual from each population (popSize) of size n_founders will be substituted with a random individual
#            in the entire population (not including those already in the set of n_founders). Helps reach convergence more quickly, but too large a value may hinder finding the optimum
#pcrossover  is the probability of two populations (popSize) of n_founders exchanging half of their individuals - if the unique set of individuals after swapping is less than n_founders per population
#            then a random subset is sampled from the overall population like with pmutation (i.e., not in the unique set of individuals)
#pelite      when pulling random subsets from the base population after pcrossover, pelite will only consider the pelite proportion individuals with the greatest GEBV (fitness)
#            set to 1 to consider all individuals
GA_output = genetic_algorithm(localGEBV = haploblock_obj$Haplotype_Effect_Matrix_GA, n_founders = 20, popSize = 10, maxiter = 300,
                              run = 150, selfing = FALSE, pmutation = 0.2, pcrossover = 0.8, pelite = 0.5)

#one unique set of parents - other solutions may exist, see the GA_output$GA@solution object for other potential sets of parents with the same fitness
GA_output$One_Solution
