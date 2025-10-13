######Pairwise LD#########

#loading map file from an R object
#the map file chromosomes MUST be numeric!
#First column should be SNP ID, second column should be chromosome, and third
#column should be the position
load(file = "Example_Files/gapit_map.R")

#simulate unordered map file - original file is already ordered, so this is just an example
map2 = map[sample(1:nrow(map), nrow(map)), ]

#check that the map file is ordered - using simulated unordered map
map2 = order_map(map = map2)

#loading the genotype file from an R object
load(file = "Example_Files/gapit_genos.R")

#compute ld
ld_pairs = pairwise_ld(geno)

#####make haploblocks#####

#haploblocking - using a threshold of 0.4, a tolerance of 2, 
#resetting the tolerance between successful marker additions,
#forming blocks around the highest LD pairs, and using LD between adjacent markers
#(i.e., not average LD to the block)
haploblocks = def_blocks(ld = ld_pairs, map = map, method = "flanking",
                          threshold = 0.2, tolerance = 4, tol_reset = TRUE, start = "LD")
haploblocks = block_obj_to_df(haploblocks, map)


##### compute localGEBV using marker effects #####

#the marker pecov (prediction error covariance matrix) is only needed if you want to compute
#p-values of specific haplotype effects - very much theoretical at this point
load(file = "Example_Files/gapit_marker_pecov.R")
load(file = "Example_Files/gapit_marker_effects.R")

#ideally you just set missing values to NA. This is fine for LD calculations,
#calculating marker effects, and computing localGEBV
#I'm just keeping it simple for this example and rounding
geno[,4:ncol(geno)] = round(geno[,4:ncol(geno)], 0)

haploblock_obj = compute_local_GEBV(geno = geno, marker_effects = marker_effects, haploblocks_df = haploblocks, marker_pecov = marker_pecov)

###### do some visualizations ######

#must provide the column of the marker effects as well as the column indicating the chromosome and position on the chromosome
#if the map and marker effects plot are in the same order, you can do it as it is done here, else merge the two data frames together
marker_plot = marker_effects_plot(marker_effects = marker_effects$Effect, chr = map$Chromosome, pos = map$Position)

#all of these plots utilize the haploblock_obj created
haplo_eff_plot = unique_haplo_effects_plot(haplo_obj = haploblock_obj)

#mean_line = TRUE adds a line to the block variance (scaled)
funnel_plot = block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE)

haploblock_plot = plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks)


###### select top 15 haploblocks (arbitrary) and perform the GA ######

haploblock_effects = haploblock_obj$Haploblocks

#order the blocks by block variance
haploblock_effects = haploblock_effects[order(haploblock_effects$Block_Var, decreasing = TRUE), ]

#select the top 15
haploblock_top_blocks = haploblock_effects[1:15, ]

#pull out the localGEBV corresponding to the top 15
localGEBV = haploblock_obj$Haplotype_Effect_Matrix
localGEBV = localGEBV[row.names(localGEBV) %in% haploblock_effects$Block_ID, ]

#transpose it and turn it into a matrix for the GA
localGEBV = as.data.frame(t(as.matrix(localGEBV)))


#run the GA
#n_founders is the number of parents you want to select
#popSize is the number of groups of 20 you want - more popSize takes longer per iteration, but leads to faster convergence
#maxiter is the maximum number of iterations to reach an optimum solution
#run is the maximum number of iterations without improvement before the GA stops
#selfing determines whether a parent's localGEBV can be passed down to offspring twice (inbreeding with selfing) in the fitness function
#pmutation is the probability a random individual from each population (popSize) of size n_founders will be subsituted with a random individual
#in the entire population (not including those already in the set of n_founders). Helps reach convergence more quickly, but too large a value may hinder finding the optimum
#pcrossover is the probability of two populations (popSize) of n_founders exchanging half of their individuals - if the unique set of individuals after swapping is less than n_founders per population
#then a random subset is sampled from the overall population like with pmutation (i.e., not in the unique set of individuals)
#pelite when pulling random subsets from the base population after pcrossover, pelite will only consider the pelite proportion individuals with the greatest GEBV (fitness)
#set to 1 to consider all individuals
GA_output = genetic_algorithm(localGEBV = localGEBV, n_founders = 20, popSize = 10, maxiter = 300, run = 150, selfing = FALSE, pmutation = 0.2, pcrossover = 0.8, pelite = 0.5)

#one unique set of parents - other solutions may exist, see the GA_output$GA@solution object for other potential sets of parents with the same fitness
GA_output$One_Solution