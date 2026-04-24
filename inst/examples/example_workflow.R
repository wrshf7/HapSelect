#####################################
#######  Example Workflow  ##########
#####################################

library(HapSelect)

#############################################################
####  Pairwise LD with Native Function (not recommended) ####
#############################################################


#loading map file from an R object
#the map file chromosomes MUST be numeric!
#First column should be SNP ID, second column should be chromosome, and third
#column should be the position

data("map", package = "HapSelect")

#simulate unordered map file - original file is already ordered, so this is just an example
map2 = map[sample(1:nrow(map), nrow(map)), ]

#check that the map file is ordered - using simulated unordered map
#this function will also check the map file to ensure it is appropriate and set column names
map2 = order_map(map = map2)

#loading the genotype file from an R object
#with well-curated marker data, ideally you just set missing genotype values to NA. This would be fine for LD calculations,
#calculating marker effects, and computing localGEBV

#loading the genotype file from an R object
#Non-integers represent imputed values (usually to the mean or based on probability of each allele). It is recommended these be imputed to genotype calls
#or set to missing, otherwise haplotype configurations will be weird values. NOTE: This will not affect localGEBV values.
data("geno", package = "HapSelect")


#compute ld

#first 3 columns of the genotype file should be the map file and the rest of the columns represent accessions/individuals with cells
#being the genotypes at each locus for each individual
#column names of the map information (first 3 columns) will not matter in this case, but the order of map info columns
#should be the same as the map file columns (SNP ID, Chromosome, Position)

#ld_pairs = pairwise_ld(geno, parallelize = FALSE)

#This is very slow and is just for demonstration purposes only!



#load the example file to see the structure if you do not want to run the function
data("pairwise_ld", package = "HapSelect")



####################################################
####  Pairwise LD with PLINK Call (recommended) ####
####################################################


#loading map file from an R object
#the map file chromosomes MUST be numeric!
#First column should be SNP ID, second column should be chromosome, and third
#column should be the position

data("map", package = "HapSelect")

#simulate unordered map file - original file is already ordered, so this is just an example
map2 = map[sample(1:nrow(map), nrow(map)), ]

#check that the map file is ordered - using simulated unordered map
#this function will also check the map file to ensure it is appropriate and set column names
map2 = order_map(map = map2)

#loading the genotype file from an R object
#with well-curated marker data, ideally you just set missing genotype values to NA. This would be fine for LD calculations,
#calculating marker effects, and computing localGEBV

#Non-integers represent imputed values (usually to the mean or based on probability of each allele). These need to be imputed to genotype calls or set to missing.

#first 3 columns of the genotype file should be the map file and the rest of the columns represent accessions/individuals with cells
#being the genotypes at each locus for each individual
#column names of the map information (first 3 columns) will not matter in this case, but the order of map info columns
#should be the same as the map file columns (SNP ID, Chromosome, Position)
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


#note: other programs can be utilized to generate the LD file, but make sure the columns c("Chrom", "Locus1", "Locus2", "Name1", "Name2", "LD")
#exist in the data frame. We recommend using PLINK v1.9 for LD calculations.
#If you have a PLINK binary fileset (.bed/.bim/.fam), you can also use:
#ld_pairs = plink_pairwise_ld("path/to/plink_prefix")

#Locus1 and Locus2 are numeric indices and MUST reference the order of the markers (i.e., ordered by chromosome and position in chromosome)
#A good way to do this is to assign Locus1 as 1:nrow(map) to the ordered map file. You can then join the LD data frame and
#map data frame by Name1 and SNP and them Name2 and SNP to get the indices for Locus1 and Locus2.

#########################################
########  Compute Haploblocks  ##########
#########################################

#method = c("flanking", "average") indicates whether to compare the LD threshold between the last marker in a block and the new marker ("flanking")
#or to compare the average LD between all markers currently in the block and a new marker under consideration ("average")

#threshold = 0 to 1. This indicates the r^2 (LD) value cutoff after which a block terminates (or tolerance comes into play)

#tolerance = integer value (0, 1, 2, etc.). This indicates whether a gap is allowed. Sometimes markers can be poorly positioned, genotyping errors
#occur, or other anomalies (e.g., new/rare mutation of interest on a chip) happen and some markers may be in low LD while subsequent markers are in high LD
#the tolerance tells the program how many markers out to look for the LD threshold before terminating a block

#tol_reset = TRUE indicates that the tolerance counter is reset once a marker is successfully added (i.e., each "gap"
#is allowed to have the tolerance integer value of markers). If this is set to FALSE, then this is cumulative across "gaps" in a block.

#start = c("LD", "beginning"). This indicates whether blocks are formed from the highest LD pairs (i.e., initiate anywhere in a chromosome, "LD") or whether
#a left to right approach with blocking starting at the first SNP in the chromosome (by position, "beginning") is utilized.

#parallel - if true, set up a paralellization framework. For small marker sets (i.e., a few thousand or less), this will actually make the analysis slower. Parallelization is only
#recommended for many markers. Note: This only works for windows currently and WILL increase memory requirements! Only set to TRUE if it is working while on FALSE.

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



#everything is as before, but set_missing_NA is an option that, if set to TRUE, sets any localGEBV with >= 1 missing genotype
#to NA. If set to FALSE, it will compute localGEBV based on non-missing genotypes and only be NA if all genotypes are missing.

#mean_adjust = TRUE is the default and in the vast majority of cases needed. When estimating markers, USUALLY the genotype matrix is centered.
#Thus, marker effects are relative to centered genotypes. This is always true for GBLUP back solve and I believe
#is more or less always true for Bayesian alphabet models and rrBLUP. This is important to pay attention to
#you MUST get this right to prevent biased localGEBV! Note: it will not affect block variance, but it will change localGEBV values and affect the GA step.

#A good confirmation is to reconstruct GEBV from Zu, where Z is the centered marker matrix and u are the marker effects.
#The mean should be 0 (or very close to it) and reflect GBLUP GEBV if using rrBLUP, BayesC, or GBLUP back solve methods.
#If the reconstructed GEBV mean is meaningfully away from 0, it indicates the wrong marker matrix (i.e., needs to be centered)
#is being used or the marker matrix was not centered when estimating marker effects.

#The package will internally center markers if center = TRUE. If the matrix is already centered, centering won't change the values
#or centering can be set to FALSE.

#If you want genotype/haplotype configurations to match the 0/1/2 format, provide the uncentered genotype matrix and set center = TRUE.
#Otherwise, the reported genotype/haplotype configurations will reflect centered values. For clarity, you must provide dosage format
# (usually 0/1/2 for diploid, up to 9 for polypoid). We currently support polyploidy up to a ploidy of 9.

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
#max_kb specified how far out to compute LD decay and filters the LD structure, which helps to reduce computation time
#method specifies whether to fit a GAM thin-plate regression (gam_tp), GAM cubic spline regression (gam_cr), exponential decay (exp), or LOESS (loess) curve
#span controls the level of smoothing for the LOESS method - specify between 0 and 1
#closer to 0 = little to no smoothing, closer to 1 = very smoothed and likely monotonic
#k is the number of basis functions/knots in the gam options - experiment with values between 10 and 100
#the lower the number of k, the more smoothed it is
#exponential does not use either of the smoothing values and is guaranteed to be monotonic

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
head(haploblock_obj$Haploblocks_GA)
head(haploblock_obj$Haplotype_Effect_Matrix_GA)
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
