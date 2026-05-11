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
map = HapSelect::map


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
geno = HapSelect::geno

#compute ld

#first 3 columns of the genotype file should be the map file and the rest of the columns represent accessions/individuals with cells
#being the genotypes at each locus for each individual
#column names of the map information (first 3 columns) will not matter in this case, but the order of map info columns
#should be the same as the map file columns (SNP ID, Chromosome, Position)

#ld_pairs = pairwise_ld(geno, parallelize = FALSE)

#This is very slow and is just for demonstration purposes only!



#load the example file to see the structure if you do not want to run the function
ld_pairs = HapSelect::ld_pairs


####################################################
####  Pairwise LD with PLINK Call (recommended) ####
####################################################


#loading map file from an R object
#the map file chromosomes MUST be numeric!
#First column should be SNP ID, second column should be chromosome, and third
#column should be the position

map = HapSelect::map

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
geno = HapSelect::geno

#compute ld with PLINK - requires installation of PLINK v1.9 and PATH availability
ld_pairs = plink_pairwise_ld_geno(geno = geno, ld_window = 999999, ld_window_kb = 1e6, ld_window_r2 = 0)

###########################################
##### Load LD File (Skip Computation) #####
###########################################

#load the example file to see the structure if you do not want to run the function
ld_pairs = HapSelect::ld_pairs


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
#### Compute Marker Effects  #####
#######################################

#Basic genomic prediction to obtain marker effects - This assumes one phenotype/BLUE/BLUP per individual.
#No other effects in the model allowed. For more advanced modeling, use other modeling software to obtain marker effects.

#The BLUE file should be structured such that the first column is comprised of individuals and the second column is comprised of a singular phenotype, BLUE, or de-regressed BLUP for each individual
BLUE = HapSelect::BLUE

#Basic genomic prediction using rrBLUP to solve for marker effects
# geno:      the genotype matrix previously utilized above
# BLUE:      the data frame containing the singular phenotype/BLUE/de-regressed BLUP for each individual as described above
# h2_method: scaling method to convert marker variance to additive genetic variance to compute narrow-sense heritability.
#              For most users, the default "VanRaden" method should be utilized for agricultural species. This is equivalent to GBLUP in VanRaden's method 1.
#              "VanRaden" uses 2 * sum(p*(1-p)) as a multiplier to scale the marker variance, whereas "marker_num" simply scales the marker variance by multiplying by the number of markers.
#              Using "marker_num" is generally only appropriate if the genotype matrix provided is already scaled.
# ploidy:    the integer value specifying the ploidy of the species. This must be an integer, which has the format of #L in R. If no L follows the number, R will not treat it as an integer and the internal check will fail!


marker_effects = create_marker_effects_file(geno = geno, BLUE = BLUE, h2_method = "VanRaden", ploidy = 2L)

#Basic functions to assess model performance:

#Basic n fold cross-validation
# geno:      the genotype matrix previously utilized above
# BLUE:      the data frame containing the singular phenotype/BLUE/de-regressed BLUP for each individual as described above
# h2_method: scaling method to convert marker variance to additive genetic variance to compute narrow-sense heritability.
#              For most users, the default "VanRaden" method should be utilized for agricultural species. This is equivalent to GBLUP in VanRaden's method 1.
#              "VanRaden" uses 2 * sum(p*(1-p)) as a multiplier to scale the marker variance, whereas "marker_num" simply scales the marker variance by multiplying by the number of markers.
#              Using "marker_num" is generally only appropriate if the genotype matrix provided is already scaled.
# ploidy:    the integer value specifying the ploidy of the species. This must be an integer, which has the format of #L in R. If no L follows the number, R will not treat it as an integer and the internal check will fail!
# nfold:     the number of folds in the cross-validation (e.g., how many groups the data is split into with all but one group used for training and each group used for validation once)

CV = n_fold_cross_validation(geno = geno, BLUE = BLUE, nfold = 5L, h2_method = "VanRaden", ploidy = 2L)

#Basic cross-validation
# geno:       the genotype matrix previously utilized above
# BLUE:       the data frame containing the singular phenotype/BLUE/de-regressed BLUP for each individual as described above
# h2_method:  scaling method to convert marker variance to additive genetic variance to compute narrow-sense heritability.
#               For most users, the default "VanRaden" method should be utilized for agricultural species. This is equivalent to GBLUP in VanRaden's method 1.
#               "VanRaden" uses 2 * sum(p*(1-p)) as a multiplier to scale the marker variance, whereas "marker_num" simply scales the marker variance by multiplying by the number of markers.
#               Using "marker_num" is generally only appropriate if the genotype matrix provided is already scaled.
# ploidy:     the integer value specifying the ploidy of the species. This must be an integer, which has the format of #L in R. If no L follows the number, R will not treat it as an integer and the internal check will fail!
# train_prop: the random proportion of individuals used for training (between 0 and 1) with the remainder utilized for validation
# fold:       the number of iterations to randomly sample train_prop proportion of individuals and compute prediction accuracy on the remainder validation set

CV = cross_validation(geno = geno, BLUE = BLUE, train_prop = 0.9, fold = 5L, h2_method = "VanRaden", ploidy = 2L)

#######################################
#### Compute localGEBV  #####
#######################################

#the first column in the marker effects file should be the SNP ID (same as the map file) and the second
#column should be the marker effect estimated from a model.
marker_effects = HapSelect::marker_effects



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

#The package will internally center markers if mean_adjust = TRUE. If the matrix is already centered, centering won't change the values
#or centering can be set to FALSE.

#If you want genotype/haplotype configurations to match the 0/1/2 format, provide the uncentered genotype matrix and set mean_adjust = TRUE.
#Otherwise, the reported genotype/haplotype configurations will reflect centered values. For clarity, you must provide dosage format
# (usually 0/1/2 for diploid, up to any integer for polypoid).

haploblock_obj = compute_local_GEBV(geno = geno, marker_effects = marker_effects, haploblocks_df = haploblocks,
                                    set_missing_NA = TRUE, mean_adjust = TRUE, parallel = TRUE)

##################################
####  Visualizations  #####
##################################

#must provide the column of the marker effects as well as the column indicating the chromosome and position on the chromosome
#if the map and marker effects plot are in the same order, you can do it as it is done here, else merge the two data frames together
# marker_effects: the vector of marker effects to plot
# chr:            the vector of chromosomal values for the marker effects
# pos:            the vector of positional values for the marker effects
# colors:         alternating chromosomal colors - only 2 valid R colors or hexidecimal colors should be provided
marker_plot = marker_effects_plot(marker_effects = marker_effects$Effect, chr = map$Chromosome, pos = map$Position, colors = c("#A01FF0", "#A7A8AA"))
marker_plot

#all of these plots utilize the haploblock_obj created
#plots the unique haplotype (localGEBV) effects, similar to the marker effect plot
# colors:   similar to the marker effects - alternating chromosomal colors. Two valid colors must be provided.
# pos_type: valid values are "midpoint" or "start" - specifies whether the haploblock's haplotype effects are positioned at the midpoint or at the start of the haploblock
haplo_eff_plot = unique_haplo_effects_plot(haplo_obj = haploblock_obj, colors = c("#A01FF0", "#A7A8AA"), pos_type = "midpoint")
haplo_eff_plot

#funnel plot scales the block variance to be between the value of 0 and 1, otherwise a quadratic term could scale exponentially
#x-axis contains haplotype effects and lets the user compare haplotype effects to variance
# mean_line:     if TRUE adds a line to the block variance (scaled)
# scale_colors:  the color scale dictating the coloring of haplotype effects from low to high (x-axis). Must be three valid R/hexidecimal colors.
funnel_plot = block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE, scale_colors = c("blue", "purple", "red"))
funnel_plot

#plots the location of the haploblocks on the chromosome
# block_fill:      the color of the haploblocks
# height:          controls the thickness of the chromosomes - lower to make skinnier chromosomes and increase to make thicker chromosomes.
# chrom_fill:      the color of the chromosome background: leave as NA to make it transparent.
# single_width_bp: a mostly defunct option, leave as NULL (does not affect single marker segments currently)
haploblock_plot = plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks, block_fill = "#A01FF0", chrom_fill = NA,
                                   height = 0.30,
                                   single_width_bp = NULL)
haploblock_plot

#marker density plot
# map_df:     the map file utilized throughout as described above
# bin_size:   the size of chromosome segments to count markers and graph: default is 500 kb
# height:     the thickness of the chromosomes - lower to make skinnier chromosomes and increase to make thicker chromosomes.
# chrom_fill: the color of the chromosome background: leave as NA to make it transparent.
# col_*:      color range from the lowest value to the highest value used to make a sliding color scale.
marker_density_plot = plot_marker_density(map = map, bin_size = 500e3, height = 0.3,
                                          chrom_fill = NA,
                                          col_low = "white", col_mid = "purple", col_high = "red")
marker_density_plot

#LD decay plot
#plot the LD decay utilizing the map and LD objects
#map:         map dataframe object utilized throghout
#ld:          the ld dataframe object utilized to construct haploblocks
#max_kb:      specified how far out to compute LD decay (in kb) and filters the LD structure, which helps to reduce computation time
#method:      specifies whether to fit a GAM thin-plate regression ("gam_tp"), GAM cubic spline regression ("gam_cr"),
#               exponential decay ("exp"), or LOESS ("loess") curve
#span:        controls the level of smoothing for the LOESS method - specify between 0 and 1
#               closer to 0 = little to no smoothing, closer to 1 = very smoothed and likely monotonic
#k:           is the number of basis functions/knots in the gam options - experiment with values between 10 and 100 (default 50)
#               the lower the number of k, the more smoothed it is
#               Note: exponential does not use either of the smoothing values and is guaranteed to be monotonic
#point_color: color of the individual LD data points, must provide a singular, valid R/hexidecimal color
#curve_color: color of the computed LD curve - must be a singular, valid color similar to point_color
#alpha:       the transparency of the individual LD data points - must be between 0 and 1 (default 0.2)
ld_decay_plot = plot_ld_decay(map = map, ld = ld_pairs, max_kb = 500, span = 0.3, k = 10L, method = "gam_cr", point_color = "#A7A8AA",
                              curve_color = "#A01FF0", alpha = 0.2)
ld_decay_plot


###########################################################
#### Parent Selection with the GA and Basic Simulation ####
###########################################################

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


#basic simulation of GA vs TS selected parents and PCA plot

#The output provides a simuation plot showing potential genetic plateaus of GA vs TS selected parents,
# a dataframe containing PCA information (if PCA = TRUE), and a PCA plot showing where the GS and TS selected parents (and their overlap) sit in the overall diversity (if PCA = TRUE).

# The primary inputs are the same as previous steps or are generated by previous steps, including the genotype dataframe object, marker effects data frame, map data frame, and GA object just generated
# genetic_map_position: an optional argument if the genetic map of the markers is known. The genetic map should be a vector (singular column of a data frame) of equal length to the map and in the same order. Values should be in centiMorgans (cM).
#                       If left as NULL, each chromosome is assumed to be of size max_cM_chr (default 100 cM). A value of 100 cM corresponds to an average of 1 recombination event per chromosome.
#                       If left as NULL, the first marker in a chromosome is assumed 0 cM and the last marker is assumed max_cM_chr (default 100 cM) and cM values for other markers are proportional to their physical distance between the first and last marker.
# max_cM_chr:           The maximum genetic map distance per chromosome if a genetic map is not provided (default is 100 centiMorgans; cM). A value of 100 cM corresponds to an average of 1 crossover (recombination) event per chromosome.
# num_gen:              The number of generations of recurrent truncation selection (TS) to simulate using both the GA selected and TS selected parents - default is 50.
# num_sim_reps:         The number of simulation replicates to conduct to quantify the +/- 2 standard deviations of genetic progress. Differences arise from recombination differences, which chromosome is passed to offspring, and alternative assignments of heterozygotes (see genomicSimulation documents for more detail).
# num_cross_per_gen:    The number of offspring generated from random mating between parents each generation. The default is 1000 progeny.
# num_TS_parents:       The number of truncation selected (TS) parents to utilize (based on whole-genome GEBV) and compare to the specified number of GA parents.
#                       If left as NULL, the value is assumed to be the same as the number of GA selected parents provided.
# mean_adjust:          Default is TRUE - centers the genotype matrix (see Van Raden, 2008) to prevent a mean shift (bias) in GEBV. This will not influence rankings or the rate of genetic progress, but will inflate/deflate average GEBV by a constant.
# PCA:                  If TRUE (default), conducts a PCA on the genotype matrix and highlights the parents selected from the GA, TS, and the overlap using colors specified in colors.
# colors:               Vector of length 4 of colors or hexidecimal color values recognized by R.
#                       The four values highlight GA selected parents, TS selected parents, the parents selected by both methods, and individuals not selected as parents, respectively.
# alpha:                Vector of length 4 with values between 0 and 1 to set the transparency of the four parent categories in colors. Default is no transparency for any selected parent and 50% transparency for individuals not selected by any method.
parent_sln_obj = GA_vs_TS_simulation(GA_output = GA_output, geno = geno, marker_effects = marker_effects, map = map, genetic_map_position = NULL, num_gen = 50, num_sim_reps = 30,
                                     num_cross_per_gen = 1000, num_TS_parents = NULL, mean_adjust = TRUE, max_cM_chr = 100, PCA = TRUE,
                                     colors = c("green", "#d95f02", "#A01FF0", "gray80"), alpha = c(1,1,1,0.5))

#display plots
parent_sln_obj$Simulation_Plot
parent_sln_obj$PCA_Plot
