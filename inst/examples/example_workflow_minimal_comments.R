##################################################
####  Example Workflow with Minimal Comments  ####
##################################################

library(HapSelect)

#############################################################
####  Pairwise LD with Native Function (not recommended) ####
#############################################################

map = HapSelect::map

#simulate unordered map file - original file is already ordered, so this is just an example
map2 = map[sample(1:nrow(map), nrow(map)), ]

#check that the map file is ordered - using simulated unordered map
map2 = order_map(map = map2)

#loading the genotype file from an R object
#first 3 columns of the genotype file should be the map file and the rest of the columns represent accessions/individuals with cells
#being the genotypes at each locus for each individual

geno = HapSelect::geno


#compute ld (painfully slow)
#ld_pairs = pairwise_ld(geno, parallelize = FALSE)

#This is very slow and is just for demonstration purposes only!



#load the example file to see the structure if you do not want to run the function
ld_pairs = HapSelect::ld_pairs

#We recommend using PLINK v1.9 for LD calculations.
#If you have a PLINK binary fileset (.bed/.bim/.fam), you can also use:
#ld_pairs = plink_pairwise_ld("path/to/plink_prefix")
#see below for more details



####################################################
####  Pairwise LD with PLINK Call (recommended) ####
####################################################

map = HapSelect::map
geno = HapSelect::geno


#compute ld with PLINK - requires installation of PLINK v1.9 and PATH availability
ld_pairs = plink_pairwise_ld_geno(geno = geno, ld_window = 999999, ld_window_kb = 1e6, ld_window_r2 = 0)

###########################################
##### Load LD File (Skip Computation) #####
###########################################

#load the example file to see the structure if you do not want to run the function
ld_pairs = HapSelect::ld_pairs

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
#### Compute Marker Effects  #####
#######################################

#Basic genomic prediction to obtain marker effects - This assumes one phenotype/BLUE/BLUP per individual.
#No other effects in the model allowed. For more advanced modeling, use other modeling software to obtain marker effects.

#The BLUE file should be structured such that the first column is comprised of individuals and the second column is comprised of a singular phenotype, BLUE, or de-regressed BLUP for each individual
BLUE = HapSelect::BLUE

#Basic genomic prediction using rrBLUP to solve for marker effects
marker_effects = create_marker_effects_file(geno = geno, BLUE = BLUE, h2_method = "VanRaden", ploidy = 2L)

#Basic functions to assess model performance:

#Basic n fold cross-validation
CV = n_fold_cross_validation(geno = geno, BLUE = BLUE, nfold = 5L, h2_method = "VanRaden", ploidy = 2L)

#Basic cross-validation
CV = cross_validation(geno = geno, BLUE = BLUE, train_prop = 0.9, fold = 5L, h2_method = "VanRaden", ploidy = 2L)


#######################################
#### Compute localGEBV  #####
#######################################


#the first column in the marker effects file should be the SNP ID (same as the map file) and the second
#column should be the marker effect estimated from a model.
marker_effects = HapSelect::marker_effects

#compute localGEBV
haploblock_obj = compute_local_GEBV(geno = geno, marker_effects = marker_effects, haploblocks_df = haploblocks,
                                    set_missing_NA = TRUE, mean_adjust = TRUE)
##################################
####  Visualizations  #####
##################################

#must provide the column of the marker effects as well as the column indicating the chromosome and position on the chromosome
#if the map and marker effects plot are in the same order, you can do it as it is done here, else merge the two data frames together
marker_plot = marker_effects_plot(marker_effects = marker_effects$Effect, chr = map$Chromosome, pos = map$Position, colors = c("#A01FF0", "#A7A8AA"))
marker_plot

#all of these plots utilize the haploblock_obj created
#plots the unique haplotype (localGEBV) effects, similar to the marker effect plot
haplo_eff_plot = unique_haplo_effects_plot(haplo_obj = haploblock_obj, colors = c("#A01FF0", "#A7A8AA"), pos_type = "midpoint")
haplo_eff_plot

#funnel plot scales the block variance to be between the value of 0 and 1, otherwise a quadratic term could scale exponentially
funnel_plot = block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE, scale_colors = c("blue", "purple", "red"))
funnel_plot

#plots the location of the haploblocks on the chromosome
haploblock_plot = plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks, block_fill = "#A01FF0", chrom_fill = NA,
                                   height = 0.30,
                                   single_width_bp = NULL)
haploblock_plot

#marker density plot
marker_density_plot = plot_marker_density(map = map, bin_size = 500e3, height = 0.3,
                                          chrom_fill = NA,
                                          col_low = "white", col_mid = "purple", col_high = "red")
marker_density_plot

#LD decay plot
#plot the LD decay utilizing the map and LD objects
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
nrow(haploblock_obj$Haploblocks_GA)

#2 select top 50% of haploblocks
haploblock_obj = select_top_blocks(haploblock_obj = haploblock_obj, perc_total = 0.5)
nrow(haploblock_obj$Haploblocks_GA)

#3 select the top blocks explaining at least 90% of the total block variance
haploblock_obj = select_top_blocks(haploblock_obj = haploblock_obj, perc_of_total_var = 0.9)
nrow(haploblock_obj$Haploblocks_GA)


#run the GA
GA_output = genetic_algorithm(localGEBV = haploblock_obj$Haplotype_Effect_Matrix_GA, n_founders = 20, popSize = 10, maxiter = 300,
                              run = 150, selfing = FALSE, pmutation = 0.2, pcrossover = 0.8, pelite = 0.5)

#one unique set of parents - other solutions may exist, see the GA_output$GA@solution object for other potential sets of parents with the same fitness
GA_output$One_Solution

#basic simulation of GA vs TS selected parents and PCA plot

parent_sln_obj = GA_vs_TS_simulation(GA_output = GA_output, geno = geno, marker_effects = marker_effects, map = map, genetic_map_position = NULL, num_gen = 50, num_sim_reps = 30,
                               num_cross_per_gen = 1000, num_TS_parents = NULL, mean_adjust = TRUE, max_cM_chr = 100, PCA = TRUE,
                               colors = c("green", "#d95f02", "#A01FF0", "gray80"), alpha = c(1,1,1,0.5))

parent_sln_obj$Simulation_Plot
parent_sln_obj$PCA_Plot
