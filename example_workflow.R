####load dependencies####
library(progressr)
library(furrr)
library(future)
library(parallel)
library(tidyverse)
library(cowplot)
library(GA)

#set working directory to the location of the scripts
#setwd("")

#load scripts
source("order_map_file.R")
source("pairwise_LD.R")
source("def_haploblocks.R")
source("localGEBV_calculation.R")
source("visualizations.R")
source("genetic_algorithm.R")

######Pairwise LD#########

#loading map file from an R object
#the map file chromosomes MUST be numeric!
#First column should be SNP ID, second column should be chromosome, and third
#column should be the position
load(file = "Example_Files/gapit_map.R")

#simulate unordered map file - original file is already ordered, so this is just an example
map2 = map[sample(1:nrow(map), nrow(map)), ]

#check that the map file is ordered - using simulated unordered map
#this function will also check the map file to ensure it is appropriate and set column names
map2 = order_map(map = map2)

#loading the genotype file from an R object

load(file = "Example_Files/gapit_genos.R")

#compute ld

#first 3 columns of the genotyep file should be the map file and the rest of the columns represent accessions/individuals with cells
#being the genotypes at each locus for each individual
#column names of the map information (first 3 columns) will not matter in this case, but the order of map info columns
#should be the same as the map file columns (SNP ID, Chromosome, Position)
ld_pairs = pairwise_ld(geno, parallelize = FALSE)

#ld_pairs = pairwise_ld(geno, parallelize = TRUE)


#load the example file to see the structure if you do not want to run the function
load("Example_Files/gapit_ld.R")
ld_pairs = gapit_pairwise_ld

#note: other programs can be utilized to generate the LD file, but make sure the columns c("Chrom", "Locus1", "Locus2", "Name1", "Name2", "LD")
#exist in the data frame

#Locus1 and Locus2 are numeric indices and MUST reference the order of the markers (i.e., ordered by chromosome and position in chromosome)
#A good way to do this is to assign Locus1 as 1:nrow(map) to the ordered map file. You can then join the LD data frame and
#map data frame by Name1 and SNP and them Name2 and SNP to get the indices for Locus1 and Locus2.

#####make haploblocks#####

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
#recommended for many markers.

#haploblocking - using a threshold of 0.4, a tolerance of 2 as an example, 
#resetting the tolerance between successful marker additions,
#forming blocks around the highest LD pairs, and using LD between adjacent markers (flanking)
#(i.e., not average LD to the block)

haploblocks = def_blocks(ld = ld_pairs, map = map, method = "flanking",
                          threshold = 0.2, tolerance = 4, tol_reset = TRUE, start = "LD", parallel = FALSE)

#turn the block object into a data frame
haploblocks = block_obj_to_df(haploblocks, map)


##### compute localGEBV using marker effects #####

#the marker pecov (prediction error covariance matrix) is only needed if you want to compute
#p-values of specific haplotype effects - very much theoretical at this point
#the marker pecov structure is a matrix of prediction error (co)variances between the markers.
#The rows and columns indicate markers and should be in the same order as marker effects

#the first column in the marker effects file should be the SNP ID (same as the map file) and the second
#column should be the marker effect estimated from a model.
load(file = "Example_Files/gapit_marker_pecov.R")
load(file = "Example_Files/gapit_marker_effects.R")

#ideally you just set missing values to NA. This would be fine for LD calculations,
#calculating marker effects, and computing localGEBV
#I'm just keeping it simple for this example and rounding so there are no missing values
#Non-integers represent imputed values, which must NOT be used for localGEBV!
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

#marker density plot
marker_density_plot = plot_marker_density(map_df = map, bin_size_kb = 500)

#LD decay plot
#plot the LD decay utilizing the map and LD objects
#max_kb specified how far out to compute LD decay and filters the LD structure, which helps to reduce computation time
#method specifies whether to fit a GAM thin-plate regression (gam_tp), GAM cubic spline regression (gam_cr), expontential decay (exp), or LOESS (loess) curve
#span controls the level of smoothing for the LOESS method - specify between 0 and 1
#closer to 0 = little to no smoothing, closer to 1 = very smoothed and likely monotonic
#k is the number of basis functions/knots in the gam options - experiment with values between 10 and 100
#the lower the number of k, the more smoothed it is
#exponential does not use either of the smoothing values and is guaranteed to be monotonic

ld_decay_plot = plot_ld_decay(map = map, ld = ld_pairs, max_kb = 500, span = 0.3, k = 10, method = "gam_cr")

###### select top 15 haploblocks (arbitrary) and perform the GA ######

haploblock_effects = haploblock_obj$Haploblocks

#order the blocks by block variance
haploblock_effects = haploblock_effects[order(haploblock_effects$Block_Var, decreasing = TRUE), ]

#select the top 15
haploblock_top_blocks = haploblock_effects[1:15, ]

#pull out the localGEBV corresponding to the top 15
localGEBV = haploblock_obj$Haplotype_Effect_Matrix
localGEBV = localGEBV[row.names(localGEBV) %in% haploblock_top_blocks$Block_ID, ]

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
GA_output = genetic_algorithm(localGEBV = localGEBV, n_founders = 20, popSize = 10, maxiter = 300, 
                              run = 150, selfing = FALSE, pmutation = 0.2, pcrossover = 0.8, pelite = 0.5)

#one unique set of parents - other solutions may exist, see the GA_output$GA@solution object for other potential sets of parents with the same fitness
GA_output$One_Solution