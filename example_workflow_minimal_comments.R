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
load(file = "Example_Files/gapit_map.R")

#simulate unordered map file - original file is already ordered, so this is just an example
map2 = map[sample(1:nrow(map), nrow(map)), ]

#check that the map file is ordered - using simulated unordered map
#this function will also check the map file to ensure it is appropriate and set column names
map2 = order_map(map = map2)

#loading the genotype file from an R object
load(file = "Example_Files/gapit_genos.R")

#compute ld
ld_pairs = pairwise_ld(geno)
#ld_pairs = pairwise_ld(geno, parallelize = TRUE)


#load the example file to see the structure if you do not want to run the function
#load("Example_Files/gapit_ld.R")

#####make haploblocks#####

#haploblocking - using a threshold of 0.4, a tolerance of 2, 
#resetting the tolerance between successful marker additions,
#forming blocks around the highest LD pairs, and using LD between adjacent markers (flanking)
#small example, so not parallelizing
haploblocks = def_blocks(ld = ld_pairs, map = map, method = "flanking",
                         threshold = 0.2, tolerance = 4, tol_reset = TRUE, start = "LD", parallel = FALSE)

#turn the block object into a data frame
haploblocks = block_obj_to_df(haploblocks, map)


##### compute localGEBV using marker effects #####

#the marker pecov (prediction error covariance matrix) is only needed if you want to compute
#p-values of specific haplotype effects - very much theoretical at this point

load(file = "Example_Files/gapit_marker_pecov.R")
load(file = "Example_Files/gapit_marker_effects.R")

#I'm rounding imputed values to integers - more ideal to set to NA
geno[,4:ncol(geno)] = round(geno[,4:ncol(geno)], 0)

haploblock_obj = compute_local_GEBV(geno = geno, marker_effects = marker_effects, haploblocks_df = haploblocks, marker_pecov = marker_pecov)

###### do some visualizations ######

marker_plot = marker_effects_plot(marker_effects = marker_effects$Effect, chr = map$Chromosome, pos = map$Position)
haplo_eff_plot = unique_haplo_effects_plot(haplo_obj = haploblock_obj)
funnel_plot = block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE)
haploblock_plot = plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks)

###### select top 15 haploblocks (arbitrary) and perform the GA ######

haploblock_effects = haploblock_obj$Haploblocks
haploblock_effects = haploblock_effects[order(haploblock_effects$Block_Var, decreasing = TRUE), ]
haploblock_top_blocks = haploblock_effects[1:15, ]
localGEBV = haploblock_obj$Haplotype_Effect_Matrix
localGEBV = localGEBV[row.names(localGEBV) %in% haploblock_effects$Block_ID, ]
localGEBV = as.data.frame(t(as.matrix(localGEBV)))


#GA with some example parameters
GA_output = genetic_algorithm(localGEBV = localGEBV, n_founders = 20, popSize = 10, maxiter = 300, 
                              run = 150, selfing = FALSE, pmutation = 0.2, pcrossover = 0.8, pelite = 0.5)

GA_output$One_Solution
