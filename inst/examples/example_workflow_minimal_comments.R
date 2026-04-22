#####################################
    ####  load dependencies  ####
#####################################

# Package functions and required dependencies are provided by HapSelect.

required_fns <- c(
  "order_map", "pairwise_ld", "def_blocks", "block_obj_to_df",
  "compute_local_GEBV", "marker_effects_plot", "unique_haplo_effects_plot",
  "block_var_funnel_plot", "plot_haploblocks", "plot_marker_density",
  "plot_ld_decay", "genetic_algorithm"
)

missing_fns <- required_fns[!vapply(required_fns, exists, logical(1), mode = "function")]
if (length(missing_fns) > 0) {
  stop(
    "Package functions are not loaded. Run devtools::load_all('.') ",
    "or library(FastStack) before sourcing this script."
  )
}


#################################
    ####  Pairwise LD  ####
#################################


#loading map file from an R object
data("map", package = "FastStack")

#simulate unordered map file - original file is already ordered, so this is just an example
map2 = map[sample(1:nrow(map), nrow(map)), ]
map2 = order_map(map = map2)

#loading the genotype file from an R object
data("geno", package = "FastStack")


#compute ld
#ld_pairs = pairwise_ld(geno, parallelize = FALSE)
#This is very slow and is just for demonstration purposes only!



#load the example file to see the structure if you do not want to run the function
data("pairwise_ld", package = "FastStack")
#We recommend using PLINK v1.9 for LD calculations.
#If you have a PLINK binary fileset (.bed/.bim/.fam), you can also use:
#ld_pairs = plink_pairwise_ld("path/to/plink_prefix")

#########################################
  #####  Compute Haploblocks  #####
#########################################

#haploblocking

haploblocks = def_blocks(ld = ld_pairs, map = map, method = "flanking",
                         threshold = 0.2, tolerance = 4, tol_reset = TRUE,
                         start = "LD", parallel = FALSE)

#turn the block object into a data frame
haploblocks = block_obj_to_df(haploblocks, map)

#######################################
    #### Compute localGEBV  #####
#######################################

#the marker pecov (prediction error covariance matrix) is only needed if you want to compute
#p-values of specific haplotype effects - very much theoretical at this point

#the first column in the marker effects file should be the SNP ID (same as the map file) and the second
#column should be the marker effect estimated from a model.
data("marker_effects", package = "FastStack")




#localGEBV calculation
haploblock_obj = compute_local_GEBV(geno = geno, marker_effects = marker_effects, haploblocks_df = haploblocks,
                                    set_missing_NA = TRUE, center = TRUE)
##################################
   ####  Visualizations  #####
##################################

marker_plot = marker_effects_plot(marker_effects = marker_effects$Effect, chr = map$Chromosome, pos = map$Position)
marker_plot

haplo_eff_plot = unique_haplo_effects_plot(haplo_obj = haploblock_obj)
haplo_eff_plot

funnel_plot = block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE)
funnel_plot

haploblock_plot = plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks)
haploblock_plot

marker_density_plot = plot_marker_density(map_df = map, bin_size = 500000)
marker_density_plot

ld_decay_plot = plot_ld_decay(map = map, ld = ld_pairs, max_kb = 500, span = 0.3, k = 10, method = "gam_cr")
ld_decay_plot


##############################################
   #### Parent Selection with the GA ####
##############################################

#select top 15 haploblocks (arbitrary) and perform the GA
haploblock_effects = haploblock_obj$Haploblocks
haploblock_effects = haploblock_effects[order(haploblock_effects$Block_Var, decreasing = TRUE), ]
haploblock_top_blocks = haploblock_effects[1:15, ]
localGEBV = haploblock_obj$Haplotype_Effect_Matrix
localGEBV = localGEBV[row.names(localGEBV) %in% haploblock_top_blocks$Block_ID, ]

localGEBV = as.data.frame(t(as.matrix(localGEBV)))


#run the GA
GA_output = genetic_algorithm(localGEBV = localGEBV, n_founders = 20, popSize = 10, maxiter = 300,
                              run = 150, selfing = FALSE, pmutation = 0.2, pcrossover = 0.8, pelite = 0.5)

#one unique set of parents
GA_output$One_Solution
