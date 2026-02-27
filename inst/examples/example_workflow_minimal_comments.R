####load dependencies####
# Package functions and required dependencies are provided by FastStack.

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

# Resolve bundled example data for both installed and development usage.
extdata_file <- function(name) {
  installed <- system.file("extdata", name, package = "FastStack")
  if (nzchar(installed)) {
    return(installed)
  }

  local <- file.path("inst", "extdata", name)
  if (file.exists(local)) {
    return(local)
  }

  stop(
    "Could not locate extdata file: ", name,
    ". Run from package root or install/load FastStack."
  )
}

######Pairwise LD#########

load(file = extdata_file("gapit_map.R"))

#simulate unordered map file - original file is already ordered, so this is just an example
map2 = map[sample(1:nrow(map), nrow(map)), ]

#check that the map file is ordered - using simulated unordered map
map2 = order_map(map = map2)

#loading the genotype file from an R object

load(file = extdata_file("gapit_genos.R"))

#compute ld

ld_pairs = pairwise_ld(geno, parallelize = FALSE)

#load the example file to see the structure if you do not want to run the function
load(extdata_file("gapit_ld.R"))
ld_pairs = gapit_pairwise_ld

#####make haploblocks#####

haploblocks = def_blocks(ld = ld_pairs, map = map, method = "flanking",
                         threshold = 0.2, tolerance = 4, tol_reset = TRUE, 
                         start = "LD", parallel = FALSE)

haploblocks = block_obj_to_df(haploblocks, map)


##### compute localGEBV using marker effects #####

load(file = extdata_file("gapit_marker_pecov.R"))
load(file = extdata_file("gapit_marker_effects.R"))

#with well-curated marker data, ideally you just set missing genotype values to NA. This would be fine for LD calculations,
#calculating marker effects, and computing localGEBV
#I'm just keeping it simple for this example and rounding so there are no missing values
#Non-integers represent imputed values, which must NOT be used for localGEBV!

geno[,4:ncol(geno)] = round(geno[,4:ncol(geno)], 0)

haploblock_obj = compute_local_GEBV(geno = geno, marker_effects = marker_effects, haploblocks_df = haploblocks, 
                                    marker_pecov = marker_pecov, set_missing_NA = TRUE, center = TRUE)

###### do some visualizations ######

marker_plot = marker_effects_plot(marker_effects = marker_effects$Effect, chr = map$Chromosome, pos = map$Position)
haplo_eff_plot = unique_haplo_effects_plot(haplo_obj = haploblock_obj)
funnel_plot = block_var_funnel_plot(haplo_obj = haploblock_obj, mean_line = FALSE)
haploblock_plot = plot_haploblocks(haploblock_df = haploblock_obj$Haploblocks)
marker_density_plot = plot_marker_density(map_df = map, bin_size_kb = 500)

ld_decay_plot = plot_ld_decay(map = map, ld = ld_pairs, max_kb = 500, span = 0.3, k = 10, method = "gam_cr")

###### select top 15 haploblocks (arbitrary) and perform the GA ######

haploblock_effects = haploblock_obj$Haploblocks
haploblock_effects = haploblock_effects[order(haploblock_effects$Block_Var, decreasing = TRUE), ]
haploblock_top_blocks = haploblock_effects[1:15, ]

localGEBV = haploblock_obj$Haplotype_Effect_Matrix
localGEBV = localGEBV[row.names(localGEBV) %in% haploblock_top_blocks$Block_ID, ]

localGEBV = as.data.frame(t(as.matrix(localGEBV)))


#run the GA
GA_output = genetic_algorithm(localGEBV = localGEBV, n_founders = 20, popSize = 10, maxiter = 300, 
                              run = 150, selfing = FALSE, pmutation = 0.2, pcrossover = 0.8, pelite = 0.5)

GA_output$One_Solution
