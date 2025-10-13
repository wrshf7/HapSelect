######Pairwise LD#########

#loading map file from an R object
#the map file chromosomes MUST be numeric!
#First column should be SNP ID, second column should be chromosome, and third
#column should be the position
load(file = "Example_Files/gapit_map.R")

#simulate unordered map file
map2 = map[sample(1:nrow(map), nrow(map)), ]

#check that the map file is ordered - using simulated unordered map
map2 = order_map(map = map2)

#loading the genotype file from an R object
load(file = "Example_Files/gapit_genos.R")

#compute ld
ld_pairs = pairwise_ld(geno)

#haploblocking - using a threshold of 0.4, a tolerance of 2, 
#resetting the tolerance between successful marker additions,
#forming blocks around the highest LD pairs, and using LD between adjacent markers
#(i.e., not average LD to the block)
haploblocks = def_blocks(ld = ld_pairs, map = map, method = "flanking",
                          threshold = 0.4, tolerance = 0.2, tol_reset = TRUE, start = "LD")
haploblocks = block_obj_to_df(haploblocks, map)

#compute localGEBV using markers
#the marker pecov (prediction error covariance matrix) is only needed if you want to compute
#p-values of specific haplotype effects - very much theoretical at this point
load(file = "Example_Files/gapit_marker_pecov.R")
load(file = "Example_Files/gapit_marker_effects.R")

#ideally you just set missing values to NA. This is fine for LD calculations,
#calculating marker effects, and computing localGEBV
#I'm just keeping it simple for this example and rounding
geno[,4:ncol(geno)] = round(geno[,4:ncol(geno)], 0)

haploblock_obj = compute_local_GEBV(geno = geno, marker_effects = marker_effects, haploblocks_df = haploblocks, marker_pecov = marker_pecov)
