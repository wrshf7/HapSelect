library(dplyr)

####GAPIT maize example with 3k SNP####
#pheno = read.table("Example_Files/mdp_traits.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
geno = read.table("Example_Files/mdp_genotype_test.hmp.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
individuals_geno = geno[1,12:ncol(geno)]
geno = read.table("Example_Files/mdp_genotype_test.hmp.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
map = read.table("Example_Files/mdp_SNP_information.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

pheno = read.table("Example_Files/mdp_traits_validation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#turn genotypes into dosages
geno = geno[,c(1,3,4,12:ncol(geno))]
geno[geno == "NN"] = -9
geno[geno == "GG"] = "AA"
geno[geno == "CC"] = "TT"
geno[geno == "GC" | geno == "CG" | geno == "AG" | geno == "GA" | geno == "GT" | geno == "TG" | geno == "CT" | geno == "TC" | geno == "AC" | geno == "CA" | geno == "GA" | geno == "AG"] = "AT"
geno[geno == "AA"] = 2
geno[geno == "AT" | geno == "TA"] = 1
geno[geno == "TT"] = 0
geno[,4:ncol(geno)] = lapply(geno[,4:ncol(geno)], as.numeric)

#do some cleaning and structuring
row.names(geno) = geno[,1]
geno_mat = geno[,4:ncol(geno)]
geno_mat = t(geno_mat)

#compute maf and mean for imputation and snp filtering
maf = apply(geno_mat, 2, function(x){
  num = length(x[x != -9])
  maf = (length(x[x == 0]) + 0.5 * length(x[x==1])) / num * 100
  mean_dose = mean(x[x != -9], na.rm = TRUE)
  het = length(x[x==1]) / num * 100
  return_df = data.frame(
    maf = maf,
    mean = mean_dose,
    het = het
  )
}) %>% do.call(rbind, .)

maf$homo = 100 - maf$het

#impute missing values with the mean (does not affect LD or genomic relationships)
for(i in 1:nrow(geno)){
  geno[i, which(geno[i,] == -9)] = maf[i,"mean"]
  geno_mat[which(geno_mat[,i] == -9),i] = maf[i, "mean"]
}

#finish maf calcuation then filter for maf
maf$maf = ifelse(maf$maf > 50, 100 - maf$maf, maf$maf)

geno = geno[maf$maf > 5, ]
geno_mat = geno_mat[,maf$maf > 5]

map = map[map$SNP %in% row.names(geno), ]

markers = map$SNP
individuals = pheno$taxa
individuals_geno = lapply(individuals_geno, function(x)x) %>% do.call(c,.)



#remove markers of individuals without phenos and vice-versa - no longer needed with the validation file
# individuals_geno_remove = individuals_geno %in% individuals
# geno = geno[,c(1:3,(which(individuals_geno_remove) + 3))]
# geno_mat = geno_mat[individuals_geno_remove, ]

# #organize phenotype file - no longer needed with the validation file
# pheno = pheno[pheno$Taxa %in% individuals_geno, ]
# pheno = pheno[match(pheno$Taxa, individuals_geno), ]

#map file is in the same order as the genotype file

save(pheno, file = "Example_Files/gapit_pheno.R")
save(geno, file = "Example_Files/gapit_genos.R")
save(geno_mat, file = "Example_Files/gapit_geno_mat.R")
save(map, file = "Example_Files/gapit_map.R")
save(markers, file = "Example_Files/markers.R")
save(individuals_geno, file = "Example_Files/individual_order.R")

#####original chickpea files, but have issues with population structure#####
# load("Example_Files/geno.R")
# load("Example_Files/map.R")
# phenos = read.table("Example_Files/BLUEs.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
# 
# geno_G_mat = geno
# 
# geno_G_mat = geno_G_mat[,-2]
# row.names(geno_G_mat) = geno_G_mat$V1
# geno_G_mat = geno_G_mat[,-1]
# 
# markers = row.names(geno_G_mat)
# individuals = colnames(geno_G_mat)
# 
# geno_G_mat = t(geno_G_mat)
# 
# #already in the same order, so no need to order
# #assuming that all filtering has already been done
# G_mat = Gmatrix(SNPmatrix = geno_G_mat, method = "VanRaden", missingValue = -9, maf = 0, thresh.missing = 0.5, ploidy = 2, integer = TRUE, scale = TRUE)
# G_inv = solve(G_mat)
# 
# wrong formulas - need to multiply the 1s by 0.5
# test = apply(geno_G_mat, 1, function(x){
#   x = sum(x == 0 | x == 1) / length(x) * 100
#   return(x)
# }) #%>% do.call(c,.)
# 
# wrong formulas - need to multiply the 1s by 0.5
# test2 = apply(geno_G_mat, 2, function(x){
#   x = sum(x == 0 | x == 1) / length(x) * 100
#   return(x)
# })
# 
# test2 = ifelse(test2 > 50, 100 - test2, test2)
# 
# geno_G_mat2 = geno_G_mat[,test2 > 10]
# 
# G_mat2 = Gmatrix(SNPmatrix = geno_G_mat2, method = "Yang", missingValue = -9, maf = 0, thresh.missing = 0.5, ploidy = 2, integer = 2)
# 
# vanraden2 <- function(geno_matrix) {
#   # geno_matrix: individuals x SNPs, coded 0/1/2
#   p <- colMeans(geno_matrix, na.rm = TRUE) / 2
#   Z <- sweep(geno_matrix, 2, 2 * p)             # center: M - 2p
#   denom <- sqrt(2 * p * (1 - p))                # std dev per SNP
#   Z_std <- sweep(Z, 2, denom, FUN = "/")        # scale columns
#   G <- tcrossprod(Z_std) / ncol(geno_matrix)    # average dot product
#   return(G)
# }
# 
# G_mat = vanraden2(geno_G_mat)
# G_mat2 = vanraden2(geno_G_mat2)
