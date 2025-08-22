#####load libraries and data####
library(dplyr)
library(AGHmatrix)
library(sommer)
library(mvtnorm)

load(file = "Example_Files/gapit_pheno.R")
load(file = "Example_Files/gapit_genos.R")
load(file = "Example_Files/gapit_geno_mat.R")
load(file = "Example_Files/gapit_map.R")
load(file = "Example_Files/markers.R")
load(file = "Example_Files/individual_order.R")

####make genomic relationship matrix####
#no need to filter - already filtered for maf

#make sure it is in 0/1/2 format

####Gmat using package####
G_mat = Gmatrix(SNPmatrix = geno_mat, method = "VanRaden", integer = FALSE)

#make sure names are the same as the order of the phenotype file
row.names(G_mat) = individuals_geno
colnames(G_mat) = individuals_geno

#calculate inverse
G_mat_inv = solve(G_mat)

#highly inbred, but expected with inbred maize panel.. a couple of extreme individuals
#some possible clones... off-diagonals greater than 1 or close to 1... could be the same inbreds most likely


####build Z matrix from VanRaden and also make sure G matrix matches####
#used in downstream analyses to estimate SNP effects from EBV

#allele frequencies of each marker
p = apply(geno_mat, 2, function(x){
  #reference allele frequency in the dosage - homozygotes are 2's and heterozygotes are 1's
  #sum the actual values then divide by the total number of alleles (2 * individuals)
  p = sum(x) / (2*length(x))
  return(p)
})



#compute Z (Zi = Mi - 2 * pi, where Zi is the ith column or marker of Z, Mi is the ith column
#or marker of the marker matrix in 0/1/2 format, and pi is the allele frequency of the ith column/marker i)
Z = geno_mat
for(i in 1:ncol(Z)){
  #for a given marker, take it and subtract 2 * it's ref allele freq previously calculated
  Z[,i] = Z[,i] - 2*p[i]
}

#calculate the total marker binomial variance (denominator of Van Raden)
sum2pq = sum(2*p*(1-p))

#test that Z is correct by calculating and comparing G to AGHMatrix
G_manual = (Z %*% t(Z)) / sum2pq
colnames(G_manual) = individuals_geno
row.names(G_manual) = individuals_geno

#stabilize G so it's not singular (equivalent to a small blending to an identity matrix in BLUPF90)
diag(G_manual) = diag(G_manual) + 0.01

#invert G for downstream applications
G_man_inv = solve(G_manual)

#compare all of the values in the two G matrices than coerce to a vector to summarize differences
G_similarity = lapply(as.data.frame(G_manual - G_mat), function(x) x) %>% do.call(c, .)
summary(G_similarity)
summary(diag(G_manual) - diag(G_mat))
summary(lapply((lower.tri(G_manual) - lower.tri(G_mat)), function(x) x) %>% do.call(c,.))

#slight differences on the absolute minimum in the diagonals, but otherwise the diagonal values are nearly identical
#might have to do with how some scaling or stabilization is conducted. Z is probably correct
#off-diagonals are identical


####Fit sommer models and extract components#####
#newer sommer model, but does not provide the full PECOV structure (just PEV)
#this is sorta like ASREML-R - it can handle Factor Analytics and time-series models
# sommer_model = mmes(Sim100 ~ 1,
#                     random = ~ vsm(ism(taxa), Gu = G_test),
#                     rcov = ~ units, 
#                     nIters = 10,
#                     data = pheno, 
#                     verbose = TRUE, 
#                     getPEV = TRUE
#                   )

#recently deprecated sommer model for a univariate example, which provides the whole PECOV structure, not as efficient
#but provides what I need

#Sim100 is the phenotype
#No fixed effects - only an overall intercept/mean, so ~ 1
#only random effect is the additive genetic effect. vsr specifies(taxa, Gu - G_manual): taxa are the genotype names and Gu is the covariance structure
#rcov = ~ units indicates that we are assuming an identity covariance structure for the residuals (I*sigma^2_e)
#data = pheno specifies the data frame
#nIters specifies how many maximum iterations will be run
#It will terminate before if convergence criteria are reached
#getPEV = TRUE is only needed if calculating SNP effect
sommer_model = mmer(Sim100 ~ 1,
                    random = ~ vsr(taxa, Gu = G_manual),
                    rcov = ~ units, 
                    nIters = 10,
                    data = pheno, 
                    verbose = TRUE, 
                    getPEV = TRUE,
                    method = "AI"
)


#multi-trait mode for two traits
#note that cbind(Trait1, Trait2) is how I specify there are two traits.
#In this case it will be your yield in 2 environments. You can add on a third for the third environment
#Note I'm also saying there are no fixed effects, so I'm just fitting the overall mean for each environment 
#that is indicated by ~ 1

#my only random effect is the accession the random = ~vsr() lets me fit a random effect
#the random effect is the Animal column, and it's covariance structure between accessions
#is given by my kinship matrix, G_mat
#the Gtc = unsm(2) is how I get the genetic correlation between both traits (environments)
#The unsm(2) function tells me I have my two traits and I need to estimate the variances for both
#and the covariance between them
#If you have 3 traits, it will be unsm(3)
#Because I have two traits, I have had to expand my residuals.
#The rcov = ~ vsr() is again letting me specify that I have a random effect beyond the simple: I*sigma^2_e
#What I now have if a residual variance for each trait, which is given by Gtc = diag(2)
#The diag(2) tells the model to estimate a residual variance for each trait, but my residuals should be independent across traits
#In your case I also believe there is no reason to suggest environments are related, but you can always fit Gtc = unsm(2)
#Same thing with random, if you have 3 traits, then change diag(2) or unsm(2) to diag(3) or unsm(3)

sommer_basePop = mmer(cbind(Trait1, Trait2) ~ 1,
                      random = ~ vsr(Animal, Gu = G_mat, Gtc = unsm(2)),
                      rcov = ~ vsr(units, Gtc=diag(2)), nIters = 3,
                      data = phenotypes_sommer, verbose = TRUE, 
                      method = "AI", getPEV = TRUE
)

#example where I add in a second random effect named "Block" and assume that blocks (e.g., fields, rows. columns, etc.) across traits are unrelated
sommer_basePop = mmer(
  cbind(Trait1, Trait2) ~ 1,
  random = ~ vsr(Animal, Gu = G_mat, Gtc = unsm(2)) +
    vsr(Block, Gtc = diag(2)),
  rcov   = ~ vsr(units, Gtc = diag(2)),
  nIters = 3,
  data   = phenotypes_sommer,
  verbose = TRUE, 
  method = "AI", 
  getPEV = TRUE
)










#differences between G_test and AGHMatrix Gmatrix G are very slight, slightly
#inflated variances, but lambda remains the same pev differences are mild, slight bias
#in breeding values

#extract breeding values
blup = sommer_model$U$`u:taxa`$Sim100

#extract prediction error covariance (PECOV) structure
pecov = sommer_model$PevU$`u:taxa`$Sim100

#variances and lambda
gen_var = sommer_model$sigma$`u:taxa` %>% lapply(.,function(x)x) %>% do.call(c,.)
res_var = sommer_model$sigma$units %>% lapply(.,function(x)x) %>% do.call(c,.)
lambda = res_var / gen_var

####Back calculate marker effects####
#back calculate marker effects - it will be slightly off using AGHMatrix if their Z
#is manipulated before obtaining G. Technically, I should account for the blending
#see Andre Legarra's papers (Aguilar 2019 I believe) or online book for more information

#scaled Z matrix (before I just centered) multiplied into the inverse of G - this is putting 
#G on the level of markers
A = (1/sum2pq) * t(Z) %*% G_man_inv

#marker effects
marker_effects = A %*% blup

#reconstruct the original breeding values
reconstructed_bv = Z %*% marker_effects #to sanity check with original blups
summary(reconstructed_bv - blup) #very slight differences except at extremes, might indicate distributional problems


####Obtain marker PECOV, that is PECOV(m-hat)####
#marker PECOV - the proof here is that Var(u), that is the variance of the true breeding values,
#is equal to G * sigma^2_u (standard), which is equal to Var(Zm), where m is the true marker effects
#remember, Z*m gives the blup without scaling (with scaling, it is equivalent to my "A" above)
#so, we know that var(u) = var(Zm) = G*sigma^2_u = Z*var(m)*Z'
#Let's use G*sigma^2_u = Zvar(m)Z' = ZZ'* sigma^2_m 
#(identity implied here, markers assumed independent... future exploration)
#Let G = 1/sum2pq * ZZ'
#then ZZ' = sum2pq * G
#back to the equation G*sigma^2_u = ZZ' * sigma^2_m
#G*sigma^2_u = G*sum2pq * sigma^2_m  (I replaced ZZ' with sum2pq*G)
#the G on each side "cancels" if you multiply each side by the inverse of G
#sigma^2_u = sigma^2_m * sum2pq
#sigma^2_u / sum2pq = sigma^2_m
#this is the proof to calculating the true marker variance and the one usually used
marker_gen_var = gen_var / sum2pq

#another identity:
#var(m) = var(m-hat) + PEV(m-hat), where var(m-hat) is the variance of the marker estimates
#per above, var(m) =Isigma^2_m    (we are assuming that markers are independent - more on that later)
#sigma^2_m = var(m-hat) + PECOV(m-hat)
#PECOV(m-hat) = sigma^2_m *I - var(m-hat)
#var(m-hat) = var(1/sum2pq*Z*G^-1*(u-hat))
#let 1/sum2pq * Z * G^-1 = A
#then, var(m-hat) = A * var(u-hat) * A'
#var(m-hat) = A * (var(u) - PECOV(u-hat)) * A'
#var(m-hat) = A * (G*sigma^2_u - PECOV(u-hat)) * A'

#create the diagonal matrix for the marker variance, I*sigma^2_m
marker_gen_var_identity = diag(1, length(marker_effects), length(marker_effects))

#A is defined above and in the equation as 1/sum2pq * Z * G^-1
#this is essentially the variance of the marker estimates.
#Note they are not independent, despite assuming the true effects are independent (off-diag != 0)
#hard to tell if they are close to independent, however. Non-independence is data-drive
#This is because they are driven by G and breeding values and BV are driven by the data and G
var_marker_blups = A %*% (G_mat * gen_var - pecov) %*% t(A)

#from our equation agove, PECOV(m-hat) = var(m) - var(m-hat) = I * sigma^2_m - A * (G*sigma^2_u - PECOV(u-hat)) * A'
marker_pecov = marker_gen_var * marker_gen_var_identity - var_marker_blups
#note again marker estimate error variance (sq root gives standard error) are not independent




####Example to compute haplotype effect and test its effect#####
#example to test haploblock localGEBV effect of first 3 markers
test_mean = marker_effects[1:3]
test_var = marker_pecov[1:3,1:3]

#define the haploblock dosages for localGEBV, this is just a linear contrast of the first 3 markers
#well call this vector c
haplo = c(1,0,2)

#compute haplotype estimated effect, h-hat
haplo_effect = (haplo %*% test_mean)[[1,1]]
#compute the var(h-hat) = var(c * m-hat[1:3]) = c * m_hat[1:3] * c'
#note the reverse transpose - must start with a row vector and R does not assume that
haplo_var = (t(haplo) %*% test_var %*% haplo)[[1,1]]

#compute the test statistic as the effect divided by the standard error - standard Z-test from stat 101
test_statistic = haplo_effect / sqrt(haplo_var)

#compute the p-value: convert the test statistic to the negative side and integrate from -infinity to the test statistic
#2 tailed test (don't care on direction), so we must consider the other side of the normal as well
p_val_freq_test_method = 2 * (1 - pnorm(abs(test_statistic)))

####Example to compute haplo effect and test with Monte Carlo methods####
#sample each marker effect x times utilizing the marker estiamtes and the PECOV
#note if you have a high -log10p-value requirement, then this method could be very bad
#would need millions upon millions of samples potentially
snp_effects = as.data.frame(rmvnorm(100000, test_mean, test_var))

#compute the haplotype effets for all samples (this is the c vector * effects)
haplo_samples = 1 * snp_effects$V1 + 0 * snp_effects$V2 + 2 * snp_effects$V3

#note mean and var of haplo_samples is ~ the mean and variance computed under the frequentist test
mean(haplo_samples)
var(haplo_samples)
haplo_effect
haplo_var

#pval = sum(abs(haplo_samples) >= abs(haplo_effect)) / length(haplo_samples)
#find the probability of greater or lower than 0 (different from 0) - note that <= and >= are OKAY for continuous distributions
#the inclusiong of the = with > or < does not change the integral as the limit does not change (integrating over a single point)
left_side_freq = sum(haplo_samples >= 0) / length(haplo_samples)
right_side_freq = sum(haplo_samples <= 0) / length(haplo_samples)

#p-value is 2 * whichever one is smaller (that's the side it's biased towards)
p_val_MC_method = 2 * min(left_side_freq, right_side_freq)

#or - perform a frquentist test utilizing the SD and mean
p_val_MC_freq_hybrid = 2 * (1 - pnorm(abs(mean(haplo_samples) / sd(haplo_samples))))

#-log10pvals
minuslog10pvals = c(p_val_freq_test_method,p_val_MC_method, p_val_MC_freq_hybrid)
minuslog10pvals = -1 * log10(minuslog10pvals)


marker_effects = data.frame(
  SNP = markers,
  Effect = marker_effects
)
