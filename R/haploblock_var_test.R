##################################
###### Haploblock GWAS Test ######
##################################

#need to add checks for the geno object

#head function to export
haploblock_var_test = function(haploblock_obj, geno, gen_var, threshold = 0.9){

  if(!is.numeric(gen_var) || is.null(gen_var) || is.na(gen_var)){
    stop("Please provide a valid value for the additive genetic variance. This should be the estimate from GBLUP or computed from the additive marker variance from a marker model.")
  }

  null_block_var = gen_var / nrow(haploblock_obj$Haploblocks)

  progressr::handlers("txtprogressbar")
  progressr::with_progress({
    p = progressr::progressor(steps = nrow(haploblock_obj$Haploblocks))

    #extract block IDs
    haploblocks = haploblock_obj$Haploblocks$Block_ID

    #compute effective marker number for each block
    eff_markers = purrr:::map_vec(haploblocks, function(haploblock){

      #extract markers
      markers = unlist(strsplit(haploblock_obj$Haploblocks[haploblock_obj$Haploblocks$Block_ID == haploblock, 1], ";"))

      block_geno = t(geno[geno[,1] %in% markers,4:ncol(geno)])

      if(length(markers) > 1){

        cor_mat = cor(block_geno) #scale and center - helps with eigen decomposition
        eig = eigen(cor_mat, symmetric = TRUE)$values #extract eigenvalues
        cumvar = cumsum(eig) / sum(eig) #cumulative variance explained by each eigenvalue
        df_val = min(which(cumvar > threshold))

      } else{
        df_val = 1
      }

      p()

      return(df_val)
    })
  })

  haploblock_obj$Haploblocks$Eff_Markers = eff_markers
  haploblock_obj$Haploblocks$Test_Stat = haploblock_obj$Haploblocks$Block_Var / null_block_var * haploblock_obj$Haploblocks$Eff_Markers
  haploblock_obj$Haploblocks$P_Value = stats::pchisq(haploblock_obj$Haploblocks$Test_Stat, df = 1, lower.tail = FALSE)
  haploblock_obj$Haploblocks$Bonferroni_P = stats::p.adjust(haploblock_obj$Haploblocks$P_Value, method = "bonferroni", n = length(haploblocks))
  haploblock_obj$Haploblocks$FDR = stats::p.adjust(haploblock_obj$Haploblocks$P_Value, method = "BH", n = length(haploblocks))

  return(haploblock_obj)
}
