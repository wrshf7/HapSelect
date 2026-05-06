test_that("PLINK and HapSelect LD implementations output match reasonably closely", {
  genotypes <- data.frame(
    marker = paste0("m", 1:8),
    chrom = rep(1, 8),
    pos = seq(100, 800, by = 100),
    ind1 = c(0, 0, 2, 0, 1, 2, 0, 1),
    ind2 = c(0, 1, 2, 0, 1, 1, 1, 1),
    ind3 = c(1, 1, 1, 0, 1, 0, 2, 0),
    ind4 = c(1, 2, 1, 1, NA, 0, 1, 2),
    ind5 = c(2, 2, 0, 1, 0, 0, 0, 1),
    ind6 = c(2, 1, 0, 2, 2, 1, 2, 1),
    ind7 = c(1, 0, 1, 2, 1, 2, 1, 0),
    stringsAsFactors = FALSE
  )

  HapSelect_ld <- pairwise_ld(genotypes, parallelize = FALSE)
  plink_ld <- plink_pairwise_ld(genotypes)

  expect_equal(
    HapSelect_ld[, c("Chrom", "Locus1", "Locus2", "Name1", "Name2")],
    plink_ld[, c("Chrom", "Locus1", "Locus2", "Name1", "Name2")]
  )
  expect_equal(HapSelect_ld$LD, plink_ld$LD, tolerance = 1e-6)
})
