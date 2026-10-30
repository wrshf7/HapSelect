# Tests: haploblock variance test ----------------------------------------------
#
# haploblock_var_test() compares each block's variance against an equal share of the genetic
# variance, scaled by how many independent markers the block holds.

# Two blocks over four markers, in the layout compute_local_GEBV() returns: the first column
# holds each block's markers, separated by ";"
var_test_fixture <- function() {
  geno <- data.frame(
    SNP = paste0("m", 1:4), Chromosome = 1, Position = 1:4 * 100,
    i1 = c(0, 0, 2, 1), i2 = c(1, 1, 0, 1), i3 = c(2, 2, 1, 0),
    i4 = c(0, 1, 2, 2), i5 = c(2, 1, 0, 0),
    stringsAsFactors = FALSE
  )
  haploblock_obj <- list(Haploblocks = data.frame(
    Markers = c("m1;m2", "m3;m4"), Block_ID = c("b1", "b2"), Block_Var = c(0.4, 0.1),
    stringsAsFactors = FALSE
  ))
  list(geno = geno, haploblock_obj = haploblock_obj)
}

test_that("haploblock_var_test runs with threshold left at its default", {
  f <- var_test_fixture()

  default <- haploblock_var_test(f$haploblock_obj, f$geno, gen_var = 1)
  explicit <- haploblock_var_test(f$haploblock_obj, f$geno, gen_var = 1, threshold = 0.9)

  expect_identical(default, explicit)
  expect_true(all(c("Eff_Markers", "P_Value", "FDR") %in% names(default$Haploblocks)))
})


test_that("haploblock_var_test still rejects an invalid threshold", {
  f <- var_test_fixture()

  for (bad in list(-0.1, 1.5, NA_real_, "0.9", c(0.5, 0.9))) {
    expect_error(haploblock_var_test(f$haploblock_obj, f$geno, gen_var = 1, threshold = bad),
                 "valid threshold")
  }
})


test_that("haploblock_var_test requires the genetic variance", {
  f <- var_test_fixture()

  expect_error(haploblock_var_test(f$haploblock_obj, f$geno), "additive genetic variance")
})
