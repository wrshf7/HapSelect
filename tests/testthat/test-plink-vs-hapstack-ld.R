test_that("PLINK and HapSelect LD implementations output match reasonably closely", {
  skip_if(
    inherits(tryCatch(HapSelect:::find_plink(), error = function(e) e), "error"),
    "PLINK executable not found"
  )

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
  plink_ld <- plink_pairwise_ld_geno(genotypes)

  expect_equal(
    HapSelect_ld[, c("Chrom", "Locus1", "Locus2", "Name1", "Name2")],
    plink_ld[, c("Chrom", "Locus1", "Locus2", "Name1", "Name2")]
  )
  expect_equal(HapSelect_ld$LD, plink_ld$LD, tolerance = 1e-6)
})


# Tests: R, C++ and PLINK agreement --------------------------------------------
#
# Three implementations of the same statistic: ld_func()'s cor() loop, the
# compiled pairwise_ld_cpp(), and PLINK's --r2. They should return the same r^2
# for any pair all three agree is computable.
#
# PLINK writes its results as text, so its precision is the binding constraint -
# hence 1e-6 here, against the 1e-12 used between the two implementations in
# test-ld.R. Fixtures avoid pairs with fewer than three shared individuals, where
# the three tools deliberately disagree; test-ld.R covers that case directly.

plink_agreement_fixture <- function() {
  set.seed(11)
  n_ind <- 40
  base <- sample(0:2, n_ind, replace = TRUE)
  dosages <- do.call(rbind, lapply(1:12, function(i) {
    v <- if (i %% 3 == 0) sample(0:2, n_ind, replace = TRUE) else base
    v[sample(n_ind, 2)] <- NA
    v
  }))
  colnames(dosages) <- paste0("Ind", seq_len(n_ind))
  cbind(data.frame(SNP = paste0("m", 1:12), Chromosome = 1,
                   Position = (1:12) * 1000, stringsAsFactors = FALSE),
        as.data.frame(dosages))
}

# PLINK reports one row per pair it computed; line them up by marker pair so the
# comparison does not depend on row order or on which pairs each tool emitted
ld_by_pair <- function(ld) {
  keys <- paste(pmin(ld$Name1, ld$Name2), pmax(ld$Name1, ld$Name2))
  stats::setNames(ld$LD, keys)
}


test_that("R, C++ and PLINK agree on r-squared", {
  skip_if(
    inherits(tryCatch(HapSelect:::find_plink(), error = function(e) e), "error"),
    "PLINK executable not found"
  )

  geno <- plink_agreement_fixture()

  r_ld <- suppressWarnings(pairwise_ld(geno, parallelize = FALSE))
  c_ld <- pairwise_ld_c(geno, parallelize = FALSE)
  p_ld <- plink_pairwise_ld_geno(geno)

  r_by <- ld_by_pair(r_ld)
  c_by <- ld_by_pair(c_ld)
  p_by <- ld_by_pair(p_ld)

  shared <- intersect(intersect(names(r_by), names(c_by)), names(p_by))
  expect_gt(length(shared), 30)

  expect_equal(unname(c_by[shared]), unname(r_by[shared]), tolerance = 1e-12)
  expect_equal(unname(p_by[shared]), unname(r_by[shared]), tolerance = 1e-6)
  expect_equal(unname(p_by[shared]), unname(c_by[shared]), tolerance = 1e-6)
})


test_that("a window matches PLINK's --ld-window, which counts the index marker", {
  skip_if(
    inherits(tryCatch(HapSelect:::find_plink(), error = function(e) e), "error"),
    "PLINK executable not found"
  )

  geno <- plink_agreement_fixture()
  window <- 4

  # PLINK counts the index marker in its window, so a window of w markers apart
  # is --ld-window w + 1. Off by one and the comparison silently narrows.
  c_ld <- pairwise_ld_c(geno, parallelize = FALSE, window = window)
  p_ld <- plink_pairwise_ld_geno(geno, ld_window = window + 1L)

  expect_setequal(names(ld_by_pair(c_ld)), names(ld_by_pair(p_ld)))
  expect_true(all(c_ld$Locus2 - c_ld$Locus1 <= window))

  shared <- names(ld_by_pair(c_ld))
  expect_equal(unname(ld_by_pair(p_ld)[shared]), unname(ld_by_pair(c_ld)[shared]),
               tolerance = 1e-6)
})


test_that("an r-squared floor matches PLINK's --ld-window-r2", {
  skip_if(
    inherits(tryCatch(HapSelect:::find_plink(), error = function(e) e), "error"),
    "PLINK executable not found"
  )

  geno <- plink_agreement_fixture()
  floor <- 0.3

  c_ld <- pairwise_ld_c(geno, parallelize = FALSE, min_r2 = floor)
  p_ld <- plink_pairwise_ld_geno(geno, ld_window_r2 = floor)

  expect_true(all(c_ld$LD >= floor))
  expect_setequal(names(ld_by_pair(c_ld)), names(ld_by_pair(p_ld)))
})
