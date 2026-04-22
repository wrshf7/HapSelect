# Fixtures --------------------------------------------------------------------
#
# Three-marker, four-individual, two-block dataset designed so expected
# haplotype effects can be worked out by hand (center = TRUE):
#
#   Genotypes (rows = markers, cols = individuals):
#         ind1  ind2  ind3  ind4
#     m1:   0     2     1     1      mean = 1.0
#     m2:   1     0     2     1      mean = 1.0
#     m3:   2     1     0     2      mean = 1.25
#
#   Marker effects: m1 = 0.5, m2 = -0.3, m3 = 0.8
#
#   Block B1 (m1, m2) — four distinct haplotypes:
#     ind1 "01":  (-1)*0.5 + ( 0)*(-0.3) = -0.5
#     ind2 "20":  ( 1)*0.5 + (-1)*(-0.3) =  0.8
#     ind3 "12":  ( 0)*0.5 + ( 1)*(-0.3) = -0.3
#     ind4 "11":  ( 0)*0.5 + ( 0)*(-0.3) =  0.0
#     Block_Var = mean(0.25, 0.64, 0.09, 0.00) = 0.245
#
#   Block B2 (m3) — three distinct haplotypes:
#     ind1 "2":  ( 0.75)*0.8 =  0.6
#     ind2 "1":  (-0.25)*0.8 = -0.2
#     ind3 "0":  (-1.25)*0.8 = -1.0
#     ind4 "2":  ( 0.75)*0.8 =  0.6
#     Block_Var = mean(0.36, 0.04, 1.00, 0.36) = 0.44

make_gebv_geno <- function() {
  data.frame(
    marker = c("m1", "m2", "m3"),
    chrom  = 1L,
    pos    = c(100L, 200L, 300L),
    ind1   = c(0L, 1L, 2L),
    ind2   = c(2L, 0L, 1L),
    ind3   = c(1L, 2L, 0L),
    ind4   = c(1L, 1L, 2L),
    stringsAsFactors = FALSE
  )
}

make_gebv_marker_effects <- function() {
  data.frame(
    SNP    = c("m1", "m2", "m3"),
    Effect = c(0.5, -0.3, 0.8),
    stringsAsFactors = FALSE
  )
}

make_gebv_haploblocks <- function() {
  data.frame(
    Block_ID = c("B1", "B2"),
    Block    = c("m1;m2", "m3"),
    stringsAsFactors = FALSE
  )
}


# Tests: output structure -----------------------------------------------------

test_that("compute_local_GEBV returns a named list with correct structure", {
  progressr::handlers("void")
  result <- compute_local_GEBV(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = TRUE,
    center         = TRUE
  )

  expect_type(result, "list")
  expect_named(result, c("Haploblocks", "Haplotype_ID_Matrix", "Haplotype_Effect_Matrix", "Haplotypes"))

  # Haploblocks data frame
  expect_named(result$Haploblocks, c("Block_ID", "Block", "Num_Uniq_Hap", "Block_Var", "Unique_Haplo_Block_Var"))
  expect_equal(nrow(result$Haploblocks), 2L)
  expect_equal(result$Haploblocks$Num_Uniq_Hap[result$Haploblocks$Block_ID == "B1"], 4L)
  expect_equal(result$Haploblocks$Num_Uniq_Hap[result$Haploblocks$Block_ID == "B2"], 3L)

  # Haplotype_Effect_Matrix: rows = blocks, cols = individuals
  eff <- result$Haplotype_Effect_Matrix
  expect_equal(nrow(eff), 2L)
  expect_equal(ncol(eff), 4L)
  expect_equal(row.names(eff), c("B1", "B2"))
  expect_equal(colnames(eff),  c("ind1", "ind2", "ind3", "ind4"))

  # Haplotype_ID_Matrix: same shape, values are haplotype ID strings
  ids <- result$Haplotype_ID_Matrix
  expect_equal(nrow(ids), 2L)
  expect_equal(ncol(ids), 4L)
  expect_true(all(grepl("^B[12]:\\d+$", unlist(ids))))

  # Haplotypes lookup table
  expect_named(result$Haplotypes, c("Block_ID", "Haplo_ID", "Haplotype", "Haplotype_Effect"))
  expect_equal(nrow(result$Haplotypes), 7L)  # 4 unique in B1 + 3 unique in B2
})


# Tests: haplotype effect values ----------------------------------------------

test_that("compute_local_GEBV computes correct haplotype effects with center = TRUE", {
  progressr::handlers("void")
  result <- compute_local_GEBV(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = TRUE,
    center         = TRUE
  )

  eff <- result$Haplotype_Effect_Matrix

  expect_equal(eff["B1", "ind1"], -0.5, tolerance = 1e-10)
  expect_equal(eff["B1", "ind2"],  0.8, tolerance = 1e-10)
  expect_equal(eff["B1", "ind3"], -0.3, tolerance = 1e-10)
  expect_equal(eff["B1", "ind4"],  0.0, tolerance = 1e-10)

  expect_equal(eff["B2", "ind1"],  0.6, tolerance = 1e-10)
  expect_equal(eff["B2", "ind2"], -0.2, tolerance = 1e-10)
  expect_equal(eff["B2", "ind3"], -1.0, tolerance = 1e-10)
  expect_equal(eff["B2", "ind4"],  0.6, tolerance = 1e-10)
})


# Tests: block variance -------------------------------------------------------

test_that("compute_local_GEBV computes correct block variance", {
  progressr::handlers("void")
  result <- compute_local_GEBV(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = TRUE,
    center         = TRUE
  )

  blocks <- result$Haploblocks
  expect_equal(blocks$Block_Var[blocks$Block_ID == "B1"], 0.245, tolerance = 1e-10)
  expect_equal(blocks$Block_Var[blocks$Block_ID == "B2"], 0.44,  tolerance = 1e-10)
})


# Tests: missing genotypes ----------------------------------------------------

test_that("set_missing_NA = TRUE produces NA effect for any individual with a missing genotype in the block", {
  progressr::handlers("void")
  geno           <- make_gebv_geno()
  geno$ind1[geno$marker == "m2"] <- NA_integer_

  result <- compute_local_GEBV(
    geno           = geno,
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = TRUE,
    center         = TRUE
  )

  eff <- result$Haplotype_Effect_Matrix

  # ind1 has a missing genotype in B1, so its effect there must be NA
  expect_true(is.na(eff["B1", "ind1"]))

  # ind1 has no missing genotype in B2, so its effect there must not be NA
  expect_false(is.na(eff["B2", "ind1"]))

  # Other individuals at B1 are unaffected
  expect_false(is.na(eff["B1", "ind2"]))
  expect_false(is.na(eff["B1", "ind3"]))
  expect_false(is.na(eff["B1", "ind4"]))
})


test_that("set_missing_NA = FALSE imputes missing genotypes rather than returning NA", {
  progressr::handlers("void")
  geno           <- make_gebv_geno()
  geno$ind1[geno$marker == "m2"] <- NA_integer_

  result <- compute_local_GEBV(
    geno           = geno,
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = FALSE,
    center         = TRUE
  )

  # Missing genotype is zeroed after centering rather than propagated as NA:
  # ind1 at B1: m1=0 centered→-1, m2=NA→0 after centering+impute
  # effect = (-1)*0.5 + 0*(-0.3) = -0.5
  expect_false(is.na(result$Haplotype_Effect_Matrix["B1", "ind1"]))
  expect_equal(result$Haplotype_Effect_Matrix["B1", "ind1"], -0.5, tolerance = 1e-10)
})
