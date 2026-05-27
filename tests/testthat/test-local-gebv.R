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
    Chrom    = 1L,
    stringsAsFactors = FALSE
  )
}


# Tests: output structure -----------------------------------------------------

test_that("local GEBV computation returns a named list with correct structure", {
  progressr::handlers("void")
  result <- .compute_local_block_effects(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE
  )

  expect_type(result, "list")
  expect_named(result, c("Haploblocks", "Haplotype_ID_Matrix", "Haplotype_Effect_Matrix", "Haplotypes"))

  # Haploblocks data frame
  expect_named(result$Haploblocks, c("Block_ID", "Block", "Chrom", "Num_Uniq_Hap", "Block_Var", "Unique_Haplo_Block_Var"))
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

test_that("local GEBV computation computes correct haplotype effects with center = TRUE", {
  progressr::handlers("void")
  result <- .compute_local_block_effects(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE
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

test_that("local GEBV computation computes correct block variance", {
  progressr::handlers("void")
  result <- .compute_local_block_effects(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE
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

  result <- .compute_local_block_effects(
    geno           = geno,
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE
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


# Tests: marker_pecov (haplo_test) -------------------------------------------

make_gebv_pecov <- function() {
  m <- diag(c(0.1, 0.2, 0.3))
  rownames(m) <- colnames(m) <- c("m1", "m2", "m3")
  m
}

test_that("local GEBV computation adds PEV and p-value columns when marker_pecov is supplied", {
  progressr::handlers("void")
  result <- .compute_local_block_effects(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    marker_pecov   = make_gebv_pecov(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE
  )

  expect_true("Haplotype_PEV"     %in% names(result$Haplotypes))
  expect_true("Haplotype_P_Value" %in% names(result$Haplotypes))

  # PEV is a quadratic form — always >= 0
  expect_true(all(result$Haplotypes$Haplotype_PEV >= 0, na.rm = TRUE))

  # p-values in [0, 1] where defined
  pvals <- result$Haplotypes$Haplotype_P_Value
  expect_true(all(pvals[!is.nan(pvals)] >= 0 & pvals[!is.nan(pvals)] <= 1))
})


test_that("local GEBV computation computes correct PEV values with diagonal marker_pecov", {
  progressr::handlers("void")
  result <- .compute_local_block_effects(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    marker_pecov   = make_gebv_pecov(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE
  )

  hap <- result$Haplotypes

  # B1 (m1, m2), means = (1.0, 1.0), pecov = diag(0.1, 0.2)
  #   PEV = sum(h_centered_i^2 * C_ii)
  #   "0,1" → centered (-1,  0) → 1*0.1 + 0*0.2 = 0.1
  #   "2,0" → centered ( 1, -1) → 1*0.1 + 1*0.2 = 0.3
  #   "1,2" → centered ( 0,  1) → 0*0.1 + 1*0.2 = 0.2
  #   "1,1" → centered ( 0,  0) → 0
  b1 <- hap[hap$Block_ID == "B1", ]
  expect_equal(b1$Haplotype_PEV[b1$Haplotype == "0,1"], 0.1, tolerance = 1e-10)
  expect_equal(b1$Haplotype_PEV[b1$Haplotype == "2,0"], 0.3, tolerance = 1e-10)
  expect_equal(b1$Haplotype_PEV[b1$Haplotype == "1,2"], 0.2, tolerance = 1e-10)
  expect_equal(b1$Haplotype_PEV[b1$Haplotype == "1,1"], 0.0, tolerance = 1e-10)

  # B2 (m3), mean = 1.25, pecov = diag(0.3)
  #   "2" → centered ( 0.75) → 0.5625 * 0.3 = 0.16875
  #   "1" → centered (-0.25) → 0.0625 * 0.3 = 0.01875
  #   "0" → centered (-1.25) → 1.5625 * 0.3 = 0.46875
  b2 <- hap[hap$Block_ID == "B2", ]
  expect_equal(b2$Haplotype_PEV[b2$Haplotype == "2"], 0.16875, tolerance = 1e-10)
  expect_equal(b2$Haplotype_PEV[b2$Haplotype == "1"], 0.01875, tolerance = 1e-10)
  expect_equal(b2$Haplotype_PEV[b2$Haplotype == "0"], 0.46875, tolerance = 1e-10)
})


test_that("local GEBV computation p-values satisfy 2*(1 - pnorm(|effect/sqrt(PEV)|))", {
  progressr::handlers("void")
  result <- .compute_local_block_effects(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    marker_pecov   = make_gebv_pecov(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE
  )

  hap   <- result$Haplotypes
  valid <- hap[!is.na(hap$Haplotype_PEV) & hap$Haplotype_PEV > 0, ]
  expected_pval <- 2 * (1 - pnorm(abs(valid$Haplotype_Effect / sqrt(valid$Haplotype_PEV))))
  expect_equal(valid$Haplotype_P_Value, expected_pval, tolerance = 1e-10)
})


test_that("local GEBV computation does not add PEV columns when marker_pecov is absent", {
  progressr::handlers("void")
  result <- .compute_local_block_effects(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE
  )

  expect_false("Haplotype_PEV"     %in% names(result$Haplotypes))
  expect_false("Haplotype_P_Value" %in% names(result$Haplotypes))
})


test_that("set_missing_NA = FALSE imputes missing genotypes rather than returning NA", {
  progressr::handlers("void")
  geno           <- make_gebv_geno()
  geno$ind1[geno$marker == "m2"] <- NA_integer_

  result <- .compute_local_block_effects(
    geno           = geno,
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    set_missing_NA = FALSE,
    mean_adjust    = TRUE
  )

  # Missing genotype is zeroed after centering rather than propagated as NA:
  # ind1 at B1: m1=0 centered→-1, m2=NA→0 after centering+impute
  # effect = (-1)*0.5 + 0*(-0.3) = -0.5
  expect_false(is.na(result$Haplotype_Effect_Matrix["B1", "ind1"]))
  expect_equal(result$Haplotype_Effect_Matrix["B1", "ind1"], -0.5, tolerance = 1e-10)
})


test_that("parallel local GEBV matches serial results", {
  progressr::handlers("void")
  old_options <- options(
    future.availableCores.custom = 1L,
    future.availableCores.methods = "custom"
  )
  on.exit(options(old_options), add = TRUE)

  serial_result <- .compute_local_block_effects(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    marker_pecov   = make_gebv_pecov(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE,
    parallel       = FALSE
  )

  parallel_result <- .compute_local_block_effects(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    marker_pecov   = make_gebv_pecov(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE,
    parallel       = TRUE,
    chunk_size     = 1
  )

  expect_equal(parallel_result$Haploblocks, serial_result$Haploblocks, tolerance = 1e-10)
  expect_equal(parallel_result$Haplotype_ID_Matrix, serial_result$Haplotype_ID_Matrix)
  expect_equal(parallel_result$Haplotype_Effect_Matrix, serial_result$Haplotype_Effect_Matrix, tolerance = 1e-10)
  expect_equal(parallel_result$Haplotypes, serial_result$Haplotypes, tolerance = 1e-10)
})


# Tests: wrapper equivalence --------------------------------------------------

test_that("compute_local_GEBV and compute_haplotype_effects produce identical output", {
  progressr::handlers("void")
  args <- list(
    geno           = make_gebv_geno(),
    marker_effects = make_gebv_marker_effects(),
    haploblocks_df = make_gebv_haploblocks(),
    marker_pecov   = make_gebv_pecov(),
    set_missing_NA = TRUE,
    mean_adjust    = TRUE
  )

  result_gebv  <- do.call(compute_local_GEBV,       args)
  result_haplo <- do.call(compute_haplotype_effects, args)

  expect_equal(result_gebv, result_haplo)
})
