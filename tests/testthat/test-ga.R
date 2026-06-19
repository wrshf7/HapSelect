# Fixtures --------------------------------------------------------------------

# 3-individual, 2-block matrix for localGEBV fitness unit tests.
#
# Values are chosen so scores can be worked out by hand:
#
#          B1   B2
# ind1    1.0  0.0
# ind2    0.0  1.0
# ind3    0.5  0.5
#
# Selecting all 3 individuals, no_selfing pairs: (1,2), (1,3), (2,3)
#   B1 midparents: 0.5, 0.75, 0.25  -> max 0.75
#   B2 midparents: 0.5, 0.25, 0.75  -> max 0.75
#   total = 1.5
#
# With selfing, self-pairs (1,1),(2,2),(3,3) lift the per-block max:
#   B1: (ind1 x ind1) = (1+1)/2 = 1.0  -> total B1 max = 1.0
#   B2: (ind2 x ind2) = (1+1)/2 = 1.0  -> total B2 max = 1.0
#   total = 2.0
make_lgebv_fitness_matrix <- function() {
  matrix(
    c(1.0, 0.0, 0.5,
      0.0, 1.0, 0.5),
    nrow = 3, ncol = 2,
    dimnames = list(paste0("ind", 1:3), paste0("B", 1:2))
  )
}

# 2-individual (2 chromosomes each), 2-block matrix for haplotype fitness unit tests.
#
#              B1   B2
# ind1_1       5.0  0.1
# ind1_2       4.0  0.1
# ind2_1       1.0  3.0
# ind2_2       1.0  4.0
#
# Selecting both individuals, the strategy controls which chromosome pairs are valid:
#
# OHS (cross-individual only: ind1 vs ind2)
#   B1: max(5+1, 5+1, 4+1, 4+1) = 6.0   B2: max(0.1+3, 0.1+4, 0.1+3, 0.1+4) = 4.1
#   total = 10.1
#
# Haploid_OHS (also allows same-individual cross-chromosome pairs)
#   adds (ind1_1, ind1_2) and (ind2_1, ind2_2)
#   B1: max now includes 5+4=9  -> max = 9.0
#   B2: max now includes 3+4=7  -> max = 7.0
#   total = 16.0
#
# OPV (also allows each chromosome paired with itself)
#   adds (ind1_1,ind1_1)=10, (ind1_2,ind1_2)=8, (ind2_1,ind2_1)=2, (ind2_2,ind2_2)=2
#   B1: max = 10.0   B2: adds (ind2_2,ind2_2)=8  -> max = 8.0
#   total = 18.0
make_haplotype_fitness_matrix <- function() {
  matrix(
    c(5.0, 4.0, 1.0, 1.0,
      0.1, 0.1, 3.0, 4.0),
    nrow = 4, ncol = 2,
    dimnames = list(
      c("ind1_1", "ind1_2", "ind2_1", "ind2_2"),
      paste0("B", 1:2)
    )
  )
}

# Minimal haploblock_obj for local_gebv_parent_selection integration tests.
# 5 individuals, 3 blocks.
make_lgebv_obj <- function() {
  set.seed(1)
  m <- matrix(rnorm(5 * 3), nrow = 5, ncol = 3,
              dimnames = list(paste0("ind", 1:5), paste0("B", 1:3)))
  list(Haplotype_Effect_Matrix_GA = as.data.frame(m))
}

# Minimal haploblock_obj for haplotype_parent_selection integration tests.
# 4 individuals x 2 chromosomes each, 3 blocks.
make_haplotype_obj <- function() {
  set.seed(1)
  rnames <- paste0(rep(paste0("ind", 1:4), each = 2), "_", rep(1:2, 4))
  m <- matrix(rnorm(8 * 3), nrow = 8, ncol = 3,
              dimnames = list(rnames, paste0("B", 1:3)))
  list(Haplotype_Effect_Matrix_GA = as.data.frame(m))
}


# Tests: build_haplotype_row_metadata ---------------------------------------------------

test_that("build_haplotype_row_metadata parses haplotype-format row names into individual and chromosome columns", {
  m    <- make_haplotype_fitness_matrix()
  meta <- build_haplotype_row_metadata(m)

  expect_equal(meta$individual,  c("ind1", "ind1", "ind2", "ind2"))
  expect_equal(meta$chromosome,  c("1", "2", "1", "2"))
  expect_equal(meta$row,         1:4)
})


# Tests: validate_strategy ----------------------------------------------------

test_that("validate_strategy rejects an invalid localGEBV strategy", {
  expect_error(validate_strategy("bad_strategy", type = "localGEBV"))
})

test_that("validate_strategy rejects an haplotype strategy name passed to localGEBV", {
  expect_error(validate_strategy("OHS", type = "localGEBV"))
})

test_that("validate_strategy rejects an invalid haplotype strategy", {
  expect_error(validate_strategy("selfing", type = "OHS"))
})

test_that("validate_strategy accepts all valid strategies without error", {
  expect_no_error(validate_strategy("selfing",    type = "localGEBV"))
  expect_no_error(validate_strategy("no_selfing", type = "localGEBV"))
  expect_no_error(validate_strategy("OPV",        type = "Haplotype"))
  expect_no_error(validate_strategy("Haploid_OHS",     type = "Haplotype"))
  expect_no_error(validate_strategy("OHS",        type = "Haplotype"))
})


# Tests: fitness_localGEBV ----------------------------------------------------

test_that("fitness_localGEBV no_selfing scores the best midparent value per block", {
  m  <- make_lgebv_fitness_matrix()
  fn <- fitness_localGEBV(m, maximize = TRUE, strategy = "no_selfing")

  # B1 best pair: (ind1, ind3) -> (1.0+0.5)/2 = 0.75
  # B2 best pair: (ind2, ind3) -> (1.0+0.5)/2 = 0.75
  expect_equal(fn(c(1, 2, 3)), 1.5)
})

test_that("fitness_localGEBV selfing scores higher than no_selfing when self-pairing is optimal", {
  m         <- make_lgebv_fitness_matrix()
  fn_self   <- fitness_localGEBV(m, maximize = TRUE, strategy = "selfing")
  fn_noself <- fitness_localGEBV(m, maximize = TRUE, strategy = "no_selfing")

  # Self-pairs lift each block: (ind1 x ind1) -> B1 = 1.0; (ind2 x ind2) -> B2 = 1.0
  expect_equal(fn_self(c(1, 2, 3)), 2.0)
  expect_gt(fn_self(c(1, 2, 3)), fn_noself(c(1, 2, 3)))
})


# Tests: fitness_OHS ----------------------------------------------------------

test_that("fitness_haplotype OHS only scores cross-individual chromosome pairs", {
  m    <- make_haplotype_fitness_matrix()
  meta <- build_haplotype_row_metadata(m)
  fn   <- fitness_haplotype(m, meta,  maximize = TRUE, strategy = "OHS")

  # B1 best cross-pair: ind1_1 + ind2_2 = 5+1 = 6.0
  # B2 best cross-pair: ind1_2 + ind2_2 = 0.1+4 = 4.1
  expect_equal(fn(c(1, 2)), 10.1)
})

test_that("fitness_haplotype hybrid allows same-individual cross-chromosome pairs", {
  m    <- make_haplotype_fitness_matrix()
  meta <- build_haplotype_row_metadata(m)
  fn   <- fitness_haplotype(m, meta, maximize = TRUE, strategy = "Haploid_OHS")

  # B1 best pair: ind1_1 + ind1_2 = 5+4 = 9.0 (same individual, different chromosomes)
  # B2 best pair: ind2_1 + ind2_2 = 3+4 = 7.0
  expect_equal(fn(c(1, 2)), 16.0)
})

test_that("fitness_haplotype OPV allows a chromosome paired with itself", {
  m    <- make_haplotype_fitness_matrix()
  meta <- build_haplotype_row_metadata(m)
  fn   <- fitness_haplotype(m, meta, maximize = TRUE, strategy = "OPV")

  # B1 best pair: ind1_1 x ind1_1 = 5+5 = 10.0
  # B2 best pair: ind2_2 x ind2_2 = 4+4 = 8.0
  expect_equal(fn(c(1, 2)), 18.0)
})

test_that("fitness_haplotype scores increase from most to least restrictive strategy", {
  m    <- make_haplotype_fitness_matrix()
  meta <- build_haplotype_row_metadata(m)

  no_self    <- fitness_haplotype(m, meta, maximize = TRUE, "OHS")(c(1, 2))
  self_cross <- fitness_haplotype(m, meta, maximize = TRUE, "Haploid_OHS")(c(1, 2))
  self_all   <- fitness_haplotype(m, meta, maximize = TRUE, "OPV")(c(1, 2))

  expect_lt(no_self, self_cross)
  expect_lt(self_cross, self_all)
})


# Tests: fitness_localGEBV vs fitness_haplotype -------------------------------------

test_that("fitness_haplotype scores exactly double fitness_localGEBV when chromosomes duplicate individual values", {
  m_lgebv <- make_lgebv_fitness_matrix()

  # haplotype version: expand each individual into 2 identical chromosomes
  m_haplotype <- m_lgebv[rep(seq_len(nrow(m_lgebv)), each = 2), ]
  rownames(m_haplotype) <- paste0(rep(rownames(m_lgebv), each = 2), "_", rep(1:2, nrow(m_lgebv)))
  meta_haplotype <- build_haplotype_row_metadata(m_haplotype)

  fn_lgebv <- fitness_localGEBV(m_lgebv, maximize = TRUE, strategy = "no_selfing")
  fn_haplotype   <- fitness_haplotype(m_haplotype, meta_haplotype, maximize = TRUE, strategy = "OHS")

  # localGEBV uses (a+b)/2; haplotype uses a+b — so haplotype total = localGEBV total * 2
  expect_equal(fn_haplotype(c(1, 2, 3)), fn_lgebv(c(1, 2, 3)) * 2)
})


# Tests: local_gebv_parent_selection ------------------------------------------

test_that("local_gebv_parent_selection returns correct list structure", {
  obj <- make_lgebv_obj()

  set.seed(42)
  result <- local_gebv_parent_selection(
    haploblock_obj = obj,
    strategy       = "no_selfing",
    n_founders     = 3,
    popSize        = 20,
    maxiter        = 50,
    run            = 10,
    maximize       = TRUE,
    monitor        = FALSE
  )

  expect_named(result, c("GA", "selected_founders"))
  expect_named(result$selected_founders, c("indices", "individuals"))
  expect_equal(nrow(result$selected_founders), 3)
  expect_true(all(result$selected_founders$indices %in% 1:5))
  expect_true(all(result$selected_founders$individuals %in% paste0("ind", 1:5)))
})

test_that("local_gebv_parent_selection runs without error for the selfing strategy", {
  obj <- make_lgebv_obj()

  set.seed(42)
  expect_no_error(
    local_gebv_parent_selection(
      haploblock_obj = obj,
      strategy       = "selfing",
      n_founders     = 3,
      popSize        = 20,
      maxiter        = 50,
      run            = 10,
      maximize       = TRUE,
      monitor        = FALSE
    )
  )
})


# Tests: haplotype_parent_selection -------------------------------------------------

test_that("haplotype_parent_selection returns correct list structure with individual-level indices", {
  obj <- make_haplotype_obj()

  set.seed(42)
  result <- haplotype_parent_selection(
    haploblock_obj = obj,
    strategy       = "OHS",
    n_founders     = 3,
    popSize        = 20,
    maxiter        = 50,
    run            = 10,
    maximize       = TRUE,
    monitor        = FALSE
  )

  expect_named(result, c("GA", "selected_founders"))
  expect_named(result$selected_founders, c("indices", "individuals"))
  expect_equal(nrow(result$selected_founders), 3)
  # Indices refer to unique individuals (1-4), not chromosome rows (1-8)
  expect_true(all(result$selected_founders$indices %in% 1:4))
  expect_true(all(result$selected_founders$individuals %in% paste0("ind", 1:4)))
})

test_that("haplotype_parent_selection runs without error for all three strategies", {
  obj <- make_haplotype_obj()

  for (strat in c("OPV",
                  "Haploid_OHS",
                  "OHS")) {
    set.seed(42)
    expect_no_error(
      haplotype_parent_selection(
        haploblock_obj = obj,
        strategy       = strat,
        n_founders     = 3,
        popSize        = 20,
        maxiter        = 50,
        run            = 10,
        maximize       = TRUE,
        monitor        = FALSE
      )
    )
  }
})
