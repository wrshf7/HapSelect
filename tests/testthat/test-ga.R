# Fixtures --------------------------------------------------------------------

make_ga_localgebv <- function() {
  matrix(
    c( 0.8, -0.5,  0.3,  0.1, -0.2,
      -0.3,  0.6, -0.1,  0.4,  0.5,
       0.1,  0.2,  0.7, -0.4,  0.3),
    nrow = 5, ncol = 3,
    dimnames = list(
      paste0("ind", 1:5),
      paste0("B",   1:3)
    )
  )
}


# Tests: wrapper equivalence --------------------------------------------------

test_that("ohs_parent_selection and local_gebv_parent_selection produce identical output", {
  lgebv <- make_ga_localgebv()

  set.seed(42)
  result_ohs <- ohs_parent_selection(
    localGEBV  = lgebv,
    n_founders = 3,
    popSize    = 20,
    maxiter    = 50,
    run        = 10
  )

  set.seed(42)
  result_gebv <- local_gebv_parent_selection(
    localGEBV  = lgebv,
    n_founders = 3,
    popSize    = 20,
    maxiter    = 50,
    run        = 10
  )

  expect_equal(result_ohs$One_Solution, result_gebv$One_Solution)
  expect_equal(result_ohs$GA@fitnessValue, result_gebv$GA@fitnessValue)
})
