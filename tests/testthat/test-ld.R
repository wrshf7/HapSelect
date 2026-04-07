test_that("ld_func computes pairwise LD for an in-memory chromosome set", {
  # Synthetic single-chromosome genotype matrix
  # rows are SNP markers, columns 4+ are genotype dosages per individual.
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

  # ld_func uses row names as marker IDs for Name1/Name2 in the output.
  rownames(genotypes) <- genotypes$marker

  observed <- FastStack:::ld_func(genotypes)

  # Enumerate the expected marker pairs for this marker set.
  pair_index <- utils::combn(seq_len(nrow(genotypes)), 2)
  expected_pairs <- data.frame(
    Chrom = genotypes$chrom[pair_index[1, ]],
    Locus1 = pair_index[1, ],
    Locus2 = pair_index[2, ],
    Name1 = genotypes$marker[pair_index[1, ]],
    Name2 = genotypes$marker[pair_index[2, ]],
    stringsAsFactors = FALSE
  )

  # Pull a few known pairs to make the tests easier to read
  m1_m3_ld <- observed$LD[observed$Name1 == "m1" & observed$Name2 == "m3"]
  m1_m4_ld <- observed$LD[observed$Name1 == "m1" & observed$Name2 == "m4"]
  m1_m5_ld <- observed$LD[observed$Name1 == "m1" & observed$Name2 == "m5"]
  m2_m6_ld <- observed$LD[observed$Name1 == "m2" & observed$Name2 == "m6"]
  m4_m8_ld <- observed$LD[observed$Name1 == "m4" & observed$Name2 == "m8"]
  m5_m7_ld <- observed$LD[observed$Name1 == "m5" & observed$Name2 == "m7"]
  m5_m6_ld <- observed$LD[observed$Name1 == "m5" & observed$Name2 == "m6"]
  m7_m8_ld <- observed$LD[observed$Name1 == "m7" & observed$Name2 == "m8"]
  observed_pairs <- paste(observed$Name1, observed$Name2, sep = "::")
  expected_pair_labels <- paste(expected_pairs$Name1, expected_pairs$Name2, sep = "::")

  # Check overall output structure and pair count - choose being the binomial coefficient. E.g., length of combn() output
  expect_equal(nrow(observed), choose(nrow(genotypes), 2))
  expect_named(observed, c("Chrom", "Locus1", "Locus2", "Name1", "Name2", "LD"))

  # Check pair metadata
  expect_true(all(observed$Chrom == 1))
  expect_true(all(observed$Locus1 < observed$Locus2))
  expect_true(all(observed$LD >= 0 & observed$LD <= 1))
  expect_equal(anyDuplicated(observed_pairs), 0)
  expect_equal(observed_pairs, expected_pair_labels)
  expect_equal(observed$Name1, genotypes$marker[observed$Locus1])
  expect_equal(observed$Name2, genotypes$marker[observed$Locus2])

  # Check a selection of known pairwise LD values
  expect_equal(m1_m3_ld, 1)
  expect_equal(m1_m4_ld, 0.4632353, tolerance = 1e-7)
  expect_equal(m1_m5_ld, 0, tolerance = 1e-12)
  expect_equal(m2_m6_ld, 0.8235294, tolerance = 1e-7)
  expect_equal(m4_m8_ld, 0.00147058823529412, tolerance = 1e-12)
  expect_equal(m5_m7_ld, 0.5, tolerance = 1e-12)
  expect_equal(m5_m6_ld, 0.125, tolerance = 1e-12)
  expect_equal(m7_m8_ld, 0.0875, tolerance = 1e-12)
})
