
# test the pairwise_ld function for computing calculate LD across all chromos
test_that("pairwise_ld computes pairwise LD across all chromos", {
  # Synthetic multi-chromosome genotype matrix
  genotypes <- data.frame(
    marker = paste0("m", 1:10),
    chrom = c(rep(1, 5), rep(2, 5)),
    pos = seq(100, 1000, by = 100),
    ind1 = c(0, 0, 2, 0, 1, 2, 1, 0, 1, 2),
    ind2 = c(0, 1, 2, 0, 1, 1, 0, 2, 1, 0),
    ind3 = c(1, 1, 1, 0, 1, 0, 2, 1, 0, 1),
    ind4 = c(1, 2, 1, 1, NA, 0, 1, 2, 0, 0),
    ind5 = c(2, 2, 0, 1, 0, 0, 0, 1, 2, 1),
    ind6 = c(2, 1, 0, 2, 2, 1, 2, 1, 0, 0),
    ind7 = c(1, 0, 1, 2, 1, 2, 1, 0, 1, 2),
    stringsAsFactors = FALSE
  )

  observed <- pairwise_ld(genotypes, parallelize = FALSE)

  # Structure
  expect_named(observed, c("Chrom", "Locus1", "Locus2", "Name1", "Name2", "LD"))
  expect_equal(nrow(observed), 20)  # C(5,2) per chromosome * 2 chromosomes

  # Both chromosomes present with correct row counts
  expect_equal(sort(unique(observed$Chrom)), c(1, 2))
  expect_equal(sum(observed$Chrom == 1), 10)
  expect_equal(sum(observed$Chrom == 2), 10)

  # Locus indices must reset to 1 for each chromosome, not continue from the previous
  chr1 <- observed[observed$Chrom == 1, ]
  chr2 <- observed[observed$Chrom == 2, ]

  expect_equal(min(chr1$Locus1), 1)
  expect_equal(max(chr1$Locus2), 5)
  expect_equal(min(chr2$Locus1), 1)  # would be 6 if reset was broken
  expect_equal(max(chr2$Locus2), 5)  # would be 10 if reset was broken

  # Correct markers assigned to each chromosome
  expect_true(all(chr1$Name1 %in% paste0("m", 1:5)))
  expect_true(all(chr1$Name2 %in% paste0("m", 1:5)))
  expect_true(all(chr2$Name1 %in% paste0("m", 6:10)))
  expect_true(all(chr2$Name2 %in% paste0("m", 6:10)))

  # Spot-check known LD values from each chromosome
  expect_equal(observed$LD[observed$Name1 == "m1" & observed$Name2 == "m2"], 0.25,        tolerance = 1e-7)
  expect_equal(observed$LD[observed$Name1 == "m1" & observed$Name2 == "m3"], 1.00,        tolerance = 1e-7)
  expect_equal(observed$LD[observed$Name1 == "m1" & observed$Name2 == "m5"], 0.00,        tolerance = 1e-12)
  expect_equal(observed$LD[observed$Name1 == "m4" & observed$Name2 == "m5"], 0.10344828,  tolerance = 1e-7)
  expect_equal(observed$LD[observed$Name1 == "m6" & observed$Name2 == "m7"], 0.00,        tolerance = 1e-12)
  expect_equal(observed$LD[observed$Name1 == "m7" & observed$Name2 == "m9"], 0.65625,     tolerance = 1e-7)
  expect_equal(observed$LD[observed$Name1 == "m8" & observed$Name2 == "m10"], 0.82352941, tolerance = 1e-7)
})



# Test the ld_func function for computing pairwise LD from a genotype matrix.
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

  observed <- HapSelect:::ld_func(genotypes)

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


# Tests: pairwise_ld against the R reference, pairwise_ld_r ---------------------
#
# pairwise_ld() must return exactly what pairwise_ld_r() returns at the same
# arguments. The compiled routine computes r^2 from sums rather than through
# cor(), so the two agree to machine precision rather than bit for bit - on the
# bundled wheat data the largest difference over 1.2 million pairs is 2.2e-16.
# The tolerances below leave room for that and nothing more.

# Six markers covering every case where cor(use = "pairwise.complete.obs") does
# something other than compute an ordinary correlation:
#   m1, m2  ordinary, identical to each other
#   m3      monomorphic, so r^2 is NA against anything
#   m4, m5  observed in only two individuals, and anti-correlated with each other
#   m6      observed in one individual
ld_degenerate_fixture <- function() {
  geno <- data.frame(SNP = paste0("m", 1:6), Chromosome = 1,
                     Position = (1:6) * 100, stringsAsFactors = FALSE)
  dosages <- rbind(
    c(0, 1, 2, 0, 1, 2),
    c(0, 1, 2, 0, 1, 2),
    c(1, 1, 1, 1, 1, 1),
    c(0, 1, NA, NA, NA, NA),
    c(1, 0, NA, NA, NA, NA),
    c(0, NA, NA, NA, NA, NA)
  )
  colnames(dosages) <- paste0("Ind", 1:6)
  cbind(geno, as.data.frame(dosages))
}

# Ten markers with scattered missing calls and one monomorphic marker, so the
# argument tests run over something less contrived
ld_mixed_fixture <- function() {
  set.seed(3)
  n_ind <- 12
  base <- sample(0:2, n_ind, replace = TRUE)
  dosages <- do.call(rbind, lapply(1:10, function(i) {
    v <- if (i %% 4 == 0) sample(0:2, n_ind, replace = TRUE) else base
    v[c(i %% n_ind + 1, (i * 3) %% n_ind + 1)] <- NA
    v
  }))
  dosages[3, ] <- 1
  colnames(dosages) <- paste0("Ind", seq_len(n_ind))
  cbind(data.frame(SNP = paste0("m", 1:10), Chromosome = 1,
                   Position = (1:10) * 100, stringsAsFactors = FALSE),
        as.data.frame(dosages))
}


test_that("pairwise_ld matches pairwise_ld_r on every degenerate case", {
  geno <- ld_degenerate_fixture()

  r_ld <- suppressWarnings(pairwise_ld_r(geno, parallelize = FALSE))
  c_ld <- pairwise_ld(geno, parallelize = FALSE)

  expect_equal(c_ld, r_ld, tolerance = 1e-12)

  # and the values really are the awkward ones, not an accident of the fixture:
  # two shared individuals always lie on a line, so r^2 is 1 whatever was observed
  pair <- function(d, a, b) d$LD[d$Name1 == a & d$Name2 == b]
  expect_equal(pair(r_ld, "m4", "m5"), 1)
  expect_equal(pair(c_ld, "m4", "m5"), 1)

  # a pair whose r^2 is undefined is left out rather than reported as NA: m3 is
  # monomorphic so it has no variance to correlate, and m6 is observed in one
  # individual
  expect_false(any(is.na(c_ld$LD)))
  expect_length(pair(c_ld, "m1", "m3"), 0)
  expect_length(pair(c_ld, "m1", "m6"), 0)
  expect_false(any(c_ld$Name1 == "m3" | c_ld$Name2 == "m3"))
})


test_that("pairwise_ld matches pairwise_ld_r for every argument", {
  geno <- ld_mixed_fixture()

  settings <- list(
    list(),
    list(window = 1),
    list(window = 2),
    list(min_r2 = 0.5),
    list(min_obs = 3L),
    list(min_obs = 11L),
    list(window = 3, min_r2 = 0.2, min_obs = 3L)
  )

  for (args in settings) {
    r_ld <- suppressWarnings(do.call(pairwise_ld_r, c(list(geno, parallelize = FALSE), args)))
    c_ld <- do.call(pairwise_ld, c(list(geno, parallelize = FALSE), args))
    expect_equal(c_ld, r_ld, tolerance = 1e-12,
                 info = paste(names(args), unlist(args), sep = "=", collapse = ", "))
  }
})


test_that("the defaults give every within-chromosome pair with a defined r-squared", {
  geno <- ld_mixed_fixture()

  ld <- suppressWarnings(pairwise_ld_r(geno, parallelize = FALSE))

  # the fixture has a monomorphic marker, whose pairs have no r^2 to report
  expect_lt(nrow(ld), choose(nrow(geno), 2))
  expect_false(any(is.na(ld$LD)))
  expect_false(any(ld$Name1 == "m3" | ld$Name2 == "m3"))
  expect_equal(colnames(ld), c("Chrom", "Locus1", "Locus2", "Name1", "Name2", "LD"))
})


test_that("window limits the marker distance compared", {
  geno <- ld_mixed_fixture()
  unwindowed <- pairwise_ld(geno, parallelize = FALSE)

  for (w in 1:3) {
    ld <- pairwise_ld(geno, parallelize = FALSE, window = w)

    expect_true(all(ld$Locus2 - ld$Locus1 <= w))
    # exactly the pairs of the unwindowed result that fall inside the window,
    # which is fewer than every combination because the fixture's monomorphic
    # marker has no r^2 to report
    expect_equal(nrow(ld), sum(unwindowed$Locus2 - unwindowed$Locus1 <= w))
  }
})


test_that("min_r2 drops pairs below the floor", {
  geno <- ld_mixed_fixture()

  unfiltered <- pairwise_ld(geno, parallelize = FALSE)
  ld <- pairwise_ld(geno, parallelize = FALSE, min_r2 = 0.5)

  expect_true(all(ld$LD >= 0.5))
  expect_equal(nrow(ld), sum(unfiltered$LD >= 0.5))
})


test_that("min_obs refuses pairs with too few shared individuals", {
  geno <- ld_degenerate_fixture()

  # m4 and m5 share two individuals: reported as r^2 = 1 by default, NA at 3
  default_ld <- pairwise_ld(geno, parallelize = FALSE)
  strict_ld <- pairwise_ld(geno, parallelize = FALSE, min_obs = 3L)

  pair <- function(d, a, b) d$LD[d$Name1 == a & d$Name2 == b]
  expect_equal(pair(default_ld, "m4", "m5"), 1)
  expect_length(pair(strict_ld, "m4", "m5"), 0)

  # the ordinary pair is untouched
  expect_equal(pair(strict_ld, "m1", "m2"), 1)
})


test_that("an undefined r-squared is never reported, by either implementation", {
  geno <- ld_degenerate_fixture()

  r_ld <- suppressWarnings(pairwise_ld_r(geno, parallelize = FALSE))
  c_ld <- pairwise_ld(geno, parallelize = FALSE)

  expect_false(any(is.na(r_ld$LD)))
  expect_false(any(is.na(c_ld$LD)))

  # cor() would have returned NA for these, so the row count is short of every pair
  expect_lt(nrow(c_ld), choose(nrow(geno), 2))
})


test_that("ld_func_c matches ld_func on a single chromosome", {
  geno <- ld_mixed_fixture()

  expect_equal(ld_func_c(geno), suppressWarnings(ld_func(geno)), tolerance = 1e-12)
  expect_equal(ld_func_c(geno, window = 2, min_obs = 3L),
               suppressWarnings(ld_func(geno, window = 2, min_obs = 3L)),
               tolerance = 1e-12)
})


# Tests: advising when a serial run should be parallelised -----------------------
#
# pairwise_ld() defaults to serial because parallelising costs a couple of
# seconds of worker startup that a windowed run never earns back. An unwindowed
# run on dense data does earn it back, several times over, so that case says so.

ld_advice_fixture <- function(markers_per_chr, n_chr = 2, n_ind = 300) {
  set.seed(5)
  n <- markers_per_chr * n_chr
  dosages <- matrix(sample(0:2, n * n_ind, replace = TRUE), nrow = n)
  colnames(dosages) <- paste0("Ind", seq_len(n_ind))
  cbind(
    data.frame(SNP = sprintf("m%05d", seq_len(n)),
               Chromosome = rep(seq_len(n_chr), each = markers_per_chr),
               Position = rep(seq_len(markers_per_chr), n_chr) * 100,
               stringsAsFactors = FALSE),
    as.data.frame(dosages)
  )
}


test_that("a large unwindowed run advises parallelising", {
  geno <- ld_advice_fixture(markers_per_chr = 3000)

  expect_message(advise_parallel_ld(geno, window = NULL),
                 "running serially over every within-chromosome pair")
  expect_message(advise_parallel_ld(geno, window = NULL), "parallelize = TRUE")
  # names the marker count that triggered it, so the advice can be judged
  expect_message(advise_parallel_ld(geno, window = NULL), "3,000 markers")
})


test_that("a windowed run of the same data does not", {
  geno <- ld_advice_fixture(markers_per_chr = 3000)

  # the same markers, but a window makes the work per chromosome linear rather
  # than quadratic, and it never reaches the threshold
  expect_silent(advise_parallel_ld(geno, window = 20))
})


test_that("a small run does not advise either", {
  expect_silent(advise_parallel_ld(ld_advice_fixture(100), window = NULL))
})


test_that("choosing serial explicitly is respected without comment", {
  geno <- ld_advice_fixture(markers_per_chr = 60)

  # the advice only fires when the default was left in place, so someone who has
  # already decided is not told about it on every call
  expect_silent(pairwise_ld(geno, parallelize = FALSE))
})


test_that("a single chromosome is never worth parallelising", {
  # work is split by chromosome, so one chromosome is one task however large
  geno <- ld_advice_fixture(markers_per_chr = 5000, n_chr = 1)

  expect_silent(advise_parallel_ld(geno, window = NULL))
})


test_that("malformed genotypes are left to check_ld_matrix to report", {
  expect_silent(advise_parallel_ld("not a data frame", NULL))
  expect_silent(advise_parallel_ld(data.frame(a = 1, b = 2), NULL))
})
