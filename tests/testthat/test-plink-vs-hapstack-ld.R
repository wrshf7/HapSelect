test_that("PLINK and HapSelect LD implementations output match reasonably closely", {
  # Use internal call_plink to get the path to the PLINK executable, which will throw an error if PLINK is not installed.
  call_plink <- HapSelect:::call_plink

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

  rownames(genotypes) <- genotypes$marker

  dosage_to_calls <- function(x){
    allele1 <- ifelse(is.na(x), "0", ifelse(x == 2L, "G", "A"))
    allele2 <- ifelse(is.na(x), "0", ifelse(x == 0L, "A", "G"))
    c(rbind(allele1, allele2))
  }

  write_plink_text_files <- function(geno, prefix){
    utils::write.table(
      data.frame(CHR = geno[[2]], SNP = geno[[1]], CM = 0, BP = geno[[3]], stringsAsFactors = FALSE),
      file = paste0(prefix, ".map"),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )

    geno_matrix <- as.matrix(geno[, -(1:3)])
    ped <- do.call(rbind, lapply(seq_len(ncol(geno_matrix)), function(i){
      c(paste0("F", i), paste0("I", i), "0", "0", "0", "-9", dosage_to_calls(geno_matrix[, i]))
    }))

    utils::write.table(
      ped,
      file = paste0(prefix, ".ped"),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
  }

  tmp_dir <- tempfile("plink_compare_ld_")
  dir.create(tmp_dir)
  text_prefix <- file.path(tmp_dir, "fixture")
  bed_prefix <- file.path(tmp_dir, "fixture_bin")

  write_plink_text_files(genotypes, text_prefix)

  make_bed_output <- call_plink(
    c("--file", text_prefix, "--make-bed", "--out", bed_prefix),
    stdout = TRUE,
    stderr = TRUE
  )
  make_bed_status <- attr(make_bed_output, "status")
  if(!is.null(make_bed_status) && make_bed_status != 0){
    fail(paste("PLINK --make-bed failed:\n", paste(make_bed_output, collapse = "\n")))
  }

  HapSelect_ld <- pairwise_ld(genotypes, parallelize = FALSE)
  plink_ld <- plink_pairwise_ld(bed_prefix)

  expect_equal(
    HapSelect_ld[, c("Chrom", "Locus1", "Locus2", "Name1", "Name2")],
    plink_ld[, c("Chrom", "Locus1", "Locus2", "Name1", "Name2")]
  )
  expect_equal(HapSelect_ld$LD, plink_ld$LD, tolerance = 1e-6)
})
