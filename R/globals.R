# R CMD check performs static analysis and flags bare names used inside dplyr pipes and
# ggplot2 aes() calls as "undefined global variables" because it cannot tell they refer
# to data frame columns rather than missing objects. globalVariables() suppresses those
# false-positive NOTEs. This has no effect at runtime.
utils::globalVariables(c(
  ".", "Chr", "Chr_len", "Chrom", "Chromosome", "Count", "Cum_Pos",
  "Effect", "End", "End_Pos", "Haplotype_Effect", "LD", "Name1", "Name2",
  "PC1", "PC2", "Pos", "Pos1", "Pos2", "SNP", "Scaled_Block_Var",
  "Start", "Start_Pos", "bin", "chr_len", "chr_start", "dist_bp",
  "dist_kb", "fit", "gen", "group", "method", "se", "xmax", "xmin",
  "y", "ymax", "ymin"
))
