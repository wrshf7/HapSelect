library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())


load("Example_Files/gapit_haploblock_obj.R")
load("Example_Files/gapit_marker_effects.R")
load("Example_Files/gapit_map.R")




######Marker Effects#####
marker_effects_plot = function(marker_effects, chr, pos, colors = c("#A01FF0", "#A7A8AA")){
  # effects: numeric vector of marker effects
  # chr: chromosome IDs (numeric or character)
  # pos: marker positions (numeric)
  # colors: vector of colors to alternate by chromosome
  
  # ensure input lengths match
  if(length(marker_effects) != length(chr) | length(marker_effects) != length(pos) | length(chr) != length(pos)){
    stop("Ensure the length of the marker effects, marker position, and chromosome vectors are the same.")
  }
  
  #group data into a df
  effects_df = data.frame(
    Effect = marker_effects,
    Chr = as.factor(chr),
    Pos = pos
  ) %>% arrange(as.numeric(as.character(Chr)), Pos)
  
  #compute cumulative distance and other chromosome info for the plot
  chr_info = effects_df %>%
    group_by(Chr) %>%
    summarise(chr_len = max(Pos)) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0),
           chr_center = chr_start + chr_len / 2)
  
  effects_df = left_join(effects_df, chr_info, by = "Chr")
  effects_df$Cum_Pos = effects_df$Pos + effects_df$chr_start
  
  #color vector
  color_vec = rep(colors, length.out = length(unique(effects_df$Chr)))
  
  #create ggplot
  marker_effects_plot = ggplot(effects_df, aes(x = Cum_Pos, y = Effect, color = Chr)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = color_vec) +
    scale_x_continuous(breaks = chr_info$chr_center, labels = chr_info$Chr) +
    labs(x = "Chromosome", y = "Marker Effect") + 
    theme(legend.position = "none")
  
  return(marker_effects_plot)
}



#####Haplotype effects plot######
unique_haplo_effects_plot = function(haplo_obj, colors = c("#A01FF0", "#A7A8AA"), pos_type = c("midpoint", "start")){
  
  #plotting strategy - midpoint of start of haplotype
  pos_type = match.arg(pos_type)
  
  haploblocks = haplo_obj$Haploblocks
  haplotypes = haplo_obj$Haplotypes
  
  #get block info for each haplotype
  haplotypes = haplotypes %>%
    left_join(haploblocks, by = "Block_ID") %>%
    mutate(Chr = as.factor(Chrom))
  
  if(pos_type == "midpoint"){
    haplotypes$Pos = ( (haplotypes$Start_Pos + haplotypes$End_Pos)/2 )
  } else { haplotypes$Pos = haplotypes$Start_Pos }
    
  haplotypes = arrange(haplotypes, as.numeric(as.character(Chr)), Pos)
  
  #compute cumulative genomic position
  chr_info = haplotypes %>%
    group_by(Chr) %>%
    summarise(chr_len = max(Pos)) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0),
           chr_center = chr_start + chr_len / 2)
  
  #All in one dataframe
  haplotypes = left_join(haplotypes, chr_info, by = "Chr")
  haplotypes$Cum_Pos = haplotypes$Pos + haplotypes$chr_start
  
  #expand colors to all chromos
  color_vec = rep(colors, length.out = length(unique(haplotypes$Chr)))
  
  #ggplot
  haplo_effects_plot = ggplot(haplotypes, aes(x = Cum_Pos, y = Haplotype_Effect, color = Chr)) +
    geom_point(size = 2, alpha = 0.3) +
    scale_color_manual(values = color_vec) +
    scale_x_continuous(breaks = chr_info$chr_center, labels = chr_info$Chr) +
    labs(x = "Chromosome", y = "Unique Haplotype Effects") +
    theme(legend.position = "none")
  
  return(haplo_effects_plot)
  
}


