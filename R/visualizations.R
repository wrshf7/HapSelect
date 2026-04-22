

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
    theme_cowplot() +
    scale_color_manual(values = color_vec) +
    scale_x_continuous(breaks = chr_info$chr_center, labels = chr_info$Chr) +
    labs(x = "Chromosome", y = "Marker Effect") +
    theme(legend.position = "none")

  return(marker_effects_plot)
}



#####Haplotype effects plot######
unique_haplo_effects_plot = function(haplo_obj, colors = c("#A01FF0", "#A7A8AA"), pos_type = c("midpoint", "start")){

  #plotting strategy - midpoint or start of haplotype
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
    theme_cowplot() +
    scale_color_manual(values = color_vec) +
    scale_x_continuous(breaks = chr_info$chr_center, labels = chr_info$Chr) +
    labs(x = "Chromosome", y = "Unique localGEBV") +
    theme(legend.position = "none")

  return(haplo_effects_plot)

}

######Funnel Plot Creation######

block_var_funnel_plot = function(haplo_obj, mean_line = TRUE, scale_colors = c("blue", "purple", "red")){
  haploblocks = haplo_obj$Haploblocks
  haplotypes = haplo_obj$Haplotypes

  #scale the block variance
  haploblocks$Scaled_Block_Var = log10(haploblocks$Block_Var)
  haploblocks$Scaled_Block_Var = (haploblocks$Scaled_Block_Var - min(haploblocks$Scaled_Block_Var, na.rm = TRUE)) /
    (max(haploblocks$Scaled_Block_Var, na.rm = TRUE) - min(haploblocks$Scaled_Block_Var, na.rm = TRUE))


  haplotypes = left_join(haplotypes, haploblocks, by = "Block_ID")

  funnel_plot = ggplot(haplotypes, aes(x = Haplotype_Effect, y = Scaled_Block_Var, color = Haplotype_Effect)) +
    geom_point(alpha = 0.3, size = 2) +
    theme_cowplot() +
    scale_color_gradient2(low = scale_colors[1], mid = scale_colors[2], high = scale_colors[3], midpoint = 0, name = "Effect Size") +
    labs(x = "localGEBV", y = "Scaled Haploblock Variance")

  if(mean_line){
    funnel_plot = funnel_plot +
      geom_hline(yintercept = mean(haploblocks$Scaled_Block_Var, na.rm = TRUE), linetype = "dashed", color = "black", linewidth = 1, alpha = 0.5)
  }

  return(funnel_plot)
}


######Visualize blocks on chromosomes######
plot_haploblocks = function(haploblock_df, block_fill = "#A01FF0", chrom_fill = "grey90",
                            height = 0.30,
                            single_width_bp = NULL){

  #convert to a factor
  haploblock_df$Chrom = as.factor(haploblock_df$Chrom)
  haploblock_df$Start_Pos = haploblock_df$Start_Pos
  haploblock_df$End_Pos = haploblock_df$End_Pos



  #Extract chromosome sizes
  chrom_sizes = haploblock_df %>% group_by(Chrom) %>% summarise(Chr_len = max(End_Pos, na.rm = TRUE), .groups = "drop") %>%
    mutate(y = row_number())



  # add y to df now
  haploblock_df = haploblock_df %>% mutate(y = as.integer(as.character(
    match(as.character(Chrom), as.character(chrom_sizes$Chrom))
  )))

  #split into multi-region blocks vs single marker blocks
  block_regions = haploblock_df[haploblock_df$Start_Pos != haploblock_df$End_Pos, ]
  single_markers = haploblock_df[haploblock_df$Start_Pos == haploblock_df$End_Pos, ]

  # set single marker width if not provided (small fraction of typical chromosome)
  if (is.null(single_width_bp)) {
    # pick small width relative to median chromosome length
    med_len = median(chrom_sizes$Chr_len, na.rm = TRUE)
    single_width_bp = pmin(med_len * 0.00005, 1)  # ~0.05% of median chr length by default
  }


  # ---- compute rectangle coords for plotting ----
  # chromosome backbone rects:
  chrom_sizes = chrom_sizes %>%
    mutate(ymin = y - height, ymax = y + height)

  # multi-block rect coords
  if (nrow(block_regions) > 0) {
    block_regions = block_regions %>%
      mutate(xmin = Start_Pos, xmax = End_Pos, ymin = y - height, ymax = y + height)
  } else {
    block_regions = block_regions %>% mutate(xmin = numeric(0), xmax = numeric(0),
                                            ymin = numeric(0), ymax = numeric(0))
  }

  # single-block (make narrow rects so they look distinct)
  if (nrow(single_markers) > 0) {
    single_markers = single_markers %>%
      mutate(
        xmin = pmax(0, Start_Pos - single_width_bp / 2),
        xmax = Start_Pos + single_width_bp / 2,
        ymin = y - height, ymax = y + height
      )
  } else {
    single_markers = single_markers %>% mutate(xmin = numeric(0), xmax = numeric(0),
                                              ymin = numeric(0), ymax = numeric(0))
  }


  haploblock_plot = ggplot() +
    #chromosome backbone
    theme_cowplot() +

    #interior of the chromosome filled
    geom_rect(data = chrom_sizes,
              aes(xmin = 0, xmax = Chr_len,
                  ymin = ymin, ymax = ymax),
              fill = chrom_fill, color = NA) +

    #multi snp blocks
    # shaded multi-SNP blocks (no border; we will add boundary lines)
    { if (nrow(block_regions) > 0)
      geom_rect(data = block_regions,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = block_fill, color = NA, alpha = 0.9)
      else NULL } +

    # block boundaries (vertical lines at start/end - only start right now, but deprecating)
    # { if (nrow(block_regions) > 0)
    #   geom_segment(data = block_regions,
    #                aes(x = xmin, xend = xmin, y = ymin, yend = ymax),
    #                color = "black", linewidth = 0.1)
    #   else NULL } +


    # single markers as narrow filled rects + black boundaries
    # { if (nrow(single_markers) > 0)
    #   geom_rect(data = single_markers,
    #             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #             fill = block_fill, color = NA, alpha = 0.4)
    #   else NULL } +
    { if (nrow(single_markers) > 0)
      geom_segment(data = single_markers,
                   aes(x = xmin, xend = xmin, y = ymin, yend = ymax),
                   color = "black", linewidth = 0.1, alpha = 0.4)
      else NULL } +

    geom_rect(data = chrom_sizes, aes(xmin = 0, xmax = Chr_len,
                                      ymin = ymin, #- 0.015,
                                      ymax = ymax),# + 0.015),
              fill = NA, color = "black"
    ) +

    scale_x_continuous("Position", expand = c(0,0)) +
    scale_y_continuous(breaks = chrom_sizes$y, labels = chrom_sizes$Chrom) +
    ylab("Chromosome")

}


###### Function to plot marker density ######


plot_marker_density = function(map_df, bin_size = 500000, height = 0.3,
                               chrom_fill = "grey95",
                               col_low = "white", col_mid = "purple", col_high = "red") {

  map_df$Chrom = as.factor(map_df$Chrom)
  map_df = map_df %>% arrange(as.numeric(as.character(Chrom)), Position)

  # Compute chromosome sizes
  chrom_sizes = map_df %>%
    group_by(Chrom) %>%
    summarise(chr_len = max(Position, na.rm = TRUE), .groups = "drop") %>%
    mutate(y = row_number(),
           ymin = y - height,
           ymax = y + height)

  # Compute bins
  density_df = map_df %>%
    mutate(Pos = Position ) %>%
    group_by(Chrom) %>%
    mutate(bin = floor(Pos / bin_size)) %>%
    group_by(Chrom, bin) %>%
    summarise(
      Start = min(Pos),
      End = max(Pos) + bin_size,
      Count = n(),
      .groups = "drop"
    ) %>%
    left_join(chrom_sizes %>% dplyr::select(Chrom, y, ymin, ymax), by = "Chrom")

  # Compute color scaling limits
  max_count = max(density_df$Count, na.rm = TRUE)
  mid_count = mean(density_df$Count, na.rm = TRUE)

  # ---- Plot ----
  ggplot() +
    theme_cowplot() +

    # Draw density first (behind)
    geom_rect(
      data = density_df,
      aes(xmin = Start, xmax = End, ymin = ymin, ymax = ymax, fill = Count),
      color = NA
    ) +

    # Chromosome outlines on top
    geom_rect(
      data = chrom_sizes,
      aes(xmin = 0, xmax = chr_len, ymin = ymin, ymax = ymax),
      fill = NA, color = "black", linewidth = 0.3
    ) +

    scale_fill_gradient2(
      name = "Marker count",
      low = col_low,
      mid = col_mid,
      high = col_high,
      midpoint = mid_count,
      limits = c(0, max_count),
      na.value = "white"
    ) +

    scale_x_continuous("Position", expand = c(0, 0)) +
    scale_y_continuous(breaks = chrom_sizes$y, labels = chrom_sizes$Chrom, expand = c(0, 0.5)) +
    ylab("Chromosome") +
    theme(legend.position = "right")
}

##### LD Decay Plot ######

plot_ld_decay = function(
    map, ld,
    max_kb = 500,
    point_color = "#A7A8AA",
    curve_color = "#A01FF0",
    alpha = 0.2,
    span = 0.3,
    method = c("gam_tp", "gam_cr", "exp", "loess"),
    k = 50
) {

  # Expect these columns from your objects
  if (!all(c("SNP", "Chromosome", "Position") %in% names(map))) {
    stop("Map must contain: SNP, Chromosome, Position")
  }
  if (!all(c("Name1", "Name2", "LD") %in% names(ld))) {
    stop("LD file must contain: Name1, Name2, LD")
  }

  method = match.arg(method)

  # merge Name1 and Name2 to map to get positions
  ld2 = ld %>%
    rename(SNP = Name1) %>%
    left_join(map, by = "SNP") %>%
    rename(Chrom1 = Chromosome, Pos1 = Position, Name1 = SNP) %>%
    rename(SNP = Name2) %>%
    left_join(map, by = "SNP") %>%
    rename(Chrom2 = Chromosome, Pos2 = Position)

  ld2 = ld2[ld2$Chrom1 == ld2$Chrom2, ]

  # compute distance (bp → kb), using absolute difference
  ld2 = ld2 %>%
    mutate(
      dist_bp = abs(Pos2 - Pos1),
      dist_kb = dist_bp / 1000
    ) %>%
    filter(!is.na(dist_kb)) %>%
    filter(dist_kb <= max_kb)

  # plot LD decay
  ld_decay_plot = ggplot(ld2, aes(x = dist_kb, y = LD, group = 1)) +
    geom_point(
      shape = 21,              # open circles
      color = point_color,
      fill = NA,
      alpha = alpha,
      stroke = 0.3
    ) +
    labs(
      x = "Distance (kb)",
      y = expression(LD *","~R^2),
      title = "LD Decay Curve"
    ) + theme_cowplot()

  if(method %in% c("gam_tp", "gam_cr")){
    formula_str = ifelse(method == "gam_tp",
                         "y ~ s(x, k = k)",
                         "y ~ s(x, bs = 'cr', k = k)"
                         )

    ld_decay_plot = ld_decay_plot + geom_smooth(
      method = "gam",
      formula = as.formula(formula_str),
      se = FALSE,
      color = curve_color,
      linewidth = 1.1
    )
  } else if(method == "exp"){
    nls_fit = nls(LD ~ a * exp(-b * dist_kb),
                   data = ld2,
                   start = list(a = max(ld2$LD), b = 0.01)
    )

    # Generate predicted curve
    grid = data.frame(dist_kb = seq(0, max(ld2$dist_kb), length.out = 500))
    grid$fit = predict(nls_fit, newdata = grid)

    # Plot
    ld_decay_plot = ld_decay_plot +
      geom_line(data = grid, aes(x = dist_kb, y = fit),
                color = curve_color,
                linewidth = 1.1)
  } else{
    ld_decay_plot = ld_decay_plot + geom_smooth(
      method = "loess", span = span, se = FALSE, color = curve_color
    )
  }
  return(ld_decay_plot)
}
