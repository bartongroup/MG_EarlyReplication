plot_correlation_matrix <- function(bg, text.size = 10, sel = NULL) {
  tab <- bg$score
  if (!is.null(sel)) tab <- tab[, sel]
  
  X <- log10(remove_zero_rows(tab) + 0.5)
  corr_mat <- cor(X)
  
  corr_mat %>% 
    as_tibble(rownames = "x") %>% 
    pivot_longer(-x, names_to = "y") %>% 
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_tile() +
    viridis::scale_fill_viridis(option = "cividis") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = text.size),
      axis.text.y = element_text(size = text.size)
    ) +
    labs(x = NULL, y = NULL, fill = "Correlation")
}


plot_clustering <- function(bg, text_size = 10, dist_method = "euclidean",
                            clust_method = "complete", colour_var = "treatment", sel = NULL) {
  tab <- bg$score
  meta <- bg$metadata
  if (!is.null(sel)) {
    tab <- tab[, sel]
    meta <- meta %>% filter(sample %in% sel)
  }

  X <- log10(remove_zero_rows(tab) + 0.5)
  
  dendr <- t(X) %>% 
    dist(method = dist_method) %>% 
    hclust(method = clust_method) %>%
    dendsort::dendsort() %>% 
    ggdendro::dendro_data()
  
  seg <- ggdendro::segment(dendr)
  mt <- meta %>%
    mutate(
      colvar = get(colour_var),
      sample = as.character(sample)
    )
  labs <- ggdendro::label(dendr) %>% 
    as_tibble() %>% 
    left_join(mt, by = c("label" = "sample")) %>% 
    mutate(colour = okabe_ito_palette[as_factor(colvar)]) %>% 
    select(label, colour)
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = text_size, colour = labs$colour),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(size = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  )
  rm(tab, X)
  ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data = seg, aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
    scale_x_continuous(breaks = seq_along(labs$label), labels = labs$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance")
}



plot_broken_fit <- function(bg, pname) {
  pc <- bg$pullcont %>%
    filter(name == pname) %>%
    select(PD, CT)
  dat <- data.frame(
    P = bg$raw[, pc$PD],
    C = bg$raw[, pc$CT]
  )
  res <- normalize_one_cis(dat, min.fragment.length = 4)
  res$rm$yp <- predict(res$seg)
  
  theme_b <- theme_bw() +
    theme(panel.grid = element_blank()) 
  g1 <- ggplot(res$rm) +
    theme_b +
    geom_point(aes(x, y), size = 0.2) +
    geom_line(aes(x, yp), colour = "red") +
    labs(x = "log t", y = "log marginal P/C", title = pname)
  
  bkg <- data.frame(
    P = dat[res$bkg.sel, "P"],
    C = dat[res$bkg.sel, "C"]
  )
  bkg$r <- bkg$P / bkg$C
  
  g2 <- ggplot(bkg) +
    theme_b +
    geom_point(aes(C, P), size = 0.2) + 
    geom_abline(intercept = 0, slope = res$cis$norm, colour = "red")
  
  g3 <- ggplot(bkg) +
    theme_b +
    geom_histogram(aes(x = r, y = ..density..), bins = 100) +
    xlim(0, 2.5) +
    labs(x = "P/C in background") +
    geom_vline(xintercept = res$cis$norm, colour = "red")
  
  cowplot::plot_grid(g1, g2, g3, ncol = 3)
}


# Demonstrate normalization between pulldown and control-background

plot_normalization <- function(bg, chrom, limits, pname, ylim = c(0, NA)) {
  pc <- bg$pullcont %>%
    filter(name == pname) %>%
    select(PD, CT)
  ms <- bind_cols(
    bg$pos,
    as_tibble(bg$raw[, c(pc$PD, pc$CT)]),
    as_tibble(bg$bkg.sel[, pname, drop = FALSE])
  ) %>%
    filter(chr == chrom & pos > limits[1] & pos < limits[2])
  names(ms)[5:7] <- c("P", "C", "B")
  cis <- bg$cis[[pname]]
  ms$CC <- ms$C * cis$norm
  ms$diff <- ms$P - ms$CC
  width <- bg$window_size / 1e6
  ms$pos <- ms$pos - 0.5 * width
  
  ggplot(ms) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_step(aes(pos, P), colour = "gold3") +
    geom_step(aes(pos, CC), colour = "lightblue4") +
    geom_point(data = ms[ms$B, ], aes(pos, CC), colour = "black", size = 0.3) +
    labs(x = "Position (Mb)", y = "Counts per bin", title = pname) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = ylim)
}


plot_pulldown_efficiency <- function(bedcnt, bg) {
  pulldown_efficiency(bedcnt, bg) %>% 
  ggplot() +
    theme_bw() +
    geom_col(aes(x = as_factor(name), y = pull_fraction, fill = as_factor(condition))) +
    labs(x = "Sample", y = "Pulldown fraction", fill = "Condition") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    coord_flip() +
    scale_fill_manual(values = okabe_ito_palette)
}


stepline_xy <- function(x, y) {
  step <- x[2] - x[1]
  n <- length(x)
  idx <- floor(seq(1, n + 0.5, 0.5))
  xx <- c(x[1] - 0.5 * step, x[idx[1:(length(idx) - 1)]] + 0.5 * step)
  yy <- y[idx]
  tibble(x = xx, y = yy)
}


# Plot signal (pulldow count - normalized control count) for all pulldowns
# provided. Creates a compact ridge plot with contrasting colour for easy
# reading.

plot_ridge_signal <- function(bg, samples, chrom, limits, title = NULL, scale = 2,
                            time_levels = c("0-40", "10-40", "40-70", "70-100", "100-130", "0-130", "0-30", "0-60")) {
  mt <- bg$pullcont %>%
    filter(name %in% samples) %>%
    mutate(sample = as.character(name))
  
  df <- bg_dat(bg, "normcount") %>% 
    filter(sample %in% samples & chr == chrom & pos >= limits[1] & pos <= limits[2]) %>% 
    group_split(sample) %>%
    map_df(function(tb) {
      value <- tb$normcount / max(tb$normcount)    # normalize to max
      sl <- stepline_xy(tb$pos, value) %>% set_names("pos", "score")
      sl$sample <- unique(tb$sample)
      sl
    })
  
  ms <- df %>%
    left_join(mt, by = "sample") %>%
    unite(timing, c(start, end), sep = "-") %>%
    mutate(timing = factor(timing, levels = time_levels)) %>%
    arrange(timing) %>%
    mutate(sample = factor(sample, levels = unique(sample))) %>%
    mutate(score = replace(score, score < 0, 0))
  
  pal <- rep(okabe_ito_palette, 4)
  xlab <- str_glue("Position along chr {chrom} (Mb)")
  ggplot(ms, aes(x = pos, y = sample, height = score, group = sample, fill = timing)) +
    theme_bw() +
    ggridges::geom_ridgeline(scale = scale) +
    scale_fill_manual(values = pal) +
    labs(x = xlab, y = "Relative count", title = title) +
    theme(legend.position = "none")
}


