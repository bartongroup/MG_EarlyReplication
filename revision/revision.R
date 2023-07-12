REVISION_DIR <- "revision/data"

F_PEAKS_FOR_TADS <- "auto isolated peaks for TADs.xls"
F_PEAK_GROWTH <- "auto isolated W+G peaks for Marek.xls"
F_VALEY_FILLING <- "valley means + minima.tsv"
F_ADJACENT_PEAK_SIMILARITY <- "peakNeighbourSimilarity.xls"


#########################################

read_peak_growth <- function() {
  f <- file.path(REVISION_DIR, F_PEAK_GROWTH)
  stopifnot(file.exists(f))
  
  read_excel(f) |> 
    clean_names() |> 
    rename(start = start_region, end = end_region) |> 
    mutate(
      chr = str_remove(chromo, "chr") |> factor(levels = CHROMOSOMES),
      peak_id = str_glue("{chr}:{start}-{end}")
    ) |> 
    arrange(chr, start) |> 
    mutate(peak_id = as_factor(peak_id)) |>
    pivot_longer(c(starts_with("w_"), starts_with("g_")), names_pattern = "([gw])_(.+)", names_to = c("type", "time_point")) |> 
    mutate(
      type = case_match(
        type,
        "g" ~ "Gaussian",
        "w" ~ "Wavelet"
      )
    ) |> 
    drop_na() |> 
    separate_wider_delim(cols = time_point, names = c("time_start", "time_end"), delim = "_", cols_remove = FALSE) |> 
    mutate(across(c(time_start, time_end), as.integer)) |> 
    mutate(time_point = fct_reorder(time_point, time_start)) |> 
    select(peak_id, chr, start, end, type, time_point, time_start, time_end, value)
}

plot_sel_peak_growth <- function(d, n_peaks, seed = 42) {
  set.seed(seed)
  sel_peaks <- d |>
    count(peak_id) |>
    filter(n == max(n)) |> 
    sample_n(n_peaks) |> 
    pull(peak_id)
  
  d |> 
    filter(peak_id %in% sel_peaks) |> 
    ggplot(aes(x = time_point, y = value, colour = type, group = type)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_line() +
    geom_point() +
    scale_colour_manual(values = okabe_ito_palette) +
    facet_wrap(~ peak_id)
}

time_lm_grad <- function(d) {
  pgn <- d |> 
    mutate(time_mid = (time_start + time_end) / 2) |> 
    select(peak_id, chr, start, end, type, time_mid, value) |> 
    nest(data = c(time_mid, value))
  
  linfit <- pgn |> 
    mutate(
      lm_fit = map(data, ~lm(value ~ time_mid, data = .)),
      tidied = map(lm_fit, tidy),
    ) |> 
    unnest(tidied) |> 
    filter(term != "(Intercept)") |> 
    select(-c(data, lm_fit))
  
  grad <- pgn |> 
    mutate(gradient = map(data, function(x) {
      n <- nrow(x) - 1
      if(n < 2)
        return(tibble(midpoint = double(0), grad = double(0)))
      tibble(
        midpoint = (x$time_mid[1:n] + x$time_mid[2:(n + 1)]) / 2,
        grad = diff(x$value) / diff(x$time_mid)
      )
    })) |>  
    unnest(gradient) |> 
    select(-data)
  
  grad_test <- grad |> 
    select(midpoint, grad) |> 
    nest(data = grad) |> 
    mutate(
      fit = map(data, ~t.test(.x$grad, mu = 0)),
      tidied = map(fit, tidy)
    ) |> 
    unnest(tidied) |> 
    mutate(p.adj = p.adjust(p.value, method = "BH")) |> 
    select(-c(data, fit))

  list(
    lm_fit = linfit,
    gradient = grad,
    grad_test = grad_test
  )
}

plot_time_lm_grad <- function(g_test, ylab1, ylab2, point_size = 2, cex = 2, ymin = NA_real_, ymax = NA_real_) {
  g1 <- g_test$lm_fit |> 
    arrange(peak_id) |> 
    mutate(id = row_number()) |> 
    ggplot(aes(x = id, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, colour = chr)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_point() +
    geom_errorbar() +
    facet_wrap(~ type, ncol = 1) +
    ylim(ymin, ymax) +
    scale_colour_manual(values = rep(okabe_ito_palette, 4)) +
    labs(x = "Peak", y = ylab1)
  
  pgm <- g_test$gradient |> 
    group_by(type, midpoint) |> 
    summarise(
      n = n(),
      m = mean(grad),
      se = sd(grad) / sqrt(n),
      tc = qt(0.975, df = n - 1),
      ci_lo = m - tc * se,
      ci_up = m + tc * se
    ) |> 
    ungroup()
  
  g2 <- g_test$gradient |> 
    ggplot(aes(x = as_factor(midpoint))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_beeswarm(aes(y = grad), cex = cex, size = point_size, colour = "grey70") +
    geom_point(data = pgm, aes(y = m), size = 2) +
    geom_errorbar(data = pgm, aes(ymin = ci_lo, ymax = ci_up), width = 0.2) +
    geom_hline(yintercept = 0, colour = "orange", alpha = 0.8) +
    ylim(ymin, ymax) +
    facet_wrap(~ type, ncol = 1) +
    labs(x = "Midpoint (min)", y = ylab2)
  
  plot_grid(g1, g2, nrow = 1, labels = c("A", "B"))
}


#################################


read_valley_filling <- function() {
  f <- file.path(REVISION_DIR, F_VALEY_FILLING)
  stopifnot(file.exists(f))
  
  read_tsv(f, show_col_types = FALSE) |> 
    clean_names() |> 
    rename(start = start_mbp, end = end_mbp) |> 
    mutate(
      chr = str_remove(chromo, "chr") |> factor(levels = CHROMOSOMES),
      peak_id = str_glue("{chr}:{start}-{end}")
    ) |> 
    arrange(chr, start) |> 
    mutate(peak_id = as_factor(peak_id)) |>
    pivot_longer(c(starts_with("mean_"), starts_with("min_")), names_pattern = "(mean|min)_(.+)", names_to = c("type", "time_point")) |> 
    separate_wider_delim(cols = time_point, names = c("time_start", "time_end"), delim = "_", cols_remove = FALSE) |> 
    mutate(across(c(time_start, time_end), as.integer)) |> 
    mutate(time_point = fct_reorder(time_point, time_start)) |> 
    select(peak_id, chr, start, end, type, time_point, time_start, time_end, value)
}


###################################

read_adjacent_peak_similarity <- function() {
  f <- file.path(REVISION_DIR, F_ADJACENT_PEAK_SIMILARITY)
  stopifnot(file.exists(f))
  
  read_excel(f) |> 
    clean_names() |> 
    mutate(peak_id = row_number()) |> 
    mutate(peak_id = str_glue("P{peak_id}") |> as_factor()) |> 
    pivot_longer(-peak_id, names_to = "time_point") |> 
    mutate(time_point = str_remove(time_point, "tm_")) |> 
    separate_wider_delim(cols = time_point, names = c("time_start", "time_end"), delim = "_", cols_remove = FALSE) |> 
    mutate(across(c(time_start, time_end), as.integer)) |> 
    mutate(time_point = fct_reorder(time_point, time_start)) |> 
    add_column(
      type = "Activation",
      chr = "1",
      start = 1,
      end = 1
    ) |> 
    drop_na()
}

time_point_repeated_anova <- function(d) {
  aov(value ~ time_point + Error(peak_id), data = d) |> 
    tidy()
}
