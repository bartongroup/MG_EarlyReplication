# Read signal file
read_gbg <- function() {
  read_tsv(F_GBG, show_col_types = FALSE)
}

# Read TAD data
read_tads <- function() {
  f <- file.path(REVISION_DIR, F_TADS)
  stopifnot(file.exists(f))
  
  read_tsv(f, show_col_types = FALSE) |> 
    mutate(
      chr = factor(chr, levels = CHROMOSOMES),
      tad_start = start / 1e6,
      tad_end = end / 1e6
    ) |> 
    mutate(tad_id = sprintf("TAD_%04d", row_number()), .before = 1)
}

# Read the Excel file with peaks for TADs
read_peaks_for_tads <- function() {
  f <- file.path(REVISION_DIR, F_PEAKS_FOR_TADS)
  stopifnot(file.exists(f))
  
  read_excel(f) |> 
    clean_names() |> 
    rename(start = region_start, end = region_end, peak_start = peak_start_mbp, peak_end = peak_end_mbp) |> 
    mutate(
      chr = str_remove(chromo, "chr") |> factor(levels = CHROMOSOMES),
      peak_id = str_glue("{chr}:{start}-{end}")
    ) |> 
    arrange(chr, start) |> 
    mutate(peak_id = as_factor(peak_id))
}

plot_peak_tads <- function(pft, tads, gbg, peakid, margin = 1, time_point = "TM_100-130", text_size = 8) {
  pk <- pft |> 
    filter(peak_id == peakid)
  delta <- pk$peak_end - pk$peak_start
  r1 <- pk$peak_start - margin * delta
  r2 <- pk$peak_end + margin * delta
  
  sig <- gbg |> 
    filter(chr == pk$chr & pos > r1 & pos < r2) |> 
    select(pos, signal := !!time_point) |> 
    filter(signal > 0)
  
  mx <- max(sig$signal)
  tad <- tads |> 
    filter(chr == pk$chr & tad_end > r1 & tad_start < r2) |> 
    mutate(
      ymin = -0.05 * mx,
      ymax = -1,
      fill = as.factor(seq_along(tad_start) %% 2)
    )
  
  sig |> 
    ggplot() +
    th +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = text_size)
    ) +
    geom_col(aes(x = pos, y = signal), fill = "orange", colour = "orange") +
    labs(x = str_glue("chr{pk$chr}"), y = NULL) +
    geom_rect(data = tad, aes(xmin = tad_start, xmax = tad_end, ymin = ymin, ymax = ymax, fill = fill)) +
    geom_segment(data = pk, aes(x = peak_start, xend = peak_end, y = 0, yend = 0), linewidth = 1.2) +
    scale_fill_manual(values = okabe_ito_palette[c(2,4)]) +
    guides(fill = "none") +
    coord_cartesian(xlim = c(r1, r2))
}

plot_tads_sel <- function(pft, tads, gbg, n_sel = NULL, margin = 1, time_point = "TM_100-130", seed = 42, ncol = 6) {
  if(is.null(n_sel)) {
    ids <- pft$peak_id
  } else {
    set.seed(seed)
    ids <- pft |>
      sample_n(n_sel) |> 
      arrange(chr, peak_start) |> 
      pull(peak_id)
  }
  pl <- map(ids, ~plot_peak_tads(pft, tads, gbg, .x, margin = margin))
  plot_grid(plotlist = pl, ncol = ncol)
}

# Find overlaps between peaks and tads
get_overlaps <- function(pft, tads) {
  r_peaks <- pft |> 
    mutate(start = 1e6 * peak_start, end = 1e6 * peak_end) |> 
    select(seqname = chr, start, end, peak_id) |> 
    GRanges()
  r_tads <- tads |> 
    select(seqname = chr, start, end, tad_id) |> 
    GRanges()
  
  hits <- findOverlaps(r_peaks, r_tads)
  peak_hits <- r_peaks[queryHits(hits)]
  tad_hits <- r_tads[subjectHits(hits)]
  
  ov <- pintersect(peak_hits, tad_hits)
  ov_frac_peaks <- width(ov) / width(peak_hits)
  ov_frac_tads <- width(ov) / width(tad_hits)

  tb_peaks <- pft |> select(peak_id, chr, peak_start, peak_end)
  tb_tads <- tads |> select(tad_id, tad_start, tad_end)
  
  bind_cols(
    as_tibble(mcols(peak_hits)),
    as_tibble(mcols(tad_hits)),
    tibble(overlap_peaks = ov_frac_peaks, overlap_tads = ov_frac_tads)
  ) |> 
    left_join(tb_peaks, by = "peak_id") |> 
    left_join(tb_tads, by = "tad_id") |> 
    mutate(overlap_score = overlap_peaks * overlap_tads)
}

# Randomly generated regions based on the TAD length distribution
random_tads <- function(tads, max_len = 2e6, rand_excess = 2) {
  # Chromosomal coverage
  chr_cover <- tads |> 
    group_by(chr) |> 
    summarise(
      chr_start = min(start),
      chr_end = max(end),
      chr_len = chr_end - chr_start
    )
  
  # TAD length distribution, limited to max_len
  len_dist <- tads |> 
    mutate(len = end - start) |> 
    filter(len <= max_len) |> 
    pull(len)
  mean_len <- mean(len_dist)
  
  # Generate random TADs
  map(1:nrow(chr_cover), function(i) {
    r <- chr_cover[i, ]
    # Select twice the mean number (ugly trick)
    n <- as.integer(r$chr_len / mean_len * rand_excess)
    tds <- sample(len_dist, n)
    tibble(
      start = r$chr_start + cumsum(c(0, tds)),
      end = c(start[-1], 1e9)
    ) |> 
      filter(start <= r$chr_end) |> 
      add_column(chr = r$chr)
  }) |> 
    list_rbind() |> 
    mutate(
      tad_id = sprintf("TAD_%04d", row_number()), .before = 1,
      tad_start = start / 1e6,
      tad_end = end / 1e6
    )
}


bootstrap_overlaps <- function(pft, tads, seed = 42, nb = 100, ...) {
  # A peak score is the max score across all overlapping TADs
  ov <- get_overlaps(pft, tads)
  scores <- ov |> 
    group_by(peak_id) |> 
    summarise(score = max(overlap_score)) |> 
    arrange(peak_id)
  
  set.seed(seed)
  R <- matrix(nrow = nrow(scores), ncol = nb)
  pb <- progress_bar$new(total = nb)
  for(i in 1:nb) {
    rtads <- random_tads(tads, ...)
    rscore <- get_overlaps(pft, rtads) |> 
      group_by(peak_id) |> 
      summarise(score = max(overlap_score)) |> 
      arrange(peak_id)
    R[, i] <- rscore$score
    pb$tick()
  }
  
  colnames(R) <- 1:nb
  rownames(R) <- scores$peak_id
  R |>
    as_tibble(rownames = "peak_id") |> 
    pivot_longer(-peak_id) |> 
    left_join(scores, by = "peak_id") |> 
    group_by(peak_id) |> 
    summarise(
      score = score[1],
      p = 1 - length(which(score > value)) / n(),
      p_adj = p.adjust(p, method = "BH")
    ) |> 
    right_join(ov, by = "peak_id")
}