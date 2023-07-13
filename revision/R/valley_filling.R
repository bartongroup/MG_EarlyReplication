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
