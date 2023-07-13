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

plot_peak_values <- function(d, ylab) {
  d |> 
    ggplot(aes(x = time_point, y = value)) +
    th +
    geom_quasirandom() +
    labs(x = "Time point (min)", y = ylab)
}

time_point_repeated_anova <- function(d) {
  aov(value ~ time_point + Error(peak_id), data = d) |> 
    tidy()
}
