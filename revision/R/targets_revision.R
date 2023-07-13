targets_revision <- function(){
  
  tads <- list(
    tar_target(gbg, read_gbg()),
    tar_target(tads, read_tads()),
    tar_target(pft, read_peaks_for_tads()),
    tar_target(fig_tads, plot_tads_sel(pft, tads, gbg, ncol = 6)),
    tar_target(tad_overlaps, get_overlaps(pft, tads)),
    tar_target(tad_boot, bootstrap_overlaps(pft, tads, seed = 42, nb = 1000))
  )
  
  peak_growth <- list(
    tar_target(pg, read_peak_growth()),
    tar_target(n_pg, pg$peak_id |> unique() |> length()),
    tar_target(fig_pg_sel_growth, plot_sel_peak_growth(pg, n_peaks = 16, ylab = "Peak width (kb)")),
    tar_target(pg_lm_grad, time_lm_grad(pg)),
    tar_target(fig_pg_growth, plot_time_lm_grad(pg_lm_grad, "Linear fit slope (kb / min)", "Gradient (kb / min)", ymin = -2, cex = 1, point_size = 0.5))
  )
  
  valley_filling <- list(
    tar_target(vf, read_valley_filling()),
    tar_target(n_vf, vf$peak_id |> unique() |> length()),
    tar_target(fig_vf_sel_filling, plot_sel_peak_growth(vf, n_peaks = 16, ylab = "Valley signal")),
    tar_target(vf_lm_grad, time_lm_grad(vf)),
    tar_target(fig_vf_filling, plot_time_lm_grad(vf_lm_grad, "Linear fit slope (signal / min)", "Gradient (signal / min)", ymin = -0.8, ymax = 1, cex = 0.5, point_size = 0.5))
  )
  
  peak_activation_order <- list(
    tar_target(aps, read_adjacent_peak_similarity()),
    tar_target(n_aps, aps$peak_id |> unique() |> length()),
    tar_target(fig_aps, plot_peak_values(aps, ylab = "Adjacent height similarity")),
    tar_target(aps_test, time_point_repeated_anova(aps))
  )
  
  c(
    tads,
    peak_growth,
    valley_filling,
    peak_activation_order
  )
}
  
