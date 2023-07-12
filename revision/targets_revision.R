targets_revision <- function(){
  
  peak_growth <- list(
    tar_target(pg, read_peak_growth()),
    tar_target(fig_pg_sel_growth, plot_sel_peak_growth(pg, n_peaks = 16)),
    tar_target(pg_lm_grad, time_lm_grad(pg)),
    tar_target(fig_pg_growth, plot_time_lm_grad(pg_lm_grad, "Linear fit slope (kb / min)", "Gradient (kb / min)", ymin = -2, cex = 1, point_size = 0.5))
  )
  
  valley_filling <- list(
    tar_target(vf, read_valley_filling()),
    tar_target(vf_lm_grad, time_lm_grad(vf)),
    tar_target(fig_vf_growth, plot_time_lm_grad(vf_lm_grad, "Linear fit slope (signal / min)", "Gradient (signal / min)", ymin = -0.8, ymax = 1, cex = 0.5, point_size = 0.5))
  )
  
  peak_activation_order <- list(
    tar_target(aps, read_adjacent_peak_similarity()),
    tar_target(aps_test, time_point_repeated_anova(aps))
  )
  
  c(
    peak_growth,
    valley_filling,
    peak_activation_order
  )
}
  
