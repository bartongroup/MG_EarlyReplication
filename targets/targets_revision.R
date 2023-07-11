targets_revision <- function(){
  
  peak_growth <- list(
    tar_target(pg, read_peak_growth()),
    tar_target(fig_pg_sel_growth, plot_sel_peak_growth(pg, n_peaks = 16)),
    tar_target(pg_test, growth_test(pg)),
    tar_target(fig_pg_growth, plot_growth(pg_test, "Linear fit slope (kb / min)", "Gradient (kb / min)", ymin = -2, cex = 1, point_size = 0.5))
  )
  
  valley_filling <- list(
    tar_target(vf, read_valley_filling()),
    tar_target(vf_test, growth_test(vf)),
    tar_target(fig_vf_growth, plot_growth(vf_test, "Linear fit slope ( / min)", "Gradient ( / min)", ymin = -0.8, ymax = 1, cex = 0.5, point_size = 0.5))
  )
  
  
  c(
    peak_growth,
    valley_filling
  )
}
  