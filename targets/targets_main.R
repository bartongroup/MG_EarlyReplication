targets_main <- function() {
  
  prepare_data <- list(
    tar_target(metadata, create_metadata(file.path(DATA_DIR, "metadata.txt"))),
    tar_target(bedcount, read_bed_counts(metadata)),
    tar_target(normfac, set_names(bedcount$count / 1e6, as.character(bedcount$sample))),
    tar_target(controls, filter(metadata, type == "CT")$sample),
    tar_target(pulldowns, filter(metadata, type == "PD")$sample)
  )
  
  process_bedgraphs <- tar_map(
    values = BEDGRAPHS,
    names = SUFFIX,
    tar_target(bg, process_bg(metadata, WINDOW_SIZE, normfac, bedcount)),
    tar_target(gbg, group_bg(bg, BAD_SAMPLES)),
    tar_target(sav_gbg, save_tab(gbg, GSEL, SUFFIX))
  )
  
  qc <- list(
    tar_target(bowtie_report, report_bowtie(metadata)),
    tar_target(fig_correlation_matrix_controls, plot_correlation_matrix(bg_50k, sel = controls)),
    tar_target(fig_correlation_matrix_pulldowns, plot_correlation_matrix(bg_50k, sel = pulldowns)), 
    tar_target(fig_clustering_controls, plot_clustering(bg_50k, sel = controls, text_size = 8)),
    tar_target(fig_clustering_pulldowns, plot_clustering(bg_50k, sel = pulldowns, text_size = 8)),
    
    tar_target(fig_broken_fit_1, plot_broken_fit(bg_50k, "E2:TM.2_10-40")),
    tar_target(fig_broken_fit_2, plot_broken_fit(bg_50k, "E2:TM.2_70-100")),
    tar_target(fig_broken_fit_3, plot_broken_fit(bg_50k, "E2:TM.2_0-130")),
    tar_target(fig_cis_example_1, plot_normalization(bg_50k, "1", c(40, 60), "E2:TM.2_10-40", ylim = c(0, 22000))),
    tar_target(fig_cis_example_2, plot_normalization(bg_50k, "1", c(40, 60), "E2:TM.2_70-100", ylim = c(0, 22000))),
    tar_target(fig_cis_example_3, plot_normalization(bg_50k, "1", c(40, 60), "E2:TM.2_0-130", ylim = c(0, 22000))),
    
    tar_target(fig_pull_eff, plot_pulldown_efficiency(bedcount, bg_50k))
  )
  
  figures <- list(
    tar_target(fig_ridge_1, plot_ridge_signal(gbg_10k, GSEL, "1", c(30, 34))),
    tar_target(fig_ridge_2, plot_ridge_signal(gbg_10k, GSEL, "21", c(28, 30)))
  )
  

  c(
    prepare_data,
    process_bedgraphs,
    qc,
    figures
  )
  
}
