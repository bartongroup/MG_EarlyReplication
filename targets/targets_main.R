targets_main <- function() {
  
  prepare_data <- list(
    tar_target(metadata, create_metadata(file.path(DATA_DIR, "metadata.txt"))),
    tar_target(bedcount, read_bed_counts(metadata)),
    tar_target(normfac, set_names(bedcount$count / 1e6, as.character(bedcount$sample))),
    tar_target(controls, filter(metadata, type == "CT")$sample),
    tar_target(pulldowns, filter(metadata, type == "PD")$sample)
  )
  
  load_bedgraphs <- list(
    tar_target(bg50, process_bg(metadata, 50000, normfac, bedcount)),
    tar_target(bg10, process_bg(metadata, 10000, normfac, bedcount))
  )
  
  group_data <- list(
    tar_target(grouping, sample_grouping(bg10, BAD_SAMPLES)),
    tar_target(gbg50, group_bg(bg50, grouping$group_name)),
    tar_target(gbg10, group_bg(bg10, grouping$group_name)),
    tar_target(gsel, c("TM_0-40", "TM_10-40", "TM_40-70", "TM_70-100", "TM_100-130", "TM_0-130"))
  )
  
  qc <- list(
    tar_target(bowtie_report, report_bowtie(metadata)),
    tar_target(fig_correlation_matrix_controls, plot_correlation_matrix(bg50, sel = controls)),
    tar_target(fig_correlation_matrix_pulldowns, plot_correlation_matrix(bg50, sel = pulldowns)), 
    tar_target(fig_clustering_controls, plot_clustering(bg50, sel = controls, text_size = 8)),
    tar_target(fig_clustering_pulldowns, plot_clustering(bg50, sel = pulldowns, text_size = 8)),
    
    tar_target(fig_broken_fit_1, plot_broken_fit(bg50, "E2:TM.2_10-40")),
    tar_target(fig_broken_fit_2, plot_broken_fit(bg50, "E2:TM.2_70-100")),
    tar_target(fig_broken_fit_3, plot_broken_fit(bg50, "E2:TM.2_0-130")),
    tar_target(fig_cis_example_1, plot_normalization(bg50, "1", c(40, 60), "E2:TM.2_10-40", ylim = c(0, 22000))),
    tar_target(fig_cis_example_2, plot_normalization(bg50, "1", c(40, 60), "E2:TM.2_70-100", ylim = c(0, 22000))),
    tar_target(fig_cis_example_3, plot_normalization(bg50, "1", c(40, 60), "E2:TM.2_0-130", ylim = c(0, 22000))),
    
    tar_target(fig_pull_eff, plot_pulldown_efficiency(bedcount, bg50))
  )
  
  figures <- list(
    tar_target(fig_ridge_1, plot_ridge_signal(gbg10, gsel, "1", c(30, 34))),
    tar_target(fig_ridge_2, plot_ridge_signal(gbg10, gsel, "21", c(28, 30)))
  )
  
  save_data <- list(
    tar_target(sav_gbg10, gbg10 %>% bg_tab(what = "normcount", samples = gsel) %>% write_tsv("tab/gbg10.tsv")),
    tar_target(sav_gbg50, gbg50 %>% bg_tab(what = "normcount", samples = gsel) %>% write_tsv("tab/gbg50.tsv"))
  )
  
  c(
    prepare_data,
    load_bedgraphs,
    group_data,
    qc,
    figures,
    save_data
  )
  
}