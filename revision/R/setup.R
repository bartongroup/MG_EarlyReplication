REVISION_DIR <- "revision/data"

# TADS from Kang et al. 2020
# See http://www.compbio.dundee.ac.uk/user/mgierlinski/earlyrep/doc/analysis.5.html for details
F_TADS <- "U2OS_tads_rep1_tp360_p0.01.tsv"
F_GBG <- "tab/gbg_50k.tsv"

F_PEAKS_FOR_TADS <- "auto isolated peaks for TADs.xls"
F_PEAK_GROWTH <- "auto isolated W+G peaks for Marek.xls"
F_VALEY_FILLING <- "valley means + minima.tsv"
F_ADJACENT_PEAK_SIMILARITY <- "peakNeighbourSimilarity.xls"

th <- ggplot2::theme_bw() +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

