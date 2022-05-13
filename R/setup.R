DATA_DIR <- "dna_seq"

TIMINGS <- data.frame(
  timing = 0:7,
  start = c(0, 10, 40, 70, 100, 0, 0, 0),
  end = c(40, 40, 70, 100, 130, 130, 60, 30)
)

CHROMOSOMES <- c(1:22, "X", "Y")

BAD_SAMPLES <- c("E1:TM.1_0-40")  # failed sample with almost no data

create_metadata <- function(meta_file) {
  read_tsv(meta_file, show_col_types = FALSE) %>% 
    left_join(TIMINGS, by = "timing") %>% 
    mutate(type = factor(type, levels = c("PD", "CT"))) %>% 
    mutate(sample = paste0(experiment, ":", type, "_", treatment, ".", replicate, "_", start, "-", end)) %>% 
    mutate(condition = paste(experiment, treatment, sep = "_")) %>% 
    select(-c(timing, sample_no)) %>% 
    select(experiment, sample, raw_sample, treatment, type, condition, start, end, replicate, sorting) %>% 
    mutate(data_path = DATA_DIR) %>% 
    arrange(treatment, type, end, start, replicate)
}

okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "grey80", "grey30", "black")
