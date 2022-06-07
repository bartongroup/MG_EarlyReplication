# Parse one bowtie report file
parse_bowtie_report <- function(file) {
  fs <- readLines(file)
  map_dfr(fs, function(f) {
    x <- unlist(str_split(trimws(f), "\\s+"))
    if (str_detect(f, "reads")) {
      desc <- "total"
    } else if (str_detect(f, "aligned exactly 1 time")) {
      desc <- "unique"
    } else if (str_detect(f, "aligned >1 times")) {
      desc <- "multi"
    } else {
      return(NULL)
    }
    tibble(description = desc, value = x[1])
  })
}


# Parse all bowtie report files, return data frame with results
report_bowtie <- function(meta) {
  map_dfr(1:nrow(meta), function(i) {
    r <- meta[i, ]
    rep_file <- file.path(r$data_path, "bam", str_glue("{r$raw_sample}.report"))
    stopifnot(file.exists(rep_file))
    parse_bowtie_report(rep_file) %>% 
      add_column(sample = r$sample)
  }) %>% 
    mutate(
      value = round(as.integer(value) / 1e6, 3),
      sample = factor(sample, levels = meta$sample)
    ) %>% 
    pivot_wider(id_cols = sample, names_from = description, values_from = value) %>% 
    mutate(
      mapped_perc = round(100 * (multi + unique) / total, 2),
      unique_perc = round(100 * unique / total, 2),
      multi_perc = round(100 * multi / total, 2)
    ) %>% 
    left_join(select(meta, sample, condition), by = "sample")
    set_names(c("Sample", "N total", "N unique", "N multi", "% mapped", "% unique", "% multi", "Condition"))
}

read_bed_counts <- function(meta) {
  map_dfr(1:nrow(meta), function(i) {
    r <- meta[i, ]
    bed_file <- file.path(r$data_path, "bed", str_glue("{r$raw_sample}.count"))
    stopifnot(file.exists(bed_file))
    read_tsv(bed_file, col_names =  "count", show_col_types = FALSE) %>% 
      add_column(sample = r$sample)
  }) %>% 
    mutate(sample = factor(sample, levels = meta$sample))
}  
  

process_bg <- function(meta, window_size, normf, bedcnt) {
  read_bg_files(meta, window_size) %>%
    clean_bg_sigma() %>%
    normalise_bg_counts(normf, bedcnt) %>%
    stats_bg() %>% 
    pull_control_bg() %>% 
    subtract_background(what = "raw")
}

get_mode <- function(x) {
  d <- density(x)
  d$x[d$y == max(d$y)]
}


sample_grouping <- function(bg, bad_samples) {
  # add group integer indices
  mi <- bg$metadata %>%
    group_by(treatment, type, start, end) %>% 
    mutate(group = cur_group_id()) %>% 
    ungroup()
  # create group names
  gn <- bg$metadata %>%
    group_by(treatment, type, start, end) %>%
    summarise() %>% ungroup %>%
    mutate(group_name = paste0(type, "_", treatment, "_", start, "-", end)) %>%
    rownames_to_column(var = "group") %>%
    mutate(group = as.integer(group)) %>%
    select(group,  group_name)
  
  mi <- mi %>% left_join(gn, by = "group") %>%
    select(-group) %>%
    # reorder sample the same way as in original metadata
    mutate(sample = factor(sample, levels = bg$metadata$sample)) %>%
    arrange(sample) %>%
    mutate(sample = as.character(sample))
  
  badSam <- function(bg, bad) {
    p <- bg$pullcont %>% filter(name %in% bad)
    c(p$PD, p$CT)
  }
  
  mi %>%
    mutate(group_name = if_else(sample %in% badSam(bg, bad_samples) | type == "G1", as.character(NA), group_name))
}


# Group data by summing replicate counts according to groups vector.

group_bg <- function(bg, bad_samples) {
  grouping <- sample_grouping(bg, bad_samples)
  groups <- grouping$group_name
  groups[is.na(groups)] = "BAD"
  raw <- bg$raw %>%
    t() %>%
    rowsum(groups) %>%
    t()
  raw <- raw[, -which(colnames(raw) == "BAD")]
  meta <- bg$metadata %>%
    add_column(group = groups) %>%
    filter(group != "BAD") %>% 
    group_by(group) %>%
    summarise(treatment = treatment[1], type = type[1], condition = condition[1], start = start[1], end = end[1]) %>%
    rename(sample = group) %>%
    mutate(sample = factor(sample, levels = colnames(raw)), experiment = "E", replicate = 0) %>%
    mutate(raw_sample = sample)
  count <- bg$count %>% 
    add_column(group = groups) %>%
    filter(group != "BAD") %>% 
    group_by(group) %>%
    summarise(count = sum(count)) %>%
    rename(sample = group) %>%
    mutate(sample = factor(sample, levels = colnames(raw)))
  normfac <- count$count / 1e6
  names(normfac) <- count$sample
  
  bg$raw <- raw
  bg$metadata <- meta
  
  bgg <- bg %>% 
    normalise_bg_counts(normfac, count) %>% 
    stats_bg() %>% 
    pull_control_bg() %>% 
    subtract_background(what = "raw")
  bgg$grouping <- grouping
  
  bgg
}

# extract data set in long format

bg_dat <- function(bg, what = "raw") {
  bg$pos %>% 
    select(chr, pos) %>% 
    bind_cols(as_tibble(bg[[what]])) %>% 
    pivot_longer(-c(chr, pos), names_to = "sample", values_to = what)
}

# extract data set in wide format

bg_tab <- function(bg, what = "raw", samples = NULL) {
  tab <- bg$pos %>% 
    select(chr, pos) %>% 
    bind_cols(as_tibble(bg[[what]]))
  if (!is.null(samples)) {
    tab <- tab[, c("chr", "pos", samples)]
  }
  tab
}

save_tab <- function(bg, samples, out_file_suffix) {
  bg %>%
    bg_tab(what = "normcount", samples = samples) %>%
    write_tsv(str_glue("tab/gbg_{out_file_suffix}.tsv.gz"))
}


remove_zero_rows <- function(tab) {
  zer <- which(rowSums(tab) == 0)
  tab[-zer, ]
}


pulldown_efficiency <- function(bedcnt, bg) {
  bc <- set_names(bedcnt$count, bedcnt$sample)
  md <- bg$pullcont %>%
    mutate(condition = paste(experiment, treatment, sep = "_"))
  map_dfr(bg$pullcont$name, function(pname) {
    pc <- bg$pullcont %>% filter(name == pname)
    pd <- pc$PD
    ct <- pc$CT
    norm <- bg$cis[[pname]]$norm
    pd_cnt <- bc[[pd]]
    ct_cnt <- bc[[ct]] * norm
    pf <- (pd_cnt - ct_cnt) / pd_cnt
    cond <- md %>% filter(name == pname) %>% pull(condition)
    tibble(name = pname, pull_fraction = pf, condition = cond)
  })
}
