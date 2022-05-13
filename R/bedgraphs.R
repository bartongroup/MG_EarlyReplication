# Read all bedgraph files from the given directory Use metadata to get sample
# names.
# Returns a BG object with all data.
#   pos - GRanges object with window co-ordinates
#   raw - matrix with raw counts
#   score - normalized score (NULL, need to use normalizeBGCounts to fill)
#   normfac - normalizing factors (NULL, need to use normalizeBGCounts to fill)
#   window_size - window size (in bases)
#   sample - vector of sample names
#   condition - vector of condition names
read_bg_files <- function(meta, window_size = 50000) {

  bg_lst <- map(1:nrow(meta), function(i) {
    r <- meta[i, ]
    bg_file <- file.path(r$data_path, "bedgraph", str_glue("{r$raw_sample}.{window_size}.bedgraph"))
    stopifnot(file.exists(bg_file))
    read_tsv(bg_file, col_names = c("chr", "start", "end", all_of(r$sample)), col_types = "ciid") %>% 
      filter(chr %in% CHROMOSOMES)
  })
  
  # genomic positions from the first object
  # we assume they are all the same
  pos <- bg_lst[[1]] %>% select(-4) %>% mutate(pos = (start + end) / 2 / 1e6)
  tab <- bg_lst %>% map_dfc(~select(., 4)) %>% set_names(meta$sample) %>% as.matrix
  
  bg <- list(
    window_size = window_size,
    pos = pos,
    raw = tab,
    score = NULL,
    normfac = NULL,
    metadata = meta
  )
  class(bg) <- append(class(bg), "BG")
  bg
}


clean_bg_sigma <- function(bg, sigma_cutoff = 5) {
  controls <- filter(bg$metadata, type == "CT")$sample
  cutidx <- bg$raw %>% 
    as_tibble() %>% 
    select(all_of(controls)) %>% 
    mutate(id = row_number(), sum = rowSums(.)) %>% 
    select(id, sum) %>% 
    mutate(Z = (sum - mean(sum)) / sd(sum)) %>% 
    filter(Z > sigma_cutoff)
  bg$raw[pull(cutidx, "id"), ] <- 0
  bg
}

normalise_bg_counts <- function(bg, normfac, bedcount) {
  stopifnot(all.equal(names(normfac), as.character(bg$metadata$sample)))
  bg$score <- t(t(bg$raw) / normfac)
  bg$normfac <- normfac
  bg$count <- bedcount
  bg
}

chromosome_structure <- function(bg) {
  bg$pos %>%
    mutate(good = bg$raw %>% rowSums > 0) %>% 
    group_split(chr) %>% 
    map_df(function(tab) {
      tab$id <- with(rle(tab$good), rep(seq_along(lengths), lengths))
      chrom <- tab$chr[1]
      tab %>% 
        group_by(id) %>% 
        filter(good) %>% 
        summarise(n = n(), start = min(pos), end = max(pos)) %>% 
        mutate(length = end - start) %>% 
        filter(length > 1) %>% 
        mutate(chr = chrom) %>% 
        select(chr, start, end)
    })
}


stats_bg <- function(bg) {
  stopifnot(!(is.null(bg$score)))
  means <- list()
  ses <- list()
  for (cond in unique(bg$metadata$condition)) {
    s <- bg$metadata$sample[bg$metadata$condition == cond]
    means[[cond]] <- rowMeans(bg$score[, s, drop = FALSE])
    ses[[cond]] <- apply(bg$score[, s, drop = FALSE], 1, sd) / sqrt(length(s))
  }
  bg$mean.score <- do.call(cbind, means)
  bg$se.score <- do.call(cbind, ses)
  bg$chromosome.structure <- chromosome_structure(bg)
  bg
}


# Create a table with pulldown-control pairing

pull_control_bg <- function(bg) {
  pc <- bg$metadata %>%
    filter(type %in% c("PD", "CT") & sample != "BAD") %>%
    select(experiment, treatment, start, end, replicate, type, sample) %>%
    spread(type, sample) %>%
    mutate(name = gsub("PD_", "", PD))
  bg$pullcont <- pc
  bg
}

remove_short_segments <- function(sel, limit) {
  df <- tibble(
    id = seq_along(sel),
    sel = sel
  )
  r <- df %>%
    group_by(sel, rleid = data.table::rleid(sel)) %>%
    summarise(len = n(), start = min(id), end = max(id))
  r <- r[r$sel & r$len < limit, ]
  # there must be a better way of doing this!
  for (i in 1:nrow(r)) {
    sel[as.integer(r[i,"start"]):as.integer(r[i, "end"])] <- FALSE
  }
  sel
}


normalize_one_cis <- function(dat, min.fragment.length, psi = NULL) {
  dat$T <- dat$C + dat$P
  # group by total count, calculate marginal P/C ratio
  rm <- dat[dat$C > 0 & dat$P > 0, ] %>% group_by(T) %>% summarise(rat = sum(P) / sum(C)) %>% 
    mutate(x = log10(T), y = log10(rat)) %>% as_tibble()
  # fit broken line
  #fit <- lm.br(log10(rat) ~ log10(T), data=rm, type="TL")
  #theta <- as.numeric(fit$coef["theta"])
  fit <- lm(y ~ x, data = rm)
  if (is.null(psi)) psi <- median(rm$x, na.rm = TRUE)
  seg <- segmented::segmented(fit, seg.Z = ~x, psi = psi)
  theta <- seg$psi[2]
  
  # everything below break-point theta is considered "background" with no signal
  bkg.sel <- dat$T < 10^theta & dat$P > 0 & dat$C > 0
  
  # remove background segments shorter than limit
  if (min.fragment.length > 1) {
    bkg.sel <- remove_short_segments(bkg.sel, min.fragment.length)
  }
  
  # tried the ratio of the means (and mean of the ratios), but with new data
  # the distribution of P/C is very skewed in the background
  # cis.norm <- mean(dat[bkg.sel, "P"]) / mean(dat[bkg.sel, "C"])
  #
  # Mode works better
  cis.norm <- get_mode(dat[bkg.sel, ]$P / dat[bkg.sel, ]$C)
  
  diff <- dat$P - cis.norm * dat$C
  bkg <- diff[bkg.sel]
  bkg.sd <- sd(bkg[bkg != 0])
  cis <- list(
    norm = cis.norm,
    theta = theta,
    coefficients = seg$coefficients,
    bkg.sd = bkg.sd
  )
  
  list(
    diff = diff,
    bkg.sel = bkg.sel,
    cis = cis,
    rm = rm,
    seg = seg
  )
}


subtract_background <- function(bg, min.fragment.length = 4, what = "raw", psi = NULL) {
  pulldowns <- bg$pullcont$PD
  controls <- bg$pullcont$CT
  names <- bg$pullcont$name
  
  # calculate normalizations (takes a while)
  R <- map(1:length(pulldowns), function(i) {
    dat <- tibble(
      P = bg[[what]][, pulldowns[i]],
      C = bg[[what]][, controls[i]]
    )
    res <- normalize_one_cis(dat, min.fragment.length, psi)
    res[c("bkg.sel", "diff", "cis")]
  }) %>% set_names(names)
  
  
  bg$normcount <- R %>% map_dfc("diff") %>% as.matrix
  bg$bkg.sel <- R %>% map_dfc("bkg.sel") %>% as.matrix
  bg$cis <- R %>% map("cis")
  bg
}

