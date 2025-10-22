# Title: "Linkage decay by chromosome"
# Author: Aglaia Szukala

# packages
library(data.table)
library(dplyr)
library(ggplot2)
library(rlang)

# read files by chromosome with LD estimates
setwd("6_Data4LDscripts/")

files <- list.files(pattern = "^ld_.*\\.geno.ld$")
files

ld_list <- lapply(files, function(f) {
  dat <- fread(f)
  dat$CHR_FILE <- sub("^ld_|\\.geno$", "", f)
  dat
})

ld <- rbindlist(ld_list, use.names = TRUE, fill = TRUE)

if ("R^2" %in% names(ld)) setnames(ld, "R^2", "R2")

ld <- ld %>%
  mutate(dist_bp = abs(POS2 - POS1),
         dist_kb = dist_bp / 1000) %>%
  filter(is.finite(R2), dist_bp > 0)

grp_col <- if ("CHR_FILE" %in% names(ld)) "CHR_FILE" else "CHR"

# Bin distances and take median r² (50th percentile) per bin in each chromosome
bin_kb <- 1
ld_binned <- ld %>%
  mutate(bin_kb = floor(dist_kb / bin_kb) * bin_kb) %>%
  group_by(.data[[grp_col]], bin_kb) %>%
  summarise(median_R2 = median(R2, na.rm = TRUE),
            pairs = n(),
            .groups = "drop")

# Add an spline-like via loess (i.e. smooth the median curve)
ld_smooth <- ld_binned %>%
  group_by(.data[[grp_col]]) %>%
  group_modify(~{
    .x <- arrange(.x, bin_kb)
    if (nrow(.x) >= 7) {
      fit <- loess(median_R2 ~ bin_kb, data = .x, span = 0.3)
      .x$median_R2_smooth <- pmin(pmax(predict(fit, .x$bin_kb), 0), 1)
    } else {
      .x$median_R2_smooth <- .x$median_R2
    }
    .x
  }) %>% ungroup()

# Retrieve half decay distance
halfmax_distance <- function(x, y) {
  o <- order(x); x <- x[o]; y <- y[o]
  y <- cummin(y)                   
  hm <- 0.5 * max(y, na.rm = TRUE)
  i  <- which(y <= hm)[1]
  if (is.na(i) || i == 1) return(NA_real_)
  x0 <- x[i-1]; x1 <- x[i]; y0 <- y[i-1]; y1 <- y[i]
  x0 + (hm - y0) * (x1 - x0) / (y1 - y0)
}

# Compute vline positions + label positions
half_df <- ld_smooth %>%
  group_by(.data[[grp_col]]) %>%
  summarise(
    kb_at_halfmax = halfmax_distance(bin_kb, median_R2_smooth),
    y_max         = max(median_R2_smooth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(is.na(kb_at_halfmax), NA_character_,
                   paste0(round(kb_at_halfmax, 1), " kb")),
    # put the label near the top left of the vline
    y_label = y_max * 0.92,
    x_label = kb_at_halfmax
  )

# Plot by chromosome including spline
p <- ggplot(ld_smooth, aes(x = bin_kb, y = median_R2)) +
  geom_point(alpha = 0.6, size = 0.9) +
  geom_line(aes(y = median_R2_smooth), linewidth = 0.7) +
  geom_vline(data = half_df, aes(xintercept = kb_at_halfmax),
             linetype = "dashed") +
  geom_text(data = subset(half_df, !is.na(kb_at_halfmax)),
            aes(x = x_label, y = y_label, label = label),
            hjust = -0.05, vjust = 0, size = 3) +  # nudge a bit to the right
  facet_wrap(reformulate(grp_col), scales = "free_y") +
  labs(
    x = "Distance (kb)",
    y = expression(median~r^2),
    title = "LD decay per chromosome (median r² points, LOESS smooth, half-max marker)"
  ) +
  theme_minimal()

p

# Plot lowest & highest half-decay chromosomes
half_clean <- half_df %>% filter(!is.na(kb_at_halfmax))

# Identify target group values (lowest & highest)
min_id <- half_clean[[grp_col]][ which.min(half_clean$kb_at_halfmax) ]
max_id <- half_clean[[grp_col]][ which.max(half_clean$kb_at_halfmax) ]
targets <- unique(c(min_id, max_id))

# Filter data to two chromosomes of interest
ld_smooth_sel <- ld_smooth %>% filter(.data[[grp_col]] %in% targets)
half_sel      <- half_df      %>% filter(.data[[grp_col]] %in% targets)

extract_chr_num <- function(x) {
  as.integer(sub(".*CHR_0*([0-9]+)_len.*", "\\1", x))
}
lab_map_targets <- setNames(
  paste0("Chr ", extract_chr_num(targets)),
  targets
)

# Plot the chromosome with highest and lowest decay distance
p_extremes <- ggplot(ld_smooth_sel, aes(x = bin_kb, y = median_R2)) +
  geom_point(alpha = 0.6, size = 0.9) +
  geom_line(aes(y = median_R2_smooth), color = "red",  linewidth = 0.8) +
  geom_vline(data = half_sel, aes(xintercept = kb_at_halfmax), linetype = "dashed", color = "blue", linewidth = 0.6) +
  geom_text(data = half_sel,
            aes(x = kb_at_halfmax, y = y_max * 0.92,
                label = paste0(round(kb_at_halfmax, 1), " kb")),
            hjust = -0.05, vjust = 0, size = 4) +
  facet_wrap(reformulate(grp_col), scales = "free_y",
             labeller = as_labeller(lab_map_targets)) +  # <- remove this to keep original names
  labs(
    x = "Distance (kb)",
    y = expression(LD~(r^2))
    #title = "LD decay — chromosomes with lowest vs highest half-decay distance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),             # remove all gridlines
    axis.line = element_line(color = "black"),# keep x and y axis lines
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),           # no gray border around panels
    strip.background = element_blank(),       # clean facet strips
    strip.text = element_text(face = "bold", size= 14)  # bold facet titles
  )


p_extremes

# Get average half LD decay distance across chromosomes
half_summary <- half_df %>%
  filter(!is.na(kb_at_halfmax)) %>%
  summarise(
    mean_half_decay_kb   = mean(kb_at_halfmax, na.rm = TRUE),
    median_half_decay_kb = median(kb_at_halfmax, na.rm = TRUE),
    sd_half_decay_kb     = sd(kb_at_halfmax, na.rm = TRUE),
    min_half_decay_kb    = min(kb_at_halfmax, na.rm = TRUE),
    max_half_decay_kb    = max(kb_at_halfmax, na.rm = TRUE),
    n_chromosomes        = n()
  )

half_summary






