#### Sample Read Analysis ###
# Analyses and plots:
# 1. Sequencing depth (filtered reads) vs. number of detected species — with
#    per-group linear regression, R², p-value annotation, and Pearson/Spearman
#    correlation statistics printed to console
# 2. Rarefaction curves per sample type (read depth as x-axis) — species and strains
# 3. Species accumulation curves per sample type (samples as x-axis)
#
# Requires: filtered_data/data_filtered_matched.csv (run filter_data.R once first)
# Libraries: tidyverse, ggplot2, vegan, patchwork

library(patchwork)

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(tibble)
library(scales)
library(vegan)

# Load pre-existing filtered data (no reprocessing)
data_matched     <- read_csv("filtered_data/data_filtered_matched.csv", show_col_types = FALSE)
metadata_matched <- read_csv("data/metadata_matched.csv", show_col_types = FALSE)

# Consistent colours
sample_type_colors <- c("Indoor Air" = "#377EB8", "Wastewater" = "#8B5A2B")

# Shared minimal theme
panel_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.border       = element_rect(color = "#2C3E50", fill = NA, linewidth = 0.8),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.ticks         = element_line(color = "black"),
    axis.ticks.length  = unit(3, "pt"),
    plot.title         = element_text(face = "bold", size = 12),
    legend.position    = "bottom",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9)
  )

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 1: Sequencing depth vs. species richness — regression & correlation
# ─────────────────────────────────────────────────────────────────────────────

# Step 1: Group by sample first to get library size and species richness per sample
per_sample_stats <- data_matched %>%
  group_by(sample_ID, sample_type) %>%
  summarise(
    total_filtered_reads = dplyr::first(total_filtered_reads_in_sample),
    n_species = n_distinct(species),
    .groups = "drop"
  )

# ── Sequencing depth summary per sample type ─────────────────────────────────

depth_summary <- per_sample_stats %>%
  group_by(sample_type) %>%
  summarise(
    n             = n(),
    total_reads_sum = sum(total_filtered_reads),
    mean_reads    = mean(total_filtered_reads),
    median_reads  = median(total_filtered_reads),
    sd_reads      = sd(total_filtered_reads),
    min_reads     = min(total_filtered_reads),
    max_reads     = max(total_filtered_reads),
    mean_species  = mean(n_species),
    median_species = median(n_species),
    .groups = "drop"
  )

cat("\n── Sequencing depth per sample type ──\n")
for (i in seq_len(nrow(depth_summary))) {
  r <- depth_summary[i, ]
  cat("\n ", r$sample_type, "(n =", r$n, ")\n")
  cat("    Total reads for group: ", formatC(r$total_reads_sum, format = "f", digits = 0, big.mark = ","), "\n")
  cat("    Reads per sample — mean:",   formatC(r$mean_reads,   format = "f", digits = 0, big.mark = ","),
      " median:", formatC(r$median_reads, format = "f", digits = 0, big.mark = ","),
      " SD:",     formatC(r$sd_reads,    format = "f", digits = 0, big.mark = ","), "\n")
  cat("    Range:",
      formatC(r$min_reads, format = "f", digits = 0, big.mark = ","), "–",
      formatC(r$max_reads, format = "f", digits = 0, big.mark = ","), "\n")
  cat("    Species detected — mean:", round(r$mean_species, 1),
      " median:", r$median_species, "\n")
}

# ── Correlation & regression statistics ──────────────────────────────────────

run_stats <- function(df, label) {
  x <- df$total_filtered_reads
  y <- df$n_species
  pearson  <- cor.test(x, y, method = "pearson")
  spearman <- cor.test(x, y, method = "spearman")
  lm_fit   <- lm(n_species ~ total_filtered_reads, data = df)
  lm_sum   <- summary(lm_fit)
  r2       <- lm_sum$r.squared
  f_p      <- pf(lm_sum$fstatistic[1], lm_sum$fstatistic[2],
                 lm_sum$fstatistic[3], lower.tail = FALSE)

  cat("\n──", label, "──\n")
  cat("  n =", nrow(df), "\n")
  cat("  Pearson  r =", round(pearson$estimate, 3),
      "  p =", signif(pearson$p.value, 3), "\n")
  cat("  Spearman ρ =", round(spearman$estimate, 3),
      "  p =", signif(spearman$p.value, 3), "\n")
  cat("  Linear regression:  R² =", round(r2, 3),
      "  F-p =", signif(f_p, 3), "\n")
  cat("  Slope =", signif(coef(lm_fit)[2], 3),
      "  Intercept =", signif(coef(lm_fit)[1], 3), "\n")

  list(pearson = pearson, spearman = spearman, lm = lm_fit, r2 = r2, f_p = f_p)
}

stats_all <- run_stats(per_sample_stats,                                       "ALL SAMPLES")
stats_air <- run_stats(filter(per_sample_stats, sample_type == "Indoor Air"),  "Indoor Air")
stats_ww  <- run_stats(filter(per_sample_stats, sample_type == "Wastewater"),  "Wastewater")

# ── Annotation labels (R² + p per group) ─────────────────────────────────────

make_label <- function(stats, stype) {
  r2_txt <- formatC(stats$r2, digits = 2, format = "f")
  p_txt  <- if (stats$f_p < 0.001) "p < 0.001" else paste0("p = ", signif(stats$f_p, 2))
  data.frame(sample_type = stype, label = paste0("R\u00b2 = ", r2_txt, "\n", p_txt))
}

annot_df <- bind_rows(
  make_label(stats_air, "Indoor Air"),
  make_label(stats_ww,  "Wastewater")
)

label_pos <- per_sample_stats %>%
  group_by(sample_type) %>%
  summarise(x = min(total_filtered_reads),
            y = max(n_species), .groups = "drop") %>%
  left_join(annot_df, by = "sample_type")

plot1 <- ggplot(
  per_sample_stats,
  aes(x = total_filtered_reads, y = n_species, color = sample_type)
) +
  geom_smooth(aes(group = sample_type, fill = sample_type),
              method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(data = label_pos,
            aes(x = x, y = y, label = label, color = sample_type),
            hjust = 0, vjust = 1, size = 3.2, inherit.aes = FALSE) +
  scale_color_manual(values = sample_type_colors, name = "Sample type") +
  scale_fill_manual(values = sample_type_colors, guide = "none") +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6),
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(labels = comma,
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Sequencing depth vs. detected species per sample",
    x     = "Total filtered reads (Millions)",
    y     = "Number of identified species"
  ) +
  panel_theme

print(plot1)

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 2: Rarefaction curves per sample type (read depth on x-axis)
# ─────────────────────────────────────────────────────────────────────────────

# Function to generate rarefaction data and find saturation point
get_rare_data <- function(df, taxon_level) {
  pooled <- df %>%
    group_by(sample_type, .data[[taxon_level]]) %>%
    summarise(reads = sum(reads_aligned), .groups = "drop") %>%
    pivot_wider(names_from = all_of(taxon_level), values_from = reads, values_fill = 0) %>%
    column_to_rownames("sample_type")

  step_val <- max(1, round(min(rowSums(pooled)) / 1000))
  tidy_rare <- rarecurve(pooled, step = step_val, tidy = TRUE)

  sat <- tidy_rare %>%
    group_by(Site) %>%
    arrange(Sample) %>%
    mutate(slope_per_million = ((Species - lag(Species)) / (Sample - lag(Sample))) * 1e6) %>%
    filter(slope_per_million < 1) %>%
    dplyr::slice(1) %>%
    ungroup()

  list(data = tidy_rare, saturation = sat)
}

species_rare <- get_rare_data(data_matched, "species")
strain_rare  <- get_rare_data(data_matched, "strain")

# Define plotting function
plot_rare <- function(rare_obj, title_txt, y_lab) {
  ggplot(rare_obj$data, aes(x = Sample, y = Species, color = Site)) +
    geom_vline(data = rare_obj$saturation, aes(xintercept = Sample, color = Site),
               linetype = "dotted", linewidth = 0.6, alpha = 0.8) +
    geom_text(data = rare_obj$saturation,
              aes(x = Sample, y = 1, label = paste0(round(Sample/1e6, 1), "M")),
              angle = 90, vjust = -0.5, hjust = 0, size = 2.8, fontface = "italic",
              show.legend = FALSE) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = sample_type_colors, name = "Sample type") +
    scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
    labs(title = title_txt, x = "Number of reads (Millions)", y = y_lab) +
    panel_theme
}

plot2a <- plot_rare(species_rare, "Rarefaction curves: Species", "Number of species")
plot2b <- plot_rare(strain_rare,  "Rarefaction curves: Strains", "Number of strains")

print(plot2a)
print(plot2b)

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 3: Species accumulation curves per sample type (samples on x-axis)
# ─────────────────────────────────────────────────────────────────────────────

build_accum <- function(data, taxon_col) {
  bind_rows(lapply(unique(data$sample_type), function(stype) {
    mat <- data %>%
      filter(sample_type == stype) %>%
      group_by(sample_ID, .data[[taxon_col]]) %>%
      summarise(reads = sum(reads_aligned), .groups = "drop") %>%
      pivot_wider(names_from = all_of(taxon_col),
                  values_from = reads, values_fill = 0) %>%
      column_to_rownames("sample_ID")

    acc <- specaccum(mat, method = "random", permutations = 999)

    tibble(sample_type = stype,
           n_samples   = acc$sites,
           richness    = acc$richness,
           sd          = acc$sd)
  }))
}

make_accum_plot <- function(accum_df, title_txt, y_label) {
  # Max-richness endpoint (last point per sample type) for y-axis label
  max_endpoints_df <- accum_df %>%
    group_by(sample_type) %>%
    filter(n_samples == max(n_samples)) %>%
    ungroup()

  plateau_df <- accum_df %>%
    group_by(sample_type) %>%
    arrange(n_samples) %>%
    mutate(delta = richness - lag(richness, default = 0)) %>%
    filter(delta < 1) %>%
    dplyr::slice(1) %>%
    ungroup()

  cat("\nPlateau —", title_txt, ":\n")
  print(plateau_df %>% select(sample_type, n_samples, richness))
  cat("\nMax richness at final sample:\n")
  print(max_endpoints_df %>% select(sample_type, n_samples, richness))

  ggplot(accum_df, aes(x = n_samples, y = richness,
                        color = sample_type, fill = sample_type)) +
    geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_vline(data = plateau_df,
               aes(xintercept = n_samples, color = sample_type),
               linetype = "dashed", linewidth = 0.7, alpha = 0.75) +
    geom_hline(data = plateau_df,
               aes(yintercept = richness, color = sample_type),
               linetype = "dashed", linewidth = 0.7, alpha = 0.75) +
    # Label max species count at the right end of each curve
    geom_text(data = max_endpoints_df,
              aes(x = n_samples, y = richness,
                  label = round(richness),
                  color = sample_type),
              hjust = -0.3, vjust = 0.5, size = 3.2, fontface = "bold",
              inherit.aes = FALSE) +
    scale_color_manual(values = sample_type_colors, name = "Sample type") +
    scale_fill_manual(values = sample_type_colors, name = "Sample type") +
    scale_x_continuous(breaks = scales::pretty_breaks(),
                       expand = expansion(mult = c(0.02, 0.1)),
                       labels = comma) +
    scale_y_continuous(labels = comma) +
    labs(title = title_txt, x = "Number of samples", y = y_label) +
    panel_theme
}

accum_species_df <- build_accum(data_matched, "species")
plot3 <- make_accum_plot(accum_species_df,
                          "Species accumulation curves per sample type",
                          "Cumulative species richness")
print(plot3)

# ─────────────────────────────────────────────────────────────────────────────
# COMBINED FIGURE: 2x2 layout
# ─────────────────────────────────────────────────────────────────────────────


# Arrange 4 plots: A (species accumulation), B (depth vs richness), C & D (rarefaction)
combined_fig <- (plot3 | plot1) / (plot2a | plot2b) +
  plot_annotation(tag_levels = "A", tag_suffix = ")")

print(combined_fig)