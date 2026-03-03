## Supplementary Figure 1: Sample type preference and detection share by genus
##For the other supplementary figures/tables, please see each related script##
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(patchwork)

# Load pre-filtered data (run filter_data.R first)
data_matched <- read_csv("filtered_data/data_filtered_matched.csv", show_col_types = FALSE)
metadata_matched <- read_csv("data/metadata_matched.csv", show_col_types = FALSE)

# -- Panel A: Sample type preference (weighted mean genome completeness) --

pref_summary <- data_matched %>%
  mutate(genome_completeness = (covered_bases / reference_length) * 100) %>%
  filter(is.finite(genome_completeness), genome_completeness >= 0, genome_completeness <= 100) %>%
  mutate(sample_type = if_else(sample_type == "Indoor Air", "Indoor_Air", sample_type)) %>%
  group_by(genus, sample_type) %>%
  summarise(
    weighted_mean_comp = if_else(
      sum(RPKMF, na.rm = TRUE) > 0,
      weighted.mean(genome_completeness, w = RPKMF, na.rm = TRUE),
      mean(genome_completeness, na.rm = TRUE)
    ),
    n = n(),
    .groups = "drop"
  )

pref_wide <- pref_summary %>%
  pivot_wider(
    names_from = sample_type,
    values_from = c(weighted_mean_comp, n),
    values_fill = list(weighted_mean_comp = 0, n = 0)
  ) %>%
  mutate(
    preference = case_when(
      weighted_mean_comp_Indoor_Air > weighted_mean_comp_Wastewater ~ "Indoor Air",
      weighted_mean_comp_Indoor_Air < weighted_mean_comp_Wastewater ~ "Wastewater",
      TRUE ~ "Tie"
    )
  )

genus_order_levels <- pref_wide %>%
  arrange(desc(weighted_mean_comp_Indoor_Air - weighted_mean_comp_Wastewater)) %>%
  pull(genus)

pref_long <- pref_wide %>%
  select(genus, preference, weighted_mean_comp_Indoor_Air, weighted_mean_comp_Wastewater) %>%
  pivot_longer(
    cols = starts_with("weighted_mean_comp_"),
    names_to = "sample_type",
    values_to = "weighted_mean_comp"
  ) %>%
  mutate(
    sample_type = if_else(sample_type == "weighted_mean_comp_Indoor_Air", "Indoor Air", "Wastewater"),
    genus = factor(genus, levels = genus_order_levels)
  )

genus_hlines <- data.frame(genus = levels(pref_long$genus)) %>%
  mutate(y = seq_along(genus))

pref_plot <- ggplot(pref_long, aes(x = weighted_mean_comp, y = genus, color = preference, group = genus)) +
  geom_hline(data = genus_hlines, aes(yintercept = y), color = "black", linewidth = 0.3, alpha = 0.4) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_point(aes(shape = sample_type), size = 2.8) +
  scale_color_manual(values = c("Indoor Air" = "#377EB8", "Wastewater" = "#8B5A2B", "Tie" = "gray60")) +
  scale_shape_manual(values = c("Indoor Air" = 16, "Wastewater" = 17)) +
  scale_y_discrete(limits = genus_order_levels) +
  labs(
    title = "Sample Type Preference by Genus (Weighted Mean Completeness, RPKMF)",
    subtitle = "Weighted by RPKMF; color shows preferred sample type",
    x = "Weighted Mean Genome Completeness (%)",
    y = NULL,
    color = "Higher in",
    shape = "Sample Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# -- Panel B: Percent of detections by genus --

genus_samples <- data_matched %>%
  count(genus, sample_type, name = "n_detections") %>%
  tidyr::complete(genus, sample_type, fill = list(n_detections = 0)) %>%
  group_by(genus) %>%
  mutate(
    total_detections = sum(n_detections, na.rm = TRUE),
    pct_samples = if_else(total_detections > 0,
                          100 * n_detections / total_detections,
                          0),
    genus = factor(genus, levels = genus_order_levels)
  ) %>%
  mutate(
    pct_preference = case_when(
      pct_samples[sample_type == "Indoor Air"][1] > pct_samples[sample_type == "Wastewater"][1] ~ "Indoor Air",
      pct_samples[sample_type == "Indoor Air"][1] < pct_samples[sample_type == "Wastewater"][1] ~ "Wastewater",
      TRUE ~ "Tie"
    )
  ) %>%
  ungroup()

# -- Panel C: Number of samples detected per genus --

sample_counts <- data_matched %>%
  group_by(genus, sample_type) %>%
  summarise(n_samples = n_distinct(sample_ID), .groups = "drop") %>%
  tidyr::complete(genus, sample_type, fill = list(n_samples = 0)) %>%
  mutate(genus = factor(genus, levels = genus_order_levels))

pct_long <- genus_samples

genus_hlines_pct <- data.frame(genus = levels(pct_long$genus)) %>%
  mutate(y = seq_along(genus))

pct_plot <- ggplot(pct_long, aes(x = pct_samples, y = genus, color = pct_preference, group = genus)) +
  geom_hline(data = genus_hlines_pct, aes(yintercept = y), color = "black", linewidth = 0.3, alpha = 0.4) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_point(aes(shape = sample_type), size = 2.8) +
  scale_color_manual(values = c("Indoor Air" = "#377EB8", "Wastewater" = "#8B5A2B", "Tie" = "gray60")) +
  scale_shape_manual(values = c("Indoor Air" = 16, "Wastewater" = 17)) +
  scale_y_discrete(limits = genus_order_levels, position = "right") +
  labs(
    title = "Percent of Samples Detected by Genus",
    subtitle = "Share of detections per genus (Indoor Air vs Wastewater)",
    x = "% of detections",
    y = NULL,
    color = "Higher in",
    shape = "Sample Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.right = element_text(color = "black"),
    axis.ticks.y.right = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

count_plot <- ggplot(sample_counts,
                     aes(x = sample_type, y = genus, label = n_samples, color = sample_type)) +
  geom_text(fontface = "bold", size = 3.5) +
  scale_color_manual(values = c("Indoor Air" = "#377EB8", "Wastewater" = "#8B5A2B")) +
  scale_y_discrete(limits = genus_order_levels) +
  scale_x_discrete(
    position = "top",
    labels = c("Indoor Air" = "IA", "Wastewater" = "WW")
  ) +
  labs(x = "No. of samples detected", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x.top = element_text(color = c("#377EB8", "#8B5A2B"), face = "bold"),
    axis.title.x.top = element_text(size = 9, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 5.5, 2)
  )

supp_fig1 <- pref_plot | pct_plot | count_plot +
  plot_layout(widths = c(3, 2, 0.1))

print(supp_fig1)

# -- Species detection table (n unique samples per species per sample type) --
# Uncomment to generate and save as CSV (fill in potential_host column manually).
#
# species_table <- data_matched %>%
#   group_by(species, sample_type) %>%
#   summarise(n_samples = n_distinct(sample_ID), .groups = "drop") %>%
#   tidyr::pivot_wider(
#     names_from = sample_type,
#     values_from = n_samples,
#     values_fill = 0
#   ) %>%
#   arrange(desc(rowSums(across(where(is.numeric))))) %>%
#   mutate(potential_host = "")
# 
# # Ensure both columns exist even if a sample type had no detections
# if (!"Indoor Air" %in% names(species_table)) species_table$`Indoor Air` <- 0L
# if (!"Wastewater" %in% names(species_table)) species_table$Wastewater <- 0L
# 
# species_table <- species_table %>%
#   select(species, `Indoor Air`, Wastewater, potential_host)
# 
# cat("\n--- Species Detection Table (n samples detected) ---\n")
# print(species_table, n = Inf)
# 
# write_csv(species_table, "filtered_data/species_detection_table.csv")
# cat("✓ Saved to filtered_data/species_detection_table.csv\n")
