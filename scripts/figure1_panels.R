library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(patchwork)
library(RColorBrewer)

# Load pre-filtered data (run filter_data.R once to create these)
data_matched <- read_csv("filtered_data/data_filtered_matched.csv", show_col_types = FALSE)
metadata_matched <- read_csv("data/metadata_matched.csv", show_col_types = FALSE)

filtered_data_with_type <- data_matched %>%
  mutate(
    completeness = covered_bases / reference_length,
    completeness_bin = cut(completeness, breaks = c(0, 0.25, 0.5, 0.75, 1),
                           labels = c("0-25%", "25-50%", "50-75%", "75-100%"))
  )

panel_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "#2C3E50", fill = NA, linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(3, "pt"),
    plot.title = element_text(face = "bold", size = 12)
  )

matrix_colors <- c("Indoor Air" = "#377EB8", "Wastewater" = "#8B5A2B")

#### Panel A: Species Totals Per Sample Type ###

species_totals <- filtered_data_with_type %>%
  group_by(sample_type) %>%
  summarise(n_species = n_distinct(species), .groups = "drop")

panel_a <- ggplot(species_totals, aes(x = sample_type, y = n_species, fill = sample_type)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n_species), hjust = 1.2, color = "white", fontface = "bold") +
  labs(title = "A. Total species detected", x = NULL, y = "Species count") +
  panel_theme +
  scale_fill_manual(values = matrix_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  coord_flip() +
  theme(legend.position = "none")

#### Panel B: Completeness Distribution ###

completeness_counts <- filtered_data_with_type %>%
  filter(completeness_bin %in% c("50-75%", "75-100%")) %>%
  group_by(sample_ID, sample_type, completeness_bin) %>%
  summarise(n_detections = n(), .groups = "drop")

completeness_medians <- completeness_counts %>%
  group_by(sample_type, completeness_bin) %>%
  summarise(median_detections = median(n_detections, na.rm = TRUE), .groups = "drop")

panel_b <- ggplot(completeness_counts,
                  aes(x = completeness_bin, y = n_detections, fill = sample_type)) +
  geom_violin(trim = FALSE, alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, alpha = 0.7,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = sample_type),
              position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
              alpha = 0.2, size = 1) +
  geom_text(data = completeness_medians,
            aes(x = completeness_bin, y = median_detections, label = round(median_detections, 1),
                group = sample_type),
            position = position_dodge(width = 0.8),
            color = "black", fontface = "bold", vjust = -0.6, inherit.aes = FALSE) +
  labs(title = "D. Detections per sample (>=50% completeness)",
       x = "Completeness bin", y = "Detections per sample") +
  panel_theme +
  scale_fill_manual(values = matrix_colors, name = "Sample type") +
  scale_color_manual(values = matrix_colors, name = "Sample type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06)))

family_sums <- filtered_data_with_type %>%
  group_by(sample_ID, sample_type, match, family) %>%
  summarise(family_reads = sum(reads_aligned, na.rm = TRUE), .groups = "drop")

# Relative abundance from aligned reads.
family_rel_raw <- family_sums %>%
  group_by(sample_ID, sample_type, match) %>%
  mutate(rel_abundance = 100 * family_reads / sum(family_reads),
         rel_abundance = if_else(is.finite(rel_abundance), rel_abundance, 0)) %>%
  ungroup()

families_keep <- family_rel_raw %>%
  group_by(family) %>%
  summarise(max_rel = max(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(max_rel >= 5) %>%
  pull(family)

family_rel <- family_rel_raw %>%
  mutate(family_group = if_else(family %in% families_keep, family, "Other")) %>%
  group_by(sample_ID, sample_type, match, family_group) %>%
  summarise(family_reads = sum(family_reads, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample_ID, sample_type, match) %>%
  mutate(rel_abundance = 100 * family_reads / sum(family_reads),
         rel_abundance = if_else(is.finite(rel_abundance), rel_abundance, 0)) %>%
  ungroup()

match_order <- metadata_matched %>%
  distinct(match) %>%
  arrange(match) %>%
  pull(match)

# Use air dates for x labels.
match_dates <- metadata_matched %>%
  filter(sample_type == "Indoor Air") %>%
  distinct(match, date_parsed) %>%
  arrange(match) %>%
  mutate(label = format(date_parsed, "%b-%d"))

match_order <- match_dates$match
match_labels <- setNames(match_dates$label, match_dates$match)

family_colors <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#333666",
  "#1F78B4", "#B2DF8A", "#FB9A99", "#FDBF6F",
  "#CAB2D6"
)
family_palette <- c(head(family_colors, length(families_keep)), "#BDBDBD")
names(family_palette) <- c(families_keep, "Other")

#### Panel C: Family Relative Abundance ###

panel_c <- ggplot(family_rel,
                  aes(x = factor(match, levels = match_order), y = rel_abundance, fill = family_group)) +
  geom_col(width = 0.9) +
  facet_wrap(~ sample_type, ncol = 1) +
  labs(title = "C. Family-level relative abundance per sample",
       x = "Match (timepoint)", y = "Relative abundance (%)", fill = "Family") +
  panel_theme +
  scale_fill_manual(values = family_palette) +
  scale_x_discrete(labels = match_labels) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#### Panel D: Species Per Family ###

family_species <- filtered_data_with_type %>%
  group_by(family, sample_type) %>%
  summarise(n_species = n_distinct(species), .groups = "drop")

top_families_species <- family_species %>%
  group_by(family) %>%
  summarise(total_species = sum(n_species), .groups = "drop") %>%
  arrange(desc(total_species)) %>%
  slice_head(n = 12) %>%
  pull(family)

panel_d <- family_species %>%
  filter(family %in% top_families_species) %>%
  ggplot(aes(x = reorder(family, n_species), y = n_species, fill = sample_type)) +
  geom_col(position = "dodge") +
  geom_text(data = ~ dplyr::filter(.x, n_species >= 5),
            aes(label = n_species, group = sample_type),
            position = position_dodge(width = 0.9),
            vjust = 1.2, color = "white", size = 3) +
  geom_text(data = ~ dplyr::filter(.x, n_species < 5),
            aes(label = n_species, group = sample_type),
            position = position_dodge(width = 0.9),
            vjust = -0.4, color = "black", size = 3) +
  labs(title = "B. Species recovered per family",
       x = NULL, y = "Species count") +
  panel_theme +
  scale_fill_manual(values = matrix_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  theme(legend.position = "inside",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

a_area <- area(t = 1.5, l = 1, b = 4, r = 3)
d_area <- area(t = 1.5, l = 4, b = 4, r = 12)
c_area <- area(t = 5, l = 1, b = 10, r = 9)
b_area <- area(t = 5, l = 10, b = 10, r = 12)

# Assemble panels.
layout_design <- c(a_area, d_area, c_area, b_area)
overview_figure <- panel_a + panel_d + panel_c + panel_b + plot_layout(design = layout_design)

print(overview_figure)
