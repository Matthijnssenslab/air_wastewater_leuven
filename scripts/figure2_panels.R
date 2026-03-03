library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(patchwork)
library(ggrepel)
library(vegan)

# Load pre-filtered data (run filter_data.R once to create these)
data_matched <- read_csv("filtered_data/data_filtered_matched.csv", show_col_types = FALSE)
metadata_matched <- read_csv("data/metadata_matched.csv", show_col_types = FALSE)

# Genus groups.
respiratory_genera <- c(
  "unclassified Pneumoviridae genus",
  "Orthopneumovirus",
  "Alphainfluenzavirus",
  "Betainfluenzavirus",
  "Respirovirus",
  "Betacoronavirus",
  "Cytomegalovirus",
  "Lymphocryptovirus",
  "Simplexvirus",
  "Enterovirus",
  "Metapneumovirus"
)

enteric_genera <- c(
  "Kobuvirus",
  "Sapovirus",
  "Salivirus",
  "Norovirus",
  "unclassified Astroviridae genus",
  "Mamastrovirus",
  "Rotavirus"
)

segmented_genera <- c("Alphainfluenzavirus", "Betainfluenzavirus", "Rotavirus")

# Build completeness table.
all_sequences_data <- data_matched %>%
  filter(genus %in% c(respiratory_genera, enteric_genera)) %>%
  mutate(
    genome_completeness = (covered_bases / reference_length) * 100,
    genus_standardized = case_when(
      genus %in% c("unclassified Pneumoviridae genus", "Orthopneumovirus") ~ "Orthopneumovirus",
      genus %in% c("unclassified Astroviridae genus", "Mamastrovirus") ~ "Mamastrovirus",
      TRUE ~ genus
    ),
    genus_clean = str_replace_all(genus_standardized, "unclassified ", "unc. "),
    virus_type = case_when(
      genus_standardized %in% respiratory_genera ~ "Respiratory",
      genus_standardized %in% enteric_genera ~ "Enteric",
      TRUE ~ "Other"
    ),
    is_segmented = genus_standardized %in% segmented_genera
  ) %>%
  select(sample_ID, genus_standardized, genus_clean, species, sample_type,
         genome_completeness, virus_type, is_segmented, match) %>%
  rename(genus = genus_standardized, completeness_value = genome_completeness)

# Build labels.
genus_sample_counts <- all_sequences_data %>%
  group_by(genus_clean, sample_type, virus_type) %>%
  summarise(
    n_samples = n_distinct(sample_ID),
    n_sequences = n(),
    .groups = "drop"
  ) %>%
  group_by(genus_clean, virus_type) %>%
  summarise(
    air_samples = ifelse(any(sample_type == "Indoor Air"),
                         n_samples[sample_type == "Indoor Air"][1], 0),
    ww_samples = ifelse(any(sample_type == "Wastewater"),
                        n_samples[sample_type == "Wastewater"][1], 0),
    air_sequences = ifelse(any(sample_type == "Indoor Air"),
                           n_sequences[sample_type == "Indoor Air"][1], 0),
    ww_sequences = ifelse(any(sample_type == "Wastewater"),
                          n_sequences[sample_type == "Wastewater"][1], 0),
    .groups = "drop"
  ) %>%
  mutate(
    label_with_counts = paste0(
      genus_clean, "\n",
      "(Air: ", air_samples, "/", air_sequences, " | ",
      "WW: ", ww_samples, "/", ww_sequences, ")"
    )
  )

all_sequences_data <- all_sequences_data %>%
  left_join(genus_sample_counts, by = c("genus_clean", "virus_type"))

combined_plot_data <- all_sequences_data %>%
  filter(virus_type %in% c("Respiratory", "Enteric")) %>%
  filter(!is.na(completeness_value), completeness_value >= 0, completeness_value <= 100)

order_stats <- combined_plot_data %>%
  group_by(virus_type, genus_clean) %>%
  summarise(
    median_comp = median(completeness_value, na.rm = TRUE),
    total_detections = n(),
    .groups = "drop"
  ) %>%
  arrange(virus_type, desc(median_comp), desc(total_detections))

genus_label_map <- genus_sample_counts %>%
  mutate(
    genus_label = paste0(
      genus_clean, "\n",
      "(Air: ", air_samples, "/", air_sequences, " | ",
      "WW: ", ww_samples, "/", ww_sequences, ")"
    )
  ) %>%
  select(genus_clean, genus_label)

label_map <- setNames(genus_label_map$genus_label, genus_label_map$genus_clean)

genus_counts_long <- genus_sample_counts %>%
  pivot_longer(
    cols = c(air_samples, air_sequences, ww_samples, ww_sequences),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    sample_type = if_else(str_starts(metric, "air_"), "Indoor Air", "Wastewater"),
    metric = case_when(
      metric %in% c("air_samples", "ww_samples") ~ "n",
      TRUE ~ "ndt"
    )
  ) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(
    n = replace_na(n, 0L),
    ndt = replace_na(ndt, 0L)
  )

resp_levels <- order_stats %>%
  filter(virus_type == "Respiratory") %>%
  pull(genus_clean)

enteric_levels <- order_stats %>%
  filter(virus_type == "Enteric") %>%
  pull(genus_clean)

plot_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13, face = "bold", color ="black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    axis.line = element_line(color = "black", size = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold", color = "white"),
    strip.background = element_rect(fill = "black", color = "black", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(0.8, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Plot helper.
make_genus_plot <- function(df, levels_order, panel_title) {
  df <- df %>%
    mutate(
      genus_ordered = factor(genus_clean, levels = levels_order),
      genus_num = as.numeric(genus_ordered),
      x_dodge = genus_num + ifelse(sample_type == "Indoor Air", -0.2, 0.2)
    )
  
  df_violin <- df %>%
    group_by(genus_ordered, sample_type) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  ggplot(df, aes(x = x_dodge, y = completeness_value, fill = sample_type, group = interaction(genus_ordered, sample_type))) +
    geom_violin(
      data = df_violin,
      alpha = 0.6, trim = TRUE, scale = "width", width = 0.35
    ) +
    geom_point(
      aes(color = sample_type, shape = is_segmented),
      alpha = 0.7, size = 1.6, stroke = 0.3,
      position = position_jitter(width = 0.08, height = 0)
    ) +
    scale_fill_manual(values = c("Wastewater" = "#8B5A2B", "Indoor Air" = "#377EB8")) +
    scale_color_manual(values = c("Wastewater" = "#8B5A2B", "Indoor Air" = "#377EB8")) +
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 19),
                       name = "Virus Type",
                       labels = c("TRUE" = "Segmented", "FALSE" = "Non-segmented")) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(
      breaks = seq_along(levels_order),
      labels = label_map[levels_order]
    ) +
    labs(
      title = panel_title,
      x = "Genus",
      y = "Genome Completeness (%)",
      fill = "Sample Type",
      color = "Sample Type"
    ) +
    plot_theme
}

#### Panel A & B: Respiratory And Enteric Completeness Plots ###

respiratory_plot <- make_genus_plot(
  combined_plot_data %>% filter(virus_type == "Respiratory"),
  resp_levels,
  panel_title = "Respiratory"
)

enteric_plot <- make_genus_plot(
  combined_plot_data %>% filter(virus_type == "Enteric"),
  enteric_levels,
  panel_title = "Enteric"
)

n_resp_genera <- length(resp_levels)
n_enteric_genera <- length(enteric_levels)
total_genera <- n_resp_genera + n_enteric_genera

total_width <- 2 + (total_genera * 1.2)

top_row <- respiratory_plot | enteric_plot +
  plot_layout(
    guides = "collect",
    widths = c(n_resp_genera, n_enteric_genera)
  )

figure2_plot <- top_row &
  theme(legend.position = "bottom")

figure2_plot <- figure2_plot +
  plot_annotation(
    title = "Respiratory + Enteric Viruses: Completeness by Genus",
    subtitle = "All sequences • Genera ordered by median completeness",
    caption = "Violin width fixed per genus/sample type; points show individual sequences. Colors: Indoor Air (blue) vs Wastewater (brown)."
  )

print(figure2_plot)

#### Panel C: Exclusivity Summary Plot ###

# Exclusivity counts.
summarize_exclusivity <- function(df, level_col, level_label) {
  df %>%
    distinct(sample_type, .data[[level_col]]) %>%
    group_by(.data[[level_col]]) %>%
    summarise(
      has_air = any(sample_type == "Indoor Air"),
      has_ww = any(sample_type == "Wastewater"),
      .groups = "drop"
    ) %>%
    mutate(
      level = level_label,
      category = case_when(
        has_air & has_ww ~ "Shared",
        has_air ~ "Indoor Air only",
        has_ww ~ "Wastewater only",
        TRUE ~ "Other"
      )
    ) %>%
    filter(category != "Other")
}

excl_genus <- summarize_exclusivity(data_matched, "genus", "Genus")
excl_species <- summarize_exclusivity(data_matched, "species", "Species")

excl_summary <- bind_rows(excl_genus, excl_species) %>%
  count(level, category, name = "count") %>%
  mutate(category = factor(category, levels = c("Indoor Air only", "Shared", "Wastewater only")))


excl_plot <- ggplot(excl_summary, aes(x = level, y = count, fill = category)) +
  geom_col(width = 0.65, position = "stack") +
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4.2, hjust = 0) +
  scale_fill_manual(values = c("Indoor Air only" = "#377EB8",
                               "Shared" = "gray70",
                               "Wastewater only" = "#8B5A2B")) +
  labs(
    title = "Exclusive vs Shared Taxa by Matrix",
    subtitle = "Unique genera/species based on filtered detections (paired samples only)",
    x = NULL,
    y = "Number of unique taxa",
    fill = "Category",
    caption = "Stacked bars summarize unique vs shared taxa across matched Indoor Air and Wastewater samples."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

print(excl_plot)

#### Panel D: Nmds Plots (Genus And Species Level) ###

# Genus nmds.
genus_abundance <- data_matched %>%
  group_by(sample_ID, sample_type, genus) %>%
  summarise(abundance = sum(reads_aligned, na.rm = TRUE), .groups = "drop")

genus_wide <- genus_abundance %>%
  pivot_wider(
    names_from = genus,
    values_from = abundance,
    values_fill = list(abundance = 0)
  )

genus_matrix <- genus_wide %>%
  select(-sample_ID, -sample_type)

genus_rel <- genus_matrix / pmax(rowSums(genus_matrix), 1)

set.seed(123)
nmds_fit <- metaMDS(genus_rel, distance = "bray", k = 2, trymax = 200, autotransform = FALSE)

nmds_scores <- as.data.frame(scores(nmds_fit, display = "sites"))
nmds_scores$sample_ID <- genus_wide$sample_ID
nmds_scores$sample_type <- genus_wide$sample_type

air_match_dates <- metadata_matched %>%
  filter(sample_type == "Indoor Air") %>%
  distinct(match, date_parsed)

nmds_scores <- nmds_scores %>%
  left_join(metadata_matched %>% select(sample_ID, match), by = "sample_ID") %>%
  left_join(air_match_dates, by = "match") %>%
  mutate(
    date_parsed = as.Date(date_parsed),
    date_label = format(date_parsed, "%b-%d"),
    date_num = as.numeric(date_parsed)
  )

date_breaks_num <- scales::breaks_pretty(n = 4)(nmds_scores$date_num)
date_breaks_lbl <- format(as.Date(date_breaks_num, origin = "1970-01-01"), "%Y-%m-%d")

air_scores <- nmds_scores %>%
  filter(sample_type == "Indoor Air")

if (nrow(air_scores) > 0) {
  air_centroid <- colMeans(air_scores[, c("NMDS1", "NMDS2")], na.rm = TRUE)
  air_scores <- air_scores %>%
    mutate(dist_to_centroid = sqrt((NMDS1 - air_centroid[1])^2 + (NMDS2 - air_centroid[2])^2))
  
  outlier_air <- air_scores %>%
    arrange(desc(dist_to_centroid)) %>%
    slice(1)
  
  print(outlier_air)
}

ww_scores <- nmds_scores %>%
  filter(sample_type == "Wastewater")

if (nrow(ww_scores) > 0) {
  ww_centroid <- colMeans(ww_scores[, c("NMDS1", "NMDS2")], na.rm = TRUE)
  ww_scores <- ww_scores %>%
    mutate(dist_to_centroid = sqrt((NMDS1 - ww_centroid[1])^2 + (NMDS2 - ww_centroid[2])^2))
  
  outlier_ww <- ww_scores %>%
    arrange(desc(dist_to_centroid)) %>%
    slice(1)
  
  print(outlier_ww)
}

p_nmds <- ggplot(nmds_scores, aes(x = NMDS2, y = NMDS1, color = date_num, shape = sample_type)) +
  geom_point(size = 3, alpha = 0.85) +
  {if (exists("outlier_air") && nrow(outlier_air) == 1) {
    geom_text(data = outlier_air, aes(label = date_label),
              color = "black", size = 3, vjust = -0.8)
  } else {
    NULL
  }} +
  {if (exists("outlier_ww") && nrow(outlier_ww) == 1) {
    geom_text(data = outlier_ww, aes(label = date_label),
              color = "black", size = 3, vjust = -0.8)
  } else {
    NULL
  }} +
  stat_ellipse(type = "t", linetype = "solid", linewidth = 0.8, alpha = 0.4, color = "black") +
  scale_color_gradientn(
    colors = c("#b2182b", "#fddbc7", "#d1e5f0", "#2166ac"),
    breaks = date_breaks_num,
    labels = date_breaks_lbl
  ) +
  scale_shape_manual(values = c("Indoor Air" = 16, "Wastewater" = 17)) +
  labs(
    title = "NMDS (Genus Level, Matched Pairs)",
    subtitle = paste0("Bray-Curtis on relative abundance | Stress = ", round(nmds_fit$stress, 3)),
    x = "NMDS2",
    y = "NMDS1",
    color = "Date",
    shape = "Sample Type",
    caption = "Points colored by matched Indoor Air sampling date; ellipses show group dispersion."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(size = 12),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "bottom"
  )

print(p_nmds)

# Species nmds.
species_abundance <- data_matched %>%
  group_by(sample_ID, sample_type, species) %>%
  summarise(abundance = sum(reads_aligned, na.rm = TRUE), .groups = "drop")

species_wide <- species_abundance %>%
  pivot_wider(
    names_from = species,
    values_from = abundance,
    values_fill = list(abundance = 0)
  )

species_matrix <- species_wide %>%
  select(-sample_ID, -sample_type)

species_rel <- species_matrix / pmax(rowSums(species_matrix), 1)

set.seed(123)
nmds_species_fit <- metaMDS(species_rel, distance = "bray", k = 2, trymax = 200, autotransform = FALSE)

nmds_species_scores <- as.data.frame(scores(nmds_species_fit, display = "sites"))
nmds_species_scores$sample_ID <- species_wide$sample_ID
nmds_species_scores$sample_type <- species_wide$sample_type

nmds_species_scores <- nmds_species_scores %>%
  left_join(metadata_matched %>% select(sample_ID, match), by = "sample_ID") %>%
  left_join(air_match_dates, by = "match") %>%
  mutate(
    date_parsed = as.Date(date_parsed),
    date_label = format(date_parsed, "%b-%d"),
    date_num = as.numeric(date_parsed)
  )

date_breaks_num_sp <- scales::breaks_pretty(n = 4)(nmds_species_scores$date_num)
date_breaks_lbl_sp <- format(as.Date(date_breaks_num_sp, origin = "1970-01-01"), "%Y-%m-%d")

air_scores_sp <- nmds_species_scores %>%
  filter(sample_type == "Indoor Air")

if (nrow(air_scores_sp) > 0) {
  air_centroid_sp <- colMeans(air_scores_sp[, c("NMDS1", "NMDS2")], na.rm = TRUE)
  air_scores_sp <- air_scores_sp %>%
    mutate(dist_to_centroid = sqrt((NMDS1 - air_centroid_sp[1])^2 + (NMDS2 - air_centroid_sp[2])^2))
  
  outlier_air_sp <- air_scores_sp %>%
    arrange(desc(dist_to_centroid)) %>%
    slice(1)
  
  print(outlier_air_sp)
}

ww_scores_sp <- nmds_species_scores %>%
  filter(sample_type == "Wastewater")

if (nrow(ww_scores_sp) > 0) {
  ww_centroid_sp <- colMeans(ww_scores_sp[, c("NMDS1", "NMDS2")], na.rm = TRUE)
  ww_scores_sp <- ww_scores_sp %>%
    mutate(dist_to_centroid = sqrt((NMDS1 - ww_centroid_sp[1])^2 + (NMDS2 - ww_centroid_sp[2])^2))
  
  outlier_ww_sp <- ww_scores_sp %>%
    arrange(desc(dist_to_centroid)) %>%
    slice(1)
  
  print(outlier_ww_sp)
}

p_nmds_species <- ggplot(nmds_species_scores, aes(x = NMDS2, y = NMDS1, color = date_num, shape = sample_type)) +
  geom_point(size = 3, alpha = 0.85) +
  {if (exists("outlier_air_sp") && nrow(outlier_air_sp) == 1) {
    geom_text(data = outlier_air_sp, aes(label = date_label),
              color = "black", size = 3, vjust = -0.8)
  } else {
    NULL
  }} +
  {if (exists("outlier_ww_sp") && nrow(outlier_ww_sp) == 1) {
    geom_text(data = outlier_ww_sp, aes(label = date_label),
              color = "black", size = 3, vjust = -0.8)
  } else {
    NULL
  }} +
  stat_ellipse(type = "t", linetype = "solid", linewidth = 0.8, alpha = 0.4, color = "black") +
  scale_color_gradientn(
    colors = c("#b2182b", "#fddbc7", "darkgreen", "#2166ac"),
    breaks = date_breaks_num_sp,
    labels = date_breaks_lbl_sp
  ) +
  scale_shape_manual(values = c("Indoor Air" = 16, "Wastewater" = 17)) +
  labs(
    title = "NMDS (Species Level, Matched Pairs)",
    subtitle = paste0("Bray-Curtis on relative abundance | Stress = ", round(nmds_species_fit$stress, 3)),
    x = "NMDS2",
    y = "NMDS1",
    color = "Date",
    shape = "Sample Type",
    caption = "Points colored by matched Indoor Air sampling date; ellipses show group dispersion."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(size = 12),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "right"
  )

print(p_nmds_species)

#### Panel E: Abundance Comparison (Ordered) ###

# Genus rpm comparison.
target_genera_list <- unique(c(respiratory_genera, enteric_genera))

data_targets <- data_matched %>%
  dplyr::filter(genus %in% target_genera_list) %>%
  dplyr::mutate(
    group_key = dplyr::case_when(
      genus %in% c("unclassified Pneumoviridae genus", "Orthopneumovirus") ~ "Orthopneumovirus",
      genus %in% c("unclassified Astroviridae genus", "Mamastrovirus") ~ "Mamastrovirus",
      TRUE ~ genus
    )
  )

sample_genus_abundance <- data_targets %>%
  group_by(sample_ID, sample_type, group_key) %>%
  summarise(
    total_rpm = sum(reads_per_million, na.rm = TRUE),
    .groups = "drop"
  )

lfc_stats <- sample_genus_abundance %>%
  mutate(log_rpm = log10(total_rpm)) %>%
  group_by(group_key, sample_type) %>%
  summarise(median_log = median(log_rpm), .groups = "drop") %>%
  pivot_wider(names_from = sample_type, values_from = median_log, values_fill = 0) %>%
  mutate(
    lfc_median = `Indoor Air` - `Wastewater`,
    lfc_label = sprintf("%+.2f", lfc_median),
    color_group = ifelse(lfc_median > 0, "#2166ac", "#8c510a")
  ) %>%
  arrange(lfc_median)

sample_genus_abundance <- sample_genus_abundance %>%
  mutate(group_key = factor(group_key, levels = lfc_stats$group_key))

max_val <- max(log10(sample_genus_abundance$total_rpm), na.rm = TRUE)
text_x_pos <- max_val * 1.15

p_abundance_ordered <- ggplot(sample_genus_abundance, aes(x = group_key, y = log10(total_rpm), fill = sample_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, 
               position = position_dodge2(width = 0.8, preserve = "single")) +
  geom_point(aes(color = sample_type), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.5, alpha = 0.8, show.legend = FALSE) +
  geom_text(data = lfc_stats, 
            aes(x = group_key, y = text_x_pos, label = lfc_label, color = I(color_group)), 
            inherit.aes = FALSE, size = 3.5, fontface = "bold", hjust = 1) +
  coord_flip() + 
  scale_fill_manual(values = c("Indoor Air" = "#2166ac", "Wastewater" = "#8c510a"), name = "Sample Type") +
  scale_color_manual(values = c("Indoor Air" = "#164473", "Wastewater" = "#5e3606")) +
  labs(
    title = "Viral Abundance by Genus",
    subtitle = "Ordered by Median Log Fold Change (WW-heavy bottom -> Air-heavy top)",
    y = "Log10(Total RPM)",
    x = NULL
  ) +
  scale_y_continuous(breaks = scales::breaks_width(1), limits = c(0, text_x_pos)) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

print(p_abundance_ordered)

# Dispersion tests.
bray_genus <- vegdist(genus_rel, method = "bray")
disp_genus <- betadisper(bray_genus, group = genus_wide$sample_type)
disp_genus_test <- permutest(disp_genus, permutations = 999)

cat("\nBeta dispersion (genus level):\n")
print(anova(disp_genus))
print(disp_genus_test)
print(tapply(disp_genus$distances, genus_wide$sample_type, summary))

bray_species <- vegdist(species_rel, method = "bray")
disp_species <- betadisper(bray_species, group = species_wide$sample_type)
disp_species_test <- permutest(disp_species, permutations = 999)

cat("\nBeta dispersion (species level):\n")
print(anova(disp_species))
print(disp_species_test)
print(tapply(disp_species$distances, species_wide$sample_type, summary))
