#### Figure 5: Genus-specific Viral Composition ###
# Analyzes relative abundance and RPM for key genera: Mastadenovirus, Mamastrovirus, 
# Norovirus, Enterovirus, Rotavirus, and Papillomaviridae.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(patchwork)
library(RColorBrewer)
library(viridis)

#### Data Loading ###

# Load filtered data and metadata
data_filtered <- read_csv("filtered_data/data_filtered_matched.csv", show_col_types = FALSE)
metadata_matched <- read_csv("data/metadata_matched.csv", show_col_types = FALSE)

# Prepare data for plotting
data_tsv <- data_filtered
metadata <- metadata_matched %>%
  filter(!is.na(sample_ID))

# Calculate RPM if missing
if (!"reads_per_million" %in% names(data_tsv)) {
  data_tsv <- data_tsv %>%
    mutate(
      reads_per_million = if_else(total_filtered_reads_in_sample > 0,
                                  (reads_aligned / total_filtered_reads_in_sample) * 1e6,
                                  0)
    )
}

data_filtered <- data_tsv

# ---- Metadata + matched pairs ----
metadata_dates <- metadata %>%
  mutate(
    sample_ID = trimws(sample_ID),
    date_parsed = as.Date(date, format = "%d/%m/%Y"),
    match = suppressWarnings(as.integer(match)),
    sample_type = ifelse(grepl("^SA_", sample_ID, ignore.case = TRUE), "Indoor Air", "Wastewater")
  ) %>%
  filter(!is.na(match))

# Identify pairs with both Air and WW samples
matched_pairs <- metadata_dates %>%
  group_by(match) %>%
  summarise(
    has_air = any(sample_type == "Indoor Air"),
    has_ww = any(sample_type == "Wastewater"),
    .groups = "drop"
  ) %>%
  filter(has_air & has_ww) %>%
  pull(match)

# Filter metadata to strictly matched pairs
metadata_matched <- metadata_dates %>%
  filter(match %in% matched_pairs)

# Join virus data with strict metadata
# Note: Cleaning metadata columns from data first to avoid duplication conflicts
data_matched <- data_filtered %>%
  select(-any_of(c("match", "date_parsed", "sample_type"))) %>%
  inner_join(metadata_matched %>% select(sample_ID, match, date_parsed, sample_type),
             by = "sample_ID")

# Ensure match column is integer
data_matched$match <- as.integer(data_matched$match)

# ---- Select target genera / families ----
target_genera_list <- c("Mastadenovirus", "Mamastrovirus", "Norovirus", "Enterovirus", "Rotavirus")
target_family_papillo <- "Papillomaviridae"

data_targets <- data_matched %>%
  filter(genus %in% target_genera_list | family == target_family_papillo | str_detect(genus, "apillomavirus")) %>%
  # Remove any Porcine species
  filter(!str_detect(species, regex("Porcine", ignore_case = TRUE))) %>%
  mutate(
    group_key = case_when(
      genus == "Mastadenovirus" ~ "Mastadenovirus",
      genus == "Mamastrovirus" ~ "Mamastrovirus",
      genus == "Norovirus" ~ "Norovirus",
      genus == "Enterovirus" ~ "Enterovirus",
      genus == "Rotavirus" ~ "Rotavirus",
      family == "Papillomaviridae" | str_detect(genus, "apillomavirus") ~ "Papillomavirus",
      TRUE ~ "Other"
    )
  )

# ---- Aggregation Logic ----

# Helper to determine grouping column
get_grouping_col <- function(g_key, available_cols) {
  # Default priority: Species
  col_name <- "species"
  
  has_strain <- "strain" %in% available_cols
  
  if (g_key %in% c("Norovirus", "Mamastrovirus", "Rotavirus")) {
    if (has_strain) col_name <- "strain"
  } else if (g_key == "Papillomavirus") {
    col_name <- "genus"
  } else if (g_key == "Enterovirus") {
    # "Use both" -> Composite key if strain exists
    if (has_strain) col_name <- "both" 
  }
  
  return(col_name)
}

# Add 'strain' placeholder if missing (so code doesn't break, using species as fallback)
real_strain_exists <- "strain" %in% names(data_targets)
if (!real_strain_exists) {
  data_targets$strain <- data_targets$species
}

# Create composite column for Enterovirus logic
# Avoid redundant "Species - Species" if strain is missing or identical
data_targets <- data_targets %>%
  mutate(
    species_strain = if_else(species == strain, species, paste(species, strain, sep = " - "))
  )

agg_list <- list()

for (g in unique(data_targets$group_key)) {
  
  # Determine which column to use based on key
  col_name <- "species"
  
  if (g %in% c("Norovirus", "Mamastrovirus")) {
    # Prefer strain if available (or if we populated it as fallback)
    col_name <- "strain"
  } else if (g == "Papillomavirus") {
    col_name <- "genus"
  } else if (g == "Enterovirus") {
    # Use composite
    col_name <- "species_strain"
  }
  
  # Group and aggregate
  sub_agg <- data_targets %>%
    filter(group_key == g) %>%
    # Rename target column to generic 'variant_label' for plotting
    rename(variant_label = !!sym(col_name))
  
  # Special fix for Astrovirus 1 variants (e.g. Beijing)
  if (g == "Mamastrovirus") {
    sub_agg <- sub_agg %>%
      mutate(variant_label = if_else(str_detect(variant_label, "strovirus 1"), "Human astrovirus 1", variant_label))
  }

  # ---- Harmonize Viral Taxonomy ----
  # Rename specific variants to distinct species names (e.g., Mastadenovirus, Norovirus)
  sub_agg <- sub_agg %>%
    mutate(variant_label = case_when(
      str_detect(variant_label, "Human mastadenovirus A") ~ "Mastadenovirus adami",
      str_detect(variant_label, "Human mastadenovirus B") ~ "Mastadenovirus blackbeardi",
      str_detect(variant_label, "Human mastadenovirus C") ~ "Mastadenovirus caesari",
      str_detect(variant_label, "Human mastadenovirus D") ~ "Mastadenovirus dominans",
      str_detect(variant_label, "Human mastadenovirus E") ~ "Mastadenovirus exoticum",
      str_detect(variant_label, "Human mastadenovirus F") ~ "Mastadenovirus faecale",
      str_detect(variant_label, "Human astrovirus 1|Mamastrovirus 1") ~ "Mamastrovirus humanum",
      str_detect(variant_label, "Norwalk virus") ~ "Norovirus norwalkense",
      TRUE ~ variant_label
    ))

  # Aggregate counts by sample and variant
  sub_agg <- sub_agg %>%
    group_by(sample_ID, match, sample_type, group_key, variant_label) %>%
    summarise(
      reads = sum(reads_aligned, na.rm = TRUE),
      rpm = sum(reads_per_million, na.rm = TRUE),
      .groups = "drop"
    )
  
  agg_list[[g]] <- sub_agg
}

# Combine all genus data and calculate relative abundance
agg_data <- bind_rows(agg_list) %>%
  group_by(sample_ID, group_key) %>%
  mutate(rel_abundance = reads / sum(reads) * 100) %>%
  ungroup()

# ---- Filter Low Abundance Variants ----
# Identify major variants (>5% abundance in at least one sample)
major_variants <- agg_data %>%
  group_by(group_key, variant_label) %>%
  summarise(max_val = max(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(max_val > 5)

# Group minor variants as "Other" and re-aggregate
agg_data <- agg_data %>%
  mutate(
    variant_final = if_else(
      paste(group_key, variant_label) %in% paste(major_variants$group_key, major_variants$variant_label),
      variant_label,
      "Other"
    )
  ) %>%
  group_by(sample_ID, match, sample_type, group_key, variant_final) %>%
  summarise(
    reads = sum(reads, na.rm = TRUE),
    rpm = sum(rpm, na.rm = TRUE),
    rel_abundance = sum(rel_abundance, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  rename(variant_label = variant_final)

# ---- Plotting Functions ----

plot_theme_base <- theme_minimal(base_size = 10) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
    legend.text = element_text(color = "black", size = 7),
    legend.title = element_text(color = "black", size = 8, face = "bold"),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.key.size = unit(0.3, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    plot.margin = margin(2, 2, 2, 2)
  )

# Prepare x-axis labels: Map Match ID -> Date
match_date_map <- metadata_matched %>%
  filter(sample_type == "Indoor Air") %>%
  select(match, date_parsed) %>%
  distinct() %>%
  arrange(match) %>%
  mutate(date_label = format(date_parsed, "%b %d"))

match_levels <- match_date_map$match
match_labels <- match_date_map$date_label

# Define Genus-specific color palettes
col_palettes <- list(
  "Mastadenovirus" = c("#E69F00", "#56B4E9", "#009E73","#0072B2", "#F0E442", "#D55E00", "#CC79A7", "#661100", "#882255", "#332288"),
  "Mamastrovirus" = brewer.pal(8, "Dark2"),
  "Norovirus" = brewer.pal(8, "Paired"),
  "Enterovirus" = brewer.pal(8, "Set1"),
  "Papillomavirus" = brewer.pal(10, "Set3"),
  "Rotavirus" = c("#E41A1C", "#377EB8", "#4DAF4A")
)

make_genus_column <- function(g_key, genus_title, only_rpm = FALSE) {
  df_sub <- agg_data %>% filter(group_key == g_key)
  
  if(nrow(df_sub) == 0) return(NULL)
  
  # Clean labels for legend
  df_sub <- df_sub %>%
    mutate(variant_clean = str_replace(variant_label, "Human ", "")) %>%
    mutate(variant_clean = str_trunc(variant_clean, 40))
  
  # Order factors by read count
  var_summ <- df_sub %>%
    group_by(variant_clean) %>%
    summarise(total = sum(reads)) %>%
    arrange(desc(total))
  
  non_other <- var_summ %>% filter(variant_clean != "Other") %>% pull(variant_clean)
  has_other <- "Other" %in% var_summ$variant_clean
  final_levels <- if(has_other) c(non_other, "Other") else non_other
  df_sub$variant_clean <- factor(df_sub$variant_clean, levels = final_levels)
  
  # Assign colors dynamically
  base_pal <- col_palettes[[g_key]]
  if(is.null(base_pal)) base_pal <- brewer.pal(8, "Set3")
  
  if(has_other) {
    n_colors_needed <- length(non_other)
    if(n_colors_needed > 0) {
      if (n_colors_needed > length(base_pal)) main_colors <- colorRampPalette(base_pal)(n_colors_needed)
      else main_colors <- base_pal[1:n_colors_needed]
    } else main_colors <- c()
    my_colors <- c(main_colors, "grey80")
    names(my_colors) <- c(non_other, "Other")
  } else {
    n_vars <- length(final_levels)
    if (n_vars > length(base_pal)) my_colors <- colorRampPalette(base_pal)(n_vars)
    else my_colors <- base_pal[1:n_vars]
    names(my_colors) <- final_levels
  }
  
  # Helper for single panel plot
  make_panel <- function(data, y_col, y_lab, limits = NULL, show_x = FALSE, plot_title = NULL) {
    p <- ggplot(data, aes(x = factor(match, levels = match_levels), y = !!sym(y_col), fill = variant_clean)) +
      geom_col(width = 0.8) +
      scale_fill_manual(values = my_colors, name = "Variant") +
      scale_x_discrete(labels = match_labels, drop = FALSE) +
      scale_y_continuous(expand = c(0,0), limits = limits) +
      labs(y = y_lab, title = plot_title) +
      plot_theme_base
    
    if(show_x) {
      p <- p + theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, color = "black"),
        axis.title.x = element_text(size = 9, color = "black")
      ) + labs(x = "Collection Date")
    }
    return(p)
  }
  
  # Split data
  df_air <- df_sub %>% filter(sample_type == "Indoor Air")
  df_ww <- df_sub %>% filter(sample_type == "Wastewater")
  
  # Generate Panels
  p1 <- make_panel(df_air, "rpm", "Air RPM", plot_title = genus_title)
  p3 <- make_panel(df_ww, "rpm", "WW RPM", show_x = only_rpm)
  
  if (only_rpm) {
    return(
      (p1 / p3) + 
      plot_layout(guides = 'collect', heights = c(1, 1.2)) &
      theme(legend.position = "bottom", legend.direction = "vertical")
    )
  }
  
  p2 <- make_panel(df_air, "rel_abundance", "Air %", limits = c(0, 105))
  p4 <- make_panel(df_ww, "rel_abundance", "WW %", limits = c(0, 105), show_x = TRUE)
  
  # Return combined plot column
  (p1 / p2 / p3 / p4) + 
    plot_layout(guides = 'collect', heights = c(0.5, 1, 0.5, 1.2)) &
    theme(legend.position = "bottom", legend.direction = "vertical")
}

# Generate columns
c_adeno <- make_genus_column("Mastadenovirus", "Mastadenovirus")
c_astro <- make_genus_column("Mamastrovirus", "Mamastrovirus")
c_noro <- make_genus_column("Norovirus", "Norovirus")
c_entero <- make_genus_column("Enterovirus", "Enterovirus")
c_pap <- make_genus_column("Papillomavirus", "Papillomavirus")
c_rota <- make_genus_column("Rotavirus", "Rotavirus (Species)", only_rpm = TRUE) # RPM only

# Define groupings
main_genera <- list(c_adeno, c_entero)
supp_genera <- list(c_astro, c_noro, c_pap)
supp_rota <- list(c_rota)

# Filter NULLs in case data is missing
main_genera <- main_genera[!sapply(main_genera, is.null)]
supp_genera <- supp_genera[!sapply(supp_genera, is.null)]
supp_rota <- supp_rota[!sapply(supp_rota, is.null)]

if (length(main_genera) == 0 && length(supp_genera) == 0 && length(supp_rota) == 0) {
  stop("No data found for any target genus.")
}

# ---- Main Figure (Mastadenovirus + Enterovirus) ----
if (length(main_genera) > 0) {
  fig_main <- wrap_plots(main_genera, nrow = 1) +
    plot_annotation(
      title = "Figure 5: Viral Composition & Abundance (Mastadenovirus & Enterovirus)",
      subtitle = "Top=Air RPM, 2nd=Air %, 3rd=WW RPM, Bottom=WW %",
      caption = "Filtered for >5% per sample. RPM=Reads Per Million. Matched Pairs.",
      theme = theme(text = element_text(color = "black"))
    )
  
  print(fig_main)
}

# ---- Supplementary Figure 1 (Others) ----
if (length(supp_genera) > 0) {
  fig_supp <- wrap_plots(supp_genera, nrow = 1) +
    plot_annotation(
      title = "Supplementary Figure 1: Viral Composition (Mamastrovirus, Norovirus, Papillomavirus)",
      subtitle = "Top=Air RPM, 2nd=Air %, 3rd=WW RPM, Bottom=WW %",
      caption = "Filtered for >5% per sample. RPM=Reads Per Million. Matched Pairs.",
      theme = theme(text = element_text(color = "black"))
    )
  
  print(fig_supp)
}

# ---- Supplementary Figure 2 (Rotavirus) ----
if (length(supp_rota) > 0) {
  fig_supp2 <- wrap_plots(supp_rota, nrow = 1) +
    plot_annotation(
      title = "Supplementary Figure 2: Rotavirus Abundance (by Species)",
      subtitle = "Reads Per Million (RPM) in Indoor Air vs Wastewater",
      caption = "Filtered for >5% per sample. RPM=Reads Per Million. Matched Pairs.",
      theme = theme(text = element_text(color = "black"))
    )
  
  print(fig_supp2)
}
