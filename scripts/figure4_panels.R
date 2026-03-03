library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(patchwork)

# Paths.
metadata_csv_path <- "data/metadata_matched.csv"
subtype_csv_path  <- "data/subtype_data.csv"

# Subtype palette.
subtype_levels <- c("ambiguous", "H1N1", "H1N2", "H3N2")
subtype_colors <- c(
  "ambiguous" = "gray",
  "H1N1"      = "#AF6F00",
  "H3N2"      = "#361688"
)

# Helpers.
parse_date_safe <- function(x) {
  d <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x))
  d
}

to_base_id <- function(x) {
  x <- as.character(x)
  case_when(
    str_detect(x, "^SA_\\d+") ~ str_replace(x, "^(SA_\\d+).*", "\\1"),
    str_detect(x, "^W\\d+")   ~ str_replace(x, "^(W\\d+).*", "\\1"),
    TRUE ~ x
  )
}

extract_subtype_from_parsed <- function(x) {
  x <- as.character(x)
  tok <- str_extract(x, "(H\\d+N\\d+|ambiguous)$")
  tok <- if_else(is.na(tok) | tok == "", "ambiguous", tok)
  tok
}

# Plot theme.
fancy_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(6, "pt"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 15)
  )

# Load metadata.
metadata_raw <- read_csv(metadata_csv_path, col_types = cols(.default = "c"), show_col_types = FALSE)

metadata <- metadata_raw %>%
  select(sample_ID, date, match) %>%
  filter(
    str_detect(sample_ID, "^(W|SA_|A1$)"),
    !str_detect(sample_ID, "^(SA_?NC|NC|NC_)")
  ) %>%
  mutate(
    date_parsed    = parse_date_safe(date),
    match_num      = suppressWarnings(as.integer(match)),
    sample_type    = if_else(str_detect(sample_ID, "^W"), "Wastewater", "Indoor Air"),
    base_sample_ID = to_base_id(sample_ID)
  ) %>%
  filter(!is.na(date_parsed))

if (any(is.na(metadata$match_num))) {
  ww_map <- metadata %>%
    filter(sample_type == "Wastewater", !is.na(match_num)) %>%
    distinct(match_num, date_parsed)
  if (nrow(ww_map) > 0) {
    metadata <- metadata %>%
      rowwise() %>%
      mutate(
        match_num = if_else(is.na(match_num),
                            ww_map$match_num[which.min(abs(as.numeric(date_parsed - ww_map$date_parsed)))],
                            match_num)
      ) %>%
      ungroup()
  }
}
metadata <- metadata %>% filter(!is.na(match_num))

# Drop W159.
drop_match_nums_W159 <- metadata %>%
  filter(base_sample_ID == "W159") %>%
  distinct(match_num) %>%
  pull(match_num)

if (length(drop_match_nums_W159) > 0) {
  metadata <- metadata %>% filter(!(match_num %in% drop_match_nums_W159))
}

meta_key <- metadata %>%
  distinct(base_sample_ID, sample_type, match_num)

match_ref <- metadata %>%
  group_by(match_num) %>%
  arrange(desc(sample_type == "Indoor Air"), date_parsed) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(date_label_hyphen = format(date_parsed, "%b-%d")) %>%
  arrange(date_parsed) %>%
  select(match_num, date_parsed, date_label_hyphen)

x_levels   <- match_ref$match_num
date_order <- match_ref$date_label_hyphen

# Load subtype table.
sub_parsed <- read_csv(subtype_csv_path, show_col_types = FALSE) %>%
  mutate(
    read_count = as.numeric(read_count),
    sample_ID  = as.character(sample_ID),
    sample_type = if_else(str_detect(sample_ID, "^W"), "Wastewater", "Indoor Air"),
    base_sample_ID = to_base_id(sample_ID),
    subtype_token  = extract_subtype_from_parsed(subtype),
    subtype_plot   = if_else(subtype_token %in% subtype_levels, subtype_token, "ambiguous")
  ) %>%
  filter(str_detect(sample_ID, "^(SA_|W)")) %>%
  left_join(meta_key, by = c("base_sample_ID", "sample_type")) %>%
  filter(!is.na(match_num), match_num %in% x_levels) %>%
  left_join(match_ref, by = "match_num")

sub_counts <- sub_parsed %>%
  group_by(match_num, sample_type, subtype_plot) %>%
  summarise(read_count = sum(read_count, na.rm = TRUE), .groups = "drop") %>%
  complete(match_num = x_levels,
           sample_type = c("Indoor Air", "Wastewater"),
           subtype_plot = factor(subtype_levels, levels = subtype_levels),
           fill = list(read_count = 0))

sub_totals <- sub_counts %>%
  group_by(match_num, sample_type) %>%
  summarise(total_sub_reads = sum(read_count, na.rm = TRUE), .groups = "drop")

sub_rel <- sub_counts %>%
  left_join(sub_totals, by = c("match_num", "sample_type")) %>%
  left_join(match_ref, by = "match_num") %>%
  mutate(
    rel_abundance = if_else(total_sub_reads > 0, 100 * read_count / total_sub_reads, 0),
    date_label_hyphen = factor(date_label_hyphen, levels = date_order, ordered = TRUE),
    subtype_plot = factor(as.character(subtype_plot), levels = subtype_levels)
  )

panel1_df <- sub_rel %>% filter(sample_type == "Indoor Air")
panel2_df <- sub_rel %>% filter(sample_type == "Wastewater")

total_reads <- sub_totals %>%
  left_join(match_ref, by = "match_num") %>%
  mutate(date_label_hyphen = factor(date_label_hyphen, levels = date_order, ordered = TRUE))

panel1_reads <- total_reads %>% filter(sample_type == "Indoor Air")
panel2_reads <- total_reads %>% filter(sample_type == "Wastewater")

#### Panel A: Indoor Air Read Counts ####
# Stacked bar chart showing subtype-specific read counts over time for Indoor Air samples

# Plot panels.
panel1_totals <- panel1_df %>%
  group_by(date_label_hyphen) %>%
  summarise(total = sum(read_count), .groups = "drop")

# Cap reads.
y_cap <- 200000

panel1_df_capped <- panel1_df %>%
  group_by(date_label_hyphen) %>%
  mutate(
    bar_total = sum(read_count),
    scale_factor = ifelse(bar_total > y_cap, y_cap / bar_total, 1),
    read_count_capped = read_count * scale_factor
  ) %>%
  ungroup()

exceeds_cap <- panel1_totals %>%
  filter(total > y_cap) %>%
  mutate(label = paste0(round(total / 1e6, 1), "M"))

x_breaks <- levels(panel1_df$date_label_hyphen)[seq(1, length(levels(panel1_df$date_label_hyphen)), by = 3)]

p1_reads_subtype <- ggplot(panel1_df_capped, aes(x = date_label_hyphen, y = read_count_capped, fill = subtype_plot)) +
  geom_col(position = "stack", width = 0.9, color = "black", linewidth = 0.05) +
  geom_segment(data = exceeds_cap, 
               aes(x = as.numeric(date_label_hyphen) - 0.35, xend = as.numeric(date_label_hyphen) + 0.35,
                   y = y_cap * 0.97, yend = y_cap * 0.97),
               inherit.aes = FALSE, color = "white", linewidth = 2) +
  geom_segment(data = exceeds_cap,
               aes(x = as.numeric(date_label_hyphen) - 0.35, xend = as.numeric(date_label_hyphen) + 0.35,
                   y = y_cap * 0.93, yend = y_cap * 0.93),
               inherit.aes = FALSE, color = "white", linewidth = 2) +
  geom_text(data = exceeds_cap,
            aes(x = date_label_hyphen, y = y_cap, label = label),
            inherit.aes = FALSE, vjust = -0.3, fontface = "bold", size = 3.5) +
  scale_fill_manual(values = subtype_colors, breaks = subtype_levels, drop = FALSE) +
  scale_x_discrete(breaks = x_breaks) +
  scale_y_continuous(limits = c(0, y_cap * 1.1),
                     breaks = c(0, 50000, 100000, 150000, 200000),
                     labels = c("0", "50K", "100K", "150K", "200K")) +
  fancy_theme +
  labs(x = NULL, y = "Read count") +
  theme(axis.text.x = element_blank(), legend.position = "none")

#### Panel B: Indoor Air Relative Abundance ####
# Stacked bar chart showing subtype relative abundance percentages for Indoor Air samples

p1_rel_abundance <- ggplot(panel1_df, aes(x = date_label_hyphen, y = rel_abundance, fill = subtype_plot)) +
  geom_col(position = "stack", width = 0.9, color = "black", linewidth = 0.05) +
  scale_fill_manual(values = subtype_colors, breaks = subtype_levels, drop = FALSE) +
  scale_x_discrete(breaks = x_breaks) +
  fancy_theme +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank(), legend.position = "none")

#### Panel C: Wastewater Read Counts ####
# Stacked bar chart showing subtype-specific read counts over time for Wastewater samples

p2_reads_subtype <- ggplot(panel2_df, aes(x = date_label_hyphen, y = read_count, fill = subtype_plot)) +
  geom_col(position = "stack", width = 0.9, color = "black", linewidth = 0.05) +
  scale_fill_manual(values = subtype_colors, breaks = subtype_levels, drop = FALSE) +
  scale_x_discrete(breaks = x_breaks) +
  fancy_theme +
  labs(x = NULL, y = "Read count") +
  theme(axis.text.x = element_blank(), legend.position = "none")

#### Panel D: Wastewater Relative Abundance ####
# Stacked bar chart showing subtype relative abundance percentages for Wastewater samples

p2_rel_abundance <- ggplot(panel2_df, aes(x = date_label_hyphen, y = rel_abundance, fill = subtype_plot)) +
  geom_col(position = "stack", width = 0.9, color = "black", linewidth = 0.05) +
  scale_fill_manual(values = subtype_colors, breaks = subtype_levels, drop = FALSE) +
  scale_x_discrete(breaks = x_breaks) +
  fancy_theme +
  labs(x = "Date", y = "Relative Abundance (%)") +
  theme(legend.position = "none")

figure4_plot <- (p1_reads_subtype / p1_rel_abundance / p2_reads_subtype / p2_rel_abundance) +
  plot_layout(ncol = 1, heights = c(0.5, 1, 0.5, 1), guides = "collect") &
  theme(legend.position = "bottom")

library(cowplot)

# Row labels.
label_air_reads <- ggdraw() + draw_label("Indoor Air\nReads", angle = -90, size = 11, fontface = "bold")
label_air_rel   <- ggdraw() + draw_label("Indoor Air\nRelative Abundance", angle = -90, size = 11, fontface = "bold")
label_ww_reads  <- ggdraw() + draw_label("Wastewater\nReads", angle = -90, size = 11, fontface = "bold")
label_ww_rel    <- ggdraw() + draw_label("Wastewater\nRelative Abundance", angle = -90, size = 11, fontface = "bold")

right_labels <- plot_grid(
  label_air_reads, label_air_rel, label_ww_reads, label_ww_rel,
  ncol = 1, rel_heights = c(0.5, 1, 0.5, 1)
)

figure4_final <- plot_grid(
  figure4_plot, right_labels,
  ncol = 2, rel_widths = c(1, 0.08)
)

print(figure4_final)
