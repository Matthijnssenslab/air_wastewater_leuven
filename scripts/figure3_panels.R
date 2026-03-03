library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(zoo)
library(scales)
library(patchwork)

# Load pre-filtered data (run filter_data.R once to create these)
data_matched <- read_csv("filtered_data/data_filtered_matched.csv", show_col_types = FALSE)
metadata_matched <- read_csv("data/metadata_matched.csv", show_col_types = FALSE)

# Update sample_type to match figure3's naming convention (Air instead of Indoor Air)
data_matched <- data_matched %>%
  mutate(sample_type = if_else(sample_type == "Indoor Air", "Air", sample_type))

metadata_matched <- metadata_matched %>%
  mutate(sample_type = if_else(sample_type == "Indoor Air", "Air", sample_type))

time_series_data <- data_matched %>%
  group_by(match, sample_type, genus) %>%
  summarise(
    mean_rpkmf = mean(RPKMF, na.rm = TRUE),
    mean_rpm = mean(reads_per_million, na.rm = TRUE),
    sum_reads = sum(reads_aligned, na.rm = TRUE),
    .groups = "drop"
  )

complete_matches <- time_series_data %>%
  group_by(match, genus) %>%
  summarise(
    n_matrices = n_distinct(sample_type),
    .groups = "drop"
  ) %>%
  filter(n_matrices == 2) %>%
  distinct(match) %>%
  arrange(match) %>%
  pull(match)

# Remap matches to sequential weeks.
match_remap <- data.frame(
  original_match = complete_matches,
  week_number = 1:length(complete_matches)
) %>%
  left_join(
    metadata_matched %>%
      filter(sample_type == "Air") %>%
      distinct(match, date_parsed) %>%
      rename(original_match = match),
    by = "original_match"
  )

time_series_data <- time_series_data %>%
  filter(match %in% complete_matches) %>%
  left_join(match_remap, by = c("match" = "original_match")) %>%
  select(-match) %>%
  rename(match = week_number)

# Load clinical data.
clinical_csv_path <- "data/uzleuven_pathogens_weekly_long.csv"

clinical_data <- read_csv(clinical_csv_path, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  filter(
    group %in% c("Influenzavirussen"),
    pathogen == "Influenza A virus"
  ) %>%
  group_by(date) %>%
  summarise(cases = sum(amount, na.rm = TRUE), .groups = "drop")

date_range <- metadata_matched %>%
  summarise(
    start_date = min(date_parsed, na.rm = TRUE),
    end_date = max(date_parsed, na.rm = TRUE)
  )

clinical_data_filtered <- clinical_data %>%
  filter(date >= date_range$start_date, date <= date_range$end_date)

# Load pcr data.
do_pcr_overlay <- TRUE
do_ww_ct_overlay <- TRUE

# Load pre-processed PCR data
pcr_matched <- read_csv("data/pcr_data.csv", show_col_types = FALSE)

if (nrow(pcr_matched) == 0) {
  stop("No PCR data found in data/pcr_data.csv")
}

pcr_matched_remap <- pcr_matched %>%
  filter(!is.na(match)) %>%
  left_join(match_remap, by = c("match" = "original_match")) %>%
  mutate(match = week_number) %>%
  select(cartridge, sample_ID, sample_date, tijdstip, match, ct_value, detected, detected_raw)

pcr_scaled <- pcr_matched_remap %>%
  filter(!is.na(match)) %>%
  mutate(ct_value_filled = if_else(is.na(ct_value), 40, ct_value))

if (nrow(pcr_scaled) > 0) {
  pcr_scaled <- pcr_scaled %>%
    group_by(match) %>%
    summarise(ct_value = mean(ct_value_filled, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      pcr_scaled = if (max(ct_value, na.rm = TRUE) > min(ct_value, na.rm = TRUE)) {
        (40 - ct_value) / (40 - min(ct_value, na.rm = TRUE))
      } else {
        NA_real_
      }
    )
} else {
  pcr_scaled <- pcr_scaled %>% mutate(ct_value = numeric(0), pcr_scaled = numeric(0))
}

focus_genus <- "Alphainfluenzavirus"
focus_data <- time_series_data %>%
  filter(genus == focus_genus) %>%
  arrange(match, sample_type)

if (nrow(focus_data) == 0) {
  stop("No data found for Alphainfluenzavirus.")
}

week_date_labels <- focus_data %>%
  filter(sample_type == "Air") %>%
  distinct(match, date_parsed) %>%
  mutate(label = format(date_parsed, "%b-%d")) %>%
  arrange(match)

# Align clinical weeks.
start_date <- min(metadata_matched$date_parsed, na.rm = TRUE)
origin_year <- as.integer(format(start_date, "%Y"))
week_origin <- as.Date(paste0(origin_year, "-12-17"))
if (start_date < week_origin) {
  week_origin <- as.Date(paste0(origin_year - 1, "-12-17"))
}

week_dates <- metadata_matched %>%
  filter(sample_type == "Air") %>%
  distinct(match, date_parsed) %>%
  arrange(match) %>%
  mutate(
    week_start = week_origin + 7 * floor(as.numeric(date_parsed - week_origin) / 7)
  )

clinical_weekly <- clinical_data_filtered %>%
  mutate(
    week_start = week_origin + 7 * floor(as.numeric(date - week_origin) / 7)
  ) %>%
  group_by(week_start) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

clinical_by_match <- week_dates %>%
  left_join(clinical_weekly, by = "week_start") %>%
  mutate(cases = if_else(is.na(cases), 0, cases)) %>%
  select(match, cases) %>%
  arrange(match)

# Average week 13 bins.
march_year <- if (format(week_origin, "%m") == "12") {
  as.integer(format(week_origin, "%Y")) + 1
} else {
  as.integer(format(week_origin, "%Y"))
}
march_11 <- as.Date(sprintf("%d-03-11", march_year))
march_18 <- as.Date(sprintf("%d-03-18", march_year))
week_start_11 <- week_origin + 7 * floor(as.numeric(march_11 - week_origin) / 7)
week_start_18 <- week_origin + 7 * floor(as.numeric(march_18 - week_origin) / 7)

avg_week13 <- clinical_weekly %>%
  filter(week_start %in% c(week_start_11, week_start_18)) %>%
  summarise(avg_cases = mean(cases, na.rm = TRUE), .groups = "drop")
avg_week13_cases <- if (nrow(avg_week13) == 1) avg_week13$avg_cases else NA_real_

match_map <- match_remap %>%
  select(original_match, week_number)

clinical_by_match <- clinical_by_match %>%
  left_join(match_map, by = c("match" = "original_match")) %>%
  mutate(match = week_number) %>%
  select(-week_number) %>%
  filter(!is.na(match))

weeks_with_virus <- unique(focus_data$match)
clinical_by_match <- clinical_by_match %>%
  filter(match %in% weeks_with_virus) %>%
  arrange(match)

pcr_scaled <- pcr_scaled %>%
  filter(match %in% weeks_with_virus)

if (exists("ww_ct_scaled") && nrow(ww_ct_scaled) > 0) {
  ww_ct_scaled <- ww_ct_scaled %>%
    filter(match %in% weeks_with_virus)
}

# Load pre-processed wastewater CT data
ww_ct_match <- read_csv("data/ww_ct_data.csv", show_col_types = FALSE)

if (nrow(ww_ct_match) > 0) {
  ww_ct_scaled <- ww_ct_match %>%
    filter(!is.na(match)) %>%
    left_join(match_remap, by = c("match" = "original_match")) %>%
    mutate(match = week_number) %>%
    mutate(ct_value_filled = if_else(is.na(CT_value), 40, CT_value)) %>%
    mutate(ct_value_filled = if_else(ct_value_filled >= 45, 40, ct_value_filled)) %>%
    group_by(match) %>%
    summarise(ct_value = mean(ct_value_filled, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      ww_ct_scaled = if (40 > min(ct_value, na.rm = TRUE)) {
        (40 - ct_value) / (40 - min(ct_value, na.rm = TRUE))
      } else {
        NA_real_
      }
    )
} else {
  ww_ct_scaled <- data.frame()
}

if (!is.na(avg_week13_cases) && 13 %in% weeks_with_virus) {
  if (any(clinical_by_match$match == 13)) {
    clinical_by_match <- clinical_by_match %>%
      mutate(cases = if_else(match == 13, avg_week13_cases, cases))
  } else {
    clinical_by_match <- bind_rows(
      clinical_by_match,
      data.frame(match = 13, cases = avg_week13_cases)
    ) %>%
      arrange(match)
  }
}

if (nrow(clinical_by_match) > 0) {
  clinical_by_match <- clinical_by_match %>%
    arrange(match) %>%
    mutate(
      rolling_4wk_cases = rollmean(cases, k = 4, fill = NA, align = "center", partial = TRUE)
    )
}

# Scale signals.
max_cases <- if_else(nrow(clinical_by_match) > 0,
                     max(clinical_by_match$cases, na.rm = TRUE),
                     0)
scale_factor <- if_else(max_cases > 0, 1 / max_cases, 1)

focus_log10 <- focus_data %>%
  mutate(
    transformed_rpm = case_when(
      sample_type == "Air" & mean_rpm > 0 ~ log10(mean_rpm),
      sample_type == "Air" & mean_rpm <= 0 ~ NA_real_,
      TRUE ~ mean_rpm
    )
  )

focus_log10_norm <- focus_log10 %>%
  group_by(sample_type) %>%
  mutate(
    scaled_rpm = (transformed_rpm - min(transformed_rpm, na.rm = TRUE)) /
                 (max(transformed_rpm, na.rm = TRUE) - min(transformed_rpm, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  arrange(match, sample_type)

focus_log10_rolling <- focus_log10_norm %>%
  group_by(sample_type) %>%
  arrange(match) %>%
  mutate(
    rolling_4wk = rollmean(scaled_rpm, k = 4, fill = NA, align = "center", partial = TRUE)
  ) %>%
  ungroup()

pcr_scaled <- pcr_scaled %>%
  arrange(match) %>%
  mutate(
    rolling_4wk_ct = rollmean(pcr_scaled, k = 4, fill = NA, align = "center", partial = TRUE)
  )

ww_ct_scaled <- ww_ct_scaled %>%
  arrange(match) %>%
  mutate(
    rolling_4wk_ct = rollmean(ww_ct_scaled, k = 4, fill = NA, align = "center", partial = TRUE)
  )

caption_text <- paste0(
  "Blue = Air | Brown = Wastewater | Gray line = Hospital Influenza A cases",
  " | Dark blue dashed = Indoor Air Ct | Brown dashed = WW Ct"
)

#### Panel A: Time Series Panel ####
# Main influenza A time series plot with Air/Wastewater data and clinical overlay

# Plot time series.
p_alpha_log10 <- ggplot() +
  {if (nrow(clinical_by_match) > 0) {
    list(
      geom_col(
        data = clinical_by_match,
        aes(x = match, y = cases * scale_factor),
        fill = "gray80", color = "gray50", alpha = 0.5, width = 0.6
      ),
      geom_line(
        data = clinical_by_match %>% filter(!is.na(rolling_4wk_cases)),
        aes(x = match, y = rolling_4wk_cases * scale_factor, linetype = "Clinical Cases"),
        color = "black", linewidth = 1.5
      )
    )
  } else {
    NULL
  }} +
  geom_point(
    data = focus_log10_rolling %>% select(-rolling_4wk),
    aes(x = match, y = scaled_rpm, color = sample_type),
    alpha = 0.5, size = 3.5
  ) +
  geom_line(
    data = focus_log10_rolling %>% filter(!is.na(rolling_4wk)),
    aes(x = match, y = rolling_4wk, color = sample_type),
    linewidth = 1.8
  ) +
  {if (do_pcr_overlay && nrow(pcr_scaled) > 0) {
    list(
      geom_line(
        data = pcr_scaled,
        aes(x = match, y = rolling_4wk_ct, linetype = "Indoor Air Ct"),
        color = "#1F4E79", linewidth = 1.2, alpha = 0.8
      ),
      geom_point(
        data = pcr_scaled,
        aes(x = match, y = pcr_scaled),
        color = "#1F4E79", size = 2.6, shape = 3, stroke = 1, alpha = 0.2
      )
    )
  } else {
    NULL
  }} +
  {if (do_ww_ct_overlay && nrow(ww_ct_scaled) > 0) {
    list(
      geom_line(
        data = ww_ct_scaled,
        aes(x = match, y = rolling_4wk_ct, linetype = "WW Ct"),
        color = "#8B4513", linewidth = 1.2, alpha = 0.8
      ),
      geom_point(
        data = ww_ct_scaled,
        aes(x = match, y = ww_ct_scaled),
        color = "#8B4513", size = 2.6, shape = 4, stroke = 1, alpha = 0.2
      )
    )
  } else {
    NULL
  }} +
  scale_color_manual(
    values = c("Air" = "#4393C3", "Wastewater" = "#8B4513"),
    name = "Matrix"
  ) +
  scale_linetype_manual(
    values = c("Clinical Cases" = "solid", "Indoor Air Ct" = "dashed", "WW Ct" = "dashed"),
    name = "Data"
  ) +
  scale_x_continuous(
    breaks = week_date_labels$match,
    labels = week_date_labels$label,
    expand = c(0.01, 0.01)
  ) +
  {if (nrow(clinical_by_match) > 0) {
    scale_y_continuous(
      name = "Scaled log10(RPM) (0 = min, 1 = max per matrix)",
      limits = c(0, 1.02),
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Weekly Clinical Cases (Influenza A)",
        breaks = pretty
      )
    )
  } else {
    scale_y_continuous(
      name = "Scaled log10(RPM) (0 = min, 1 = max per matrix)",
      limits = c(0, 1.02),
      breaks = seq(0, 1, 0.2)
    )
  }} +
  labs(
    title = paste0(focus_genus, ": Scaled log10(RPM) with Clinical Cases"),
    subtitle = paste0(
      "Each matrix scaled 0-1 independently (min=0, max=1)\n",
      "Points = individual week data (α=0.3) | Lines = 4-week rolling average | Gray line = clinical cases"
    ),
    x = "Date",
    caption = caption_text
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16, color = "#2b2b2b"),
    plot.subtitle = element_text(size = 11, color = "gray40", lineheight = 1.2),
    plot.caption = element_text(hjust = 0, color = "gray60", size = 10, face = "italic"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold", size = 12),
    axis.title.y.right = element_text(color = "gray60"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p_alpha_log10)

#### Panels B,C,D,E: Correlation Panels ####
# Multiple correlation plots comparing viral signals with clinical data

# Correlation panels.
if (exists("focus_log10_rolling") && nrow(focus_log10_rolling) > 0 &&
    exists("clinical_by_match") && nrow(clinical_by_match) > 0) {
  corr_theme <- theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm")
    )

  get_r <- function(df, x_col, y_col) {
    df <- df %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
    if (nrow(df) < 3 || sd(df[[x_col]]) == 0 || sd(df[[y_col]]) == 0) {
      return(NA_real_)
    }
    cor(df[[x_col]], df[[y_col]], use = "complete.obs")
  }

  get_r_p <- function(df, x_col, y_col) {
    df <- df %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
    if (nrow(df) < 3 || sd(df[[x_col]]) == 0 || sd(df[[y_col]]) == 0) {
      return(list(r = NA_real_, p = NA_real_))
    }
    r_val <- cor(df[[x_col]], df[[y_col]], use = "complete.obs")
    p_val <- cor.test(df[[x_col]], df[[y_col]])$p.value
    list(r = r_val, p = p_val)
  }

  make_corr_plot <- function(df_roll, x_col_roll, df_raw, x_col_raw, y_col, color, title_text) {
    df_roll <- df_roll %>% filter(!is.na(.data[[x_col_roll]]), !is.na(.data[[y_col]]))
    df_raw <- df_raw %>% filter(!is.na(.data[[x_col_raw]]), !is.na(.data[[y_col]]))
    if (nrow(df_roll) < 3 || nrow(df_raw) < 3) {
      return(ggplot() + theme_void() + labs(title = title_text, subtitle = "Not enough data"))
    }
    r_4wk <- get_r(df_roll, x_col_roll, y_col)
    raw_stats <- get_r_p(df_raw, x_col_raw, y_col)
    label <- paste0(
      "r_4wk = ", formatC(r_4wk, digits = 2, format = "f"),
      "\nr_raw = ", formatC(raw_stats$r, digits = 2, format = "f"),
      "\np_raw = ", formatC(raw_stats$p, digits = 2, format = "g")
    )
    x_rng <- range(df_roll[[x_col_roll]], na.rm = TRUE)
    y_rng <- range(df_roll[[y_col]], na.rm = TRUE)
    x_pos <- x_rng[1] + 0.02 * diff(x_rng)
    y_pos <- y_rng[2] - 0.08 * diff(y_rng)

    ggplot(df_roll, aes(x = .data[[x_col_roll]], y = .data[[y_col]])) +
      geom_point(color = color, alpha = 0.8, size = 3) +
      geom_smooth(method = "lm", se = TRUE, color = color, fill = color, alpha = 0.2, linewidth = 1.1) +
      annotate("text", x = x_pos, y = y_pos, label = label,
               hjust = 0, vjust = 1, size = 3.6, fontface = "bold") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, NA),
                         expand = expansion(mult = c(0, 0.05))) +
      labs(
        title = title_text,
        x = "Scaled value (4-week rolling average)",
        y = "Weekly clinical cases (Influenza A)"
      ) +
      corr_theme
  }

  corr_air <- focus_log10_rolling %>%
    filter(sample_type == "Air") %>%
    select(match, rolling_4wk) %>%
    left_join(clinical_by_match, by = "match")

  corr_ww <- focus_log10_rolling %>%
    filter(sample_type == "Wastewater") %>%
    select(match, rolling_4wk) %>%
    left_join(clinical_by_match, by = "match")

  corr_air_raw <- focus_log10_norm %>%
    filter(sample_type == "Air") %>%
    select(match, scaled_rpm) %>%
    left_join(clinical_by_match, by = "match")

  corr_ww_raw <- focus_log10_norm %>%
    filter(sample_type == "Wastewater") %>%
    select(match, scaled_rpm) %>%
    left_join(clinical_by_match, by = "match")

  corr_ct <- pcr_scaled %>%
    select(match, rolling_4wk_ct) %>%
    left_join(clinical_by_match, by = "match")

  corr_ct_raw <- pcr_scaled %>%
    select(match, pcr_scaled) %>%
    left_join(clinical_by_match, by = "match")

  corr_ww_ct <- if (do_ww_ct_overlay && nrow(ww_ct_scaled) > 0) {
    ww_ct_scaled %>%
      select(match, rolling_4wk_ct) %>%
      left_join(clinical_by_match, by = "match")
  } else {
    tibble(match = integer(), rolling_4wk_ct = numeric(), cases = numeric())
  }

  corr_ww_ct_raw <- if (do_ww_ct_overlay && nrow(ww_ct_scaled) > 0) {
    ww_ct_scaled %>%
      select(match, ww_ct_scaled) %>%
      left_join(clinical_by_match, by = "match")
  } else {
    tibble(match = integer(), ww_ct_scaled = numeric(), cases = numeric())
  }

  p_corr_air <- make_corr_plot(corr_air, "rolling_4wk", corr_air_raw, "scaled_rpm", "cases", "#4393C3", "Air RPM")
  p_corr_ww <- make_corr_plot(corr_ww, "rolling_4wk", corr_ww_raw, "scaled_rpm", "cases", "#8B4513", "Wastewater RPM")
  p_corr_ct <- make_corr_plot(corr_ct, "rolling_4wk_ct", corr_ct_raw, "pcr_scaled", "cases", "#1F4E79", "Indoor Air Ct")
  p_corr_ww_ct <- make_corr_plot(corr_ww_ct, "rolling_4wk_ct", corr_ww_ct_raw, "ww_ct_scaled", "cases", "#8B4513", "WW Ct (Leuven)")

  p_corr <- (p_corr_air | p_corr_ww | p_corr_ct | p_corr_ww_ct) +
    plot_annotation(
      title = "Correlations vs Clinical Cases (4-week rolling average)",
      subtitle = "p_raw from weekly unsmoothed values; r_4wk from rolling averages."
    )

  print(p_corr)

combined_fig3 <- p_alpha_log10 / p_corr +
  plot_layout(heights = c(2, 1))

  # Show combined figure.
  print(combined_fig3)
} else {
  warning("Skipping correlation plot (missing focus or clinical data)")
}
