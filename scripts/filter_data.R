#### Common Data Filtering And Processing Functions ###
# This script contains all shared data loading and filtering logic used across
# the figure scripts. It processes raw data and saves filtered datasets to
# filtered_data/ folder for use by individual analysis scripts.

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

#### Configuration ###

# Define relative paths within upload_to_github folder
DATA_DIR <- "data"
FILTERED_DATA_DIR <- "filtered_data"

# File paths
MAIN_DATA_FILE <- file.path(DATA_DIR, "air_wastewater2.detected_virus.combined.tax.tsv")
METADATA_MATCHED_FILE <- file.path(DATA_DIR, "metadata_matched.csv")

#### Main Data Processing Functions ###

#' Load and preprocess main viral detection data
load_main_data <- function() {
  cat("Loading main viral detection data...\n")
  
  data_tsv <- read_tsv(MAIN_DATA_FILE, 
                       col_types = cols(), 
                       show_col_types = FALSE)
  
  cat("✓ Loaded", nrow(data_tsv), "rows from main data file\n")
  return(data_tsv)
}

#' Calculate reads per million (RPM)
calculate_rpm <- function(data_tsv) {
  cat("Calculating RPM...\n")
  
  data_tsv <- data_tsv %>%
    mutate(
      reads_per_million = if_else(total_filtered_reads_in_sample > 0,
                                  (reads_aligned / total_filtered_reads_in_sample) * 1e6,
                                  0)
    )
  
  cat("✓ RPM calculated\n")
  return(data_tsv)
}

#' Apply contamination and quality filters
apply_filters <- function(data_tsv) {
  cat("Applying contamination and quality filters...\n")
  
  data_filtered <- data_tsv %>%
    filter(
      # Remove known wet-lab contaminants
      !str_detect(species, "Parvovirus NIH-CQV"),
      !str_detect(species, "Chikungunya virus"),
      (genus != "Alphamesonivirus"),
      (family != "Retroviridae"),
      !str_detect(species, "unclassified Culex Bastrovirus-like virus"),
      !str_detect(species, "unclassified Mus musculus mobilized endogenous polytropic provirus"),
      
      # Amplicon contamination: remove Human mastadenovirus C with < 20% coverage
      !(species == "Human mastadenovirus C" & (covered_bases / reference_length) * 100 < 20),
      
      # Remove low-certainty hits: >=500 bp OR >=50% completeness AND >10 reads AND >1 RPM
      ((covered_bases >= 500) | (covered_bases / reference_length) * 100 >= 50) & 
        reads_aligned > 10 & reads_per_million > 1,
      
      # Remove long references with low coverage
      !(reference_length > 100000 & covered_bases < 3000),
      
      # Remove additional known contaminants
      !(species %in% c("Circovirus-like genome DCCV-4", 
                       "Human endogenous retrovirus K", 
                       "Circovirus-like genome DCCV-13")),
      
      # Keep only relevant sample IDs
      str_detect(sample_ID, "^(W|SA_)")
    )
  
  cat("✓ Filtered from", nrow(data_tsv), "to", nrow(data_filtered), "rows\n")
  return(data_filtered)
}


#' Join filtered data with matched metadata
join_data_metadata <- function(data_filtered, metadata_matched) {
  cat("Joining data with matched metadata...\n")
  
  data_matched <- data_filtered %>%
    inner_join(metadata_matched %>% select(sample_ID, match, date_parsed, sample_type), 
               by = "sample_ID")
  
  # Ensure match column is integer
  data_matched$match <- as.integer(data_matched$match)
  
  cat("✓ Final filtered dataset has", nrow(data_matched), "rows\n")
  return(data_matched)
}

#### Master Processing Function ###

#' Process all data and save to filtered_data directory
process_and_save_filtered_data <- function(force_refresh = FALSE) {
  
  cat("\n================================================================================\n")
  cat("PROCESSING AND FILTERING AIR-WASTEWATER DATA\n")
  cat("================================================================================\n\n")
  
  # Create output directory
  dir.create(FILTERED_DATA_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Check if already processed (unless forced refresh)
  main_output <- file.path(FILTERED_DATA_DIR, "data_filtered_matched.csv")
  if (!force_refresh && file.exists(main_output)) {
    cat("✓ Filtered data already exists. Use force_refresh=TRUE to reprocess.\n")
    cat("✓ Figure scripts can load:", main_output, "\n")
    return(invisible())
  }
  
  # Step 1: Load pre-filtered metadata
  if (file.exists(file.path(DATA_DIR, "metadata_matched.csv"))) {
    cat("Loading pre-filtered metadata from /data...\n")
    metadata_matched <- read_csv(file.path(DATA_DIR, "metadata_matched.csv"), show_col_types = FALSE)
  } else {
    stop("Pre-filtered metadata (metadata_matched.csv) not found in data/.")
  }

  # Step 2: Load and process main data
  data_tsv <- load_main_data()
  data_tsv <- calculate_rpm(data_tsv)
  data_filtered <- apply_filters(data_tsv)
  
  # Step 3: Join data with metadata
  data_matched <- join_data_metadata(data_filtered, metadata_matched)
  
  # Step 4: Save resulting filtered data
  cat("\nSaving filtered dataset as CSV file...\n")
  write_csv(data_matched, file.path(FILTERED_DATA_DIR, "data_filtered_matched.csv"))
  
  cat("✓ Main dataset saved to", FILTERED_DATA_DIR, "\n")
  
  # Print summary
  cat("\n================================================================================\n")
  cat("FILTERING SUMMARY\n")
  cat("================================================================================\n")
  cat("Main data (filtered & matched):", nrow(data_matched), "detections\n")
  cat("Metadata (loaded):", nrow(metadata_matched), "samples\n")
  cat("================================================================================\n\n")
}

#### Convenience Loading Functions For Figure Scripts ###

#' Load main filtered dataset
load_filtered_data <- function() {
  csv_file <- file.path(FILTERED_DATA_DIR, "data_filtered_matched.csv")
  if (!file.exists(csv_file)) {
    stop("Filtered data not found. Run process_and_save_filtered_data() first.")
  }
  return(read_csv(csv_file, show_col_types = FALSE))
}

#' Load matched metadata
load_matched_metadata <- function() {
  csv_file <- file.path(DATA_DIR, "metadata_matched.csv")
  if (!file.exists(csv_file)) {
    stop("Matched metadata not found. Run process_and_save_filtered_data() first.")
  }
  return(read_csv(csv_file, show_col_types = FALSE))
}

#' Load PCR data
load_processed_pcr_data <- function() {
  csv_file <- file.path(DATA_DIR, "pcr_data.csv")
  if (!file.exists(csv_file)) {
    # If not found, return empty df consistent with previous behavior
    return(data.frame())
  }
  return(read_csv(csv_file, show_col_types = FALSE))
}

#' Load wastewater CT data
load_processed_ww_ct_data <- function() {
  csv_file <- file.path(DATA_DIR, "ww_ct_data.csv")
  if (!file.exists(csv_file)) {
    return(data.frame())
  }
  return(read_csv(csv_file, show_col_types = FALSE))
}

#' Load subtype data
load_processed_subtype_data <- function() {
  csv_file <- file.path(DATA_DIR, "subtype_data.csv")
  if (!file.exists(csv_file)) {
    return(data.frame())
  }
  return(read_csv(csv_file, show_col_types = FALSE))
}

#' Load resistance data
load_processed_resistance_data <- function() {
  csv_file <- file.path(DATA_DIR, "resistance_data.csv")
  if (!file.exists(csv_file)) {
    return(data.frame())
  }
  return(read_csv(csv_file, show_col_types = FALSE))
}

#### Auto Execution ###

# If this script is run directly, process the data
if (!interactive()) {
  process_and_save_filtered_data()
}