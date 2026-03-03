#### Figure 4 Part 2: Influenza A Resistance Analysis Pipeline ###
# Analyzes resistance mutations in H1N1 samples using BAM files and reference sequences
# Generates coverage plots and mutation analysis for HA, NA, M2, and PA genes

library(jsonlite)
library(dplyr)
library(ggplot2)
library(cowplot)
library(readr)

# Paths - Updated for upload_to_github folder structure.
project_dir <- getwd()  # Should be the upload_to_github directory
bam_dir <- file.path(project_dir, "bam")
h1n1_ref_dir <- file.path(project_dir, "references/H1N1pdm09_A_California_07_2009")

python_script <- file.path(project_dir, "scripts/helper_scripts/analyze_mutations.py")
depth_script <- file.path(project_dir, "scripts/helper_scripts/generate_depths.py")
python_cmd <- "python3"

# Tool paths - try to detect automatically, fall back to system PATH
samtools_cmd <- Sys.getenv("SAMTOOLS", unset = "")
if (samtools_cmd == "" || !file.exists(samtools_cmd)) {
  # Try common conda locations first
  conda_samtools <- file.path(Sys.getenv("HOME"), "anaconda3/envs/bam-filter/bin/samtools")
  if (file.exists(conda_samtools)) {
    samtools_cmd <- conda_samtools
  } else {
    samtools_cmd <- "samtools"  # Use system PATH
  }
}

bcftools_cmd <- Sys.getenv("BCFTOOLS", unset = "")
if (bcftools_cmd == "" || !file.exists(bcftools_cmd)) {
  # Try common conda locations first
  conda_bcftools <- file.path(Sys.getenv("HOME"), "anaconda3/envs/bam-filter/bin/bcftools")
  if (file.exists(conda_bcftools)) {
    bcftools_cmd <- conda_bcftools
  } else {
    bcftools_cmd <- "bcftools"  # Use system PATH
  }
}
resistance_csv <- file.path(project_dir, "data/inf_resistance/h1n1_resist.csv")
resistance_analysis_dir <- file.path(project_dir, "resistance_analysis")
dir.create(resistance_analysis_dir, showWarnings = FALSE, recursive = TRUE)
output_csv <- file.path(resistance_analysis_dir, "figure4_2_data.csv")
json_config_path <- file.path(resistance_analysis_dir, "analysis_config.json")

#### Configuration Setup ###
# Analysis configuration for H1N1 sample SA_10

# Config.
analysis_config <- list(
  samples = list(
    list(
      name = "H1N1 (SA_10)",
      virus_type = "A(H1N1)pdm09",
      bam = file.path(resistance_analysis_dir, "SA_10.trimmed.bam"),
      consensus = file.path(resistance_analysis_dir, "SA_10.consensus_from_bam.fasta"),
      segments = list(
        "HA" = list(
          gene = "HA",
          ref_gb = file.path(h1n1_ref_dir, "segment_4_NC_026433.gb"),
          cons_acc = "NC_026433.1"
        ),
        "NA" = list(
          gene = "NA",
          ref_gb = file.path(h1n1_ref_dir, "segment_6_NC_026434.gb"),
          cons_acc = "NC_026434.1"
        ),
        "M2" = list(
          gene = "M2",
          ref_gb = file.path(h1n1_ref_dir, "segment_7_NC_026431.gb"),
          cons_acc = "NC_026431.1"
        ),
        "PA" = list(
          gene = "PA",
          ref_gb = file.path(h1n1_ref_dir, "segment_3_NC_026437.gb"),
          cons_acc = "NC_026437.1"
        )
      )
    )
  )
)

#### Pipeline Execution ###
# Run consensus generation and mutation analysis pipeline
# Skip if output data already exists

if (file.exists(output_csv)) {
  message("✓ Analysis results already exist, loading from: ", output_csv)
  message("  (Delete resistance_analysis/figure4_2_data.csv to force re-analysis)")
  plot_data <- read_csv(output_csv, show_col_types = FALSE)
} else {
  message("Result CSV not found. Checking if intermediate files exist (GitHub Mode)...")

  # GitHub Mode: Skip BAM processing if we have the depth/consensus/mutation files
  # This allows running the pipeline without the heavy BAM files
  h1n1_depth <- file.path(resistance_analysis_dir, "SA_10.depth.txt")
  h1n1_mutations <- file.path(resistance_analysis_dir, "all_mutations.csv")
  
  if (file.exists(h1n1_depth) && file.exists(h1n1_mutations)) {
     message("✓ Intermediate data found! Skipping BAM processing.")
     
     # Skip directly to plotting logic or simple aggregation
     # (The full re-analysis requires python scripts that might need BAMs if not careful,
     # so we assume if all_mutations.csv exists, we can largely rely on it)
     
     # NOTE: The original script runs python to generate alignment summaries.
     # If those are also cached, we load them.
  }
  
  if (!file.exists(h1n1_depth)) {
      message("Running full pipeline (Requires BAM files)...")
  
      message("Genering configuration JSON...")
      write_json(analysis_config, json_config_path, auto_unbox = TRUE, pretty = TRUE)
      
      #### Consensus Generation ###

# Generate filtered consensus sequences from BAM files

# Build consensus.
message("Checking for existing consensus files...")

h1n1_ref_fa <- file.path(h1n1_ref_dir, "ref.fasta")
h1n1_bam_in <- file.path(bam_dir, "SA_10.sorted.MQ30.primary.bam")
h1n1_bam <- analysis_config$samples[[1]]$bam
h1n1_depth <- file.path(resistance_analysis_dir, "SA_10.depth.txt")
h1n1_vcf <- file.path(resistance_analysis_dir, "SA_10.variants.vcf.gz")
h1n1_cons <- analysis_config$samples[[1]]$consensus

# Skip consensus generation if files already exist
if (file.exists(h1n1_cons) && file.exists(h1n1_bam)) {
  message("✓ Consensus files already exist, skipping generation...")
} else {
  message("Generating consensus with strict filters...")

depth_args_h1n1 <- c(
  depth_script,
  "--bam", h1n1_bam_in, "--out", h1n1_depth, "--vcf-raw", h1n1_vcf, "--consensus", h1n1_cons,
  "--ref", h1n1_ref_fa,
  "--trim-softclips", "--trim-bam", h1n1_bam,
  "--samtools", samtools_cmd,
  "--bcftools", bcftools_cmd,
  "--min-dp", "100",
  "--min-af", "0.9"
)

  exit_code <- system2(python_cmd, depth_args_h1n1)
  if (exit_code != 0) {
    stop("H1N1 consensus generation failed. Check bcftools/samtools in PATH.")
  }
}

#### Mutation Analysis ###
# Run Python script to identify and analyze resistance mutations

message("Running Python analysis script...")

args <- c(
  python_script,
  "--config", json_config_path,
  "--resistance-csv", resistance_csv,
  "--out-csv", output_csv,
  "--samtools", samtools_cmd
)

exit_code <- system2(python_cmd, args)

if (exit_code != 0) {
  stop("Python script failed! Check terminal for errors. Ensure 'pandas' and 'biopython' are installed.")
}

if (!file.exists(output_csv)) {
  stop("Output CSV not generated.")
}

}
}

# End of conditional pipeline execution

#### Data Processing ###
# Process analysis results and export mutation data

message("Data processing complete. Loading results...")

# Export mutations.
plot_data <- read_csv(output_csv, show_col_types = FALSE, na = c(""))

mutations_all <- plot_data %>%
  filter(!is.na(aa_pos), ref_aa != cons_aa, cons_aa != "X")

available_cols <- colnames(mutations_all)
cols_to_keep <- c("sample", "gene", "aa_pos", "ref_aa", "cons_aa", "depth", 
                  "coverage_status", "cons_acc", "ref_acc")
cols_to_keep <- cols_to_keep[cols_to_keep %in% available_cols]

mutations_all <- mutations_all %>%
  select(all_of(cols_to_keep)) %>%
  arrange(sample, gene, aa_pos)

mutations_csv <- file.path(resistance_analysis_dir, "all_mutations.csv")
write_csv(mutations_all, mutations_csv)
message("All mutations exported to: ", mutations_csv)
message("Total mutations found: ", nrow(mutations_all))

#### Plotting Functions ###
# Define plot themes and panel building functions

# Plot helpers.
theme_fancy <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# Panel builder.
make_panel <- function(df, sample_id, genes_ordered, ref_labels) {
  df_sub <- df[df$sample == sample_id, ]

  df_sub$gene <- factor(df_sub$gene, levels = genes_ordered)

  plot_list_ordered <- list()
  gene_refs <- list()

  status_colors <- c(
      "Ref" = "grey85", 
      "Low Coverage" = "#F6E8A6", 
      "Other Mutation" = "grey50", 
      "Resistance Mutation" = "grey50",
      "Stop/Ambig" = "#F6E8A6"
  )
    
  point_colors <- c(
      "Resistance mutation" = "red",
      "Resistance site" = "black"
  )

  status_levels <- names(status_colors)
  point_levels <- names(point_colors)
  legend_extracted <- NULL

  dummy_status <- tibble(
    x = seq_along(status_levels),
    y = 1,
    status = factor(status_levels, levels = status_levels)
  )
  dummy_points <- tibble(
    x = seq_along(point_levels),
    y = 1,
    point_type = factor(point_levels, levels = point_levels)
  )
  p_dummy <- ggplot() +
    geom_tile(data = dummy_status, aes(x = x, y = y, fill = status)) +
    geom_point(data = dummy_points, aes(x = x, y = y, color = point_type), size = 2.5) +
    scale_fill_manual(values = status_colors, guide = guide_legend(nrow = 1, order = 1)) +
    scale_color_manual(values = point_colors, guide = guide_legend(nrow = 1, order = 2)) +
    labs(fill = "Status", color = "Resistance") +
    theme_fancy() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )
  legend_extracted <- p_dummy +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(0, 0, 0, 0)
    )
  message("Legend class: ", class(legend_extracted)[1])

  for (g in genes_ordered) {
    df_g <- df_sub[df_sub$gene == g, ]
    
    if (nrow(df_g) == 0) next
    max_aa <- max(df_g$aa_pos, na.rm = TRUE)
    max_depth <- max(df_g$depth, na.rm = TRUE)
    ref_label <- ref_labels[[paste(sample_id, g, sep = "|")]]
    
    gene_refs[[g]] <- ref_label
    
    p_cov <- ggplot(df_g, aes(x = aa_pos, y = depth)) +
      geom_area(fill = "grey80", alpha = 0.5) +
      geom_line(color = "grey30", linewidth = 0.5) +
      scale_x_continuous(breaks = c(1, max_aa), expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(breaks = c(0, max_depth), limits = c(0, max_depth)) +
      labs(y = "Depth", x = NULL, title = g) +
      coord_cartesian(xlim = c(1, max_aa), clip = "off") +
      theme_fancy() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.length.y = unit(0.15, "cm"),
        plot.margin = margin(1, 0, 0, 0),
        legend.position = "none"
      )

    df_g$yval <- 1
    df_g$point_type <- ifelse(
      df_g$res_mut_found,
      "Resistance mutation",
      ifelse(df_g$is_resistance, "Resistance site", NA)
    )
    
    mut_pos <- df_g$aa_pos[!is.na(df_g$point_type)]
    x_breaks <- unique(sort(c(1, max_aa, mut_pos)))
    
    p_mut <- ggplot(df_g, aes(x = aa_pos, y = yval)) +
      geom_tile(aes(fill = status), height = 0.8) +
      scale_fill_manual(values = status_colors) +
      scale_color_manual(values = point_colors, na.translate = FALSE) +
      scale_x_continuous(breaks = x_breaks, expand = expansion(mult = c(0, 0))) +
      
      geom_point(
        data = df_g[!is.na(df_g$point_type), ],
        aes(color = point_type),
        size = 2.5,
        y = 1
      ) +
      
      labs(x = NULL, y = NULL, fill = "Status", color = "Resistance") +
      coord_cartesian(xlim = c(1, max_aa), ylim = c(0.2, 1.5), clip = "on") +
      theme_fancy()
    
    p_mut <- p_mut +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.length.x = unit(0.15, "cm"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0, 0, 1, 0),
        legend.position = "none"
      )
      
    plot_list_ordered[[paste0(g, "_cov")]] <- p_cov
    plot_list_ordered[[paste0(g, "_mut")]] <- p_mut
  }
  
  aligned_plots <- align_plots(plotlist = plot_list_ordered, align = "v", axis = "lr")
  
  final_gene_rows <- list()
  for (i in seq_along(genes_ordered)) {
     g <- genes_ordered[i]
     
     p_c <- aligned_plots[[ (i-1)*2 + 1 ]]
     p_m <- aligned_plots[[ (i-1)*2 + 2 ]]
     
     pair <- plot_grid(p_c, p_m, ncol = 1, rel_heights = c(0.5, 0.5))
     
     ref_text <- gene_refs[[g]]
     label_panel <- ggdraw() + 
       draw_label(ref_text, angle = -90, size = 10)
     
     row_with_label <- plot_grid(
       pair, 
       label_panel, 
       ncol = 2, 
       rel_widths = c(1, 0.05)
     )
     
     final_gene_rows[[g]] <- row_with_label
  }
  
  combined_genes <- plot_grid(plotlist = final_gene_rows, ncol = 1)
  
  with_title <- plot_grid(
    ggdraw() + draw_label(sample_id, fontface = "bold", size = 16),
    combined_genes,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
  
  return(list(
    main = with_title,
    legend = legend_extracted
  ))
}

#### Figure Generation ###
# Create resistance analysis plots for H1N1 sample

# Make plots.
message("Generating plots...")

ref_labels <- list(
  "H1N1 (SA_10)|HA" = "NC_026433.1",
  "H1N1 (SA_10)|NA" = "NC_026434.1",
  "H1N1 (SA_10)|M2" = "NC_026431.1",
  "H1N1 (SA_10)|PA" = "NC_026437.1"
)

p1 <- make_panel(plot_data, "H1N1 (SA_10)", c("HA", "NA", "M2", "PA"), ref_labels)

main_plot <- p1$main
legend_grob <- p1$legend

main_plot

plot_grid(legend_grob)

final_plot <- plot_grid(
  main_plot,
  legend_grob,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

# Display the final plot
print(final_plot)
message("✓ Figure 4 Part 2 analysis complete! Plot displayed above.")
message("✓ Mutation data saved to: ", mutations_csv)

