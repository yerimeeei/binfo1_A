
# ==================================== #
# Visualization: Ribosome Density Plot #
# ==================================== #
library(ggplot2)
library(gridExtra)
library(zoo)
library(dplyr)

plot_ribosome_density <- function(siLuc_df, siLin28a_df, output_dir = "analysis_results/figures") {

  calculate_stats <- function(df) {
    mean_vals <- colMeans(df, na.rm = TRUE)
    sem_vals <- apply(df, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
    list(mean = mean_vals, sem = sem_vals)
  }

  stats_siLuc <- calculate_stats(siLuc_df)
  stats_siLin28a <- calculate_stats(siLin28a_df)

  positions <- as.integer(colnames(siLuc_df))

  p_values <- mapply(function(x, y) {
    t.test(x, y, var.equal = FALSE)$p.value
  }, as.data.frame(siLuc_df), as.data.frame(siLin28a_df))

  p_values_adj <- p.adjust(p_values, method = "BH")

  data_plot <- data.frame(
    Position = positions,
    siLuc_mean = stats_siLuc$mean,
    siLuc_lower = stats_siLuc$mean - stats_siLuc$sem,
    siLuc_upper = stats_siLuc$mean + stats_siLuc$sem,
    siLin28a_mean = stats_siLin28a$mean,
    siLin28a_lower = stats_siLin28a$mean - stats_siLin28a$sem,
    siLin28a_upper = stats_siLin28a$mean + stats_siLin28a$sem,
    neglog10_p = -log10(p_values_adj)
  )

  window <- 5
  data_plot$rolling_avg <- rollmean(data_plot$neglog10_p, k = window, fill = NA, align = "center")

  y_max <- max(c(data_plot$siLuc_mean, data_plot$siLin28a_mean), na.rm = TRUE)
  sig_positions <- data_plot$Position[p_values_adj < 0.05]

  # Plot 1: Ribosome density
  p1 <- ggplot(data_plot, aes(x = Position)) +
    geom_line(aes(y = siLuc_mean, color = "siLuc (control)"), size = 1.2) +
    geom_ribbon(aes(ymin = siLuc_lower, ymax = siLuc_upper), fill = "#1f77b4", alpha = 0.2) +
    geom_line(aes(y = siLin28a_mean, color = "siLin28a"), size = 1.2) +
    geom_ribbon(aes(ymin = siLin28a_lower, ymax = siLin28a_upper), fill = "#ff7f0e", alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_segment(data = data.frame(x = sig_positions), aes(x = x, xend = x, y = y_max * 1.03, yend = y_max * 1.06),
                 color = "red", size = 0.3, inherit.aes = FALSE) +
    scale_color_manual(values = c("siLuc (control)" = "#1f77b4", "siLin28a" = "#ff7f0e")) +
    labs(y = "Mean ribosome density\n(reads per million)",
         title = "Ribosome Density Around Start Codon") +
    theme_minimal(base_size = 13) +
    theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))

  # Plot 2: -log10(p-value) and running average
  p2 <- ggplot(data_plot, aes(x = Position)) +
    geom_line(aes(y = neglog10_p), color = "#2ca02c", size = 1) +
    geom_line(aes(y = rolling_avg), color = "#9467bd", size = 1.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
    labs(x = "Position relative to start codon (nt)", y = "-log10(p-value)") +
    theme_minimal(base_size = 13)

  # Combine plots
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_path <- file.path(output_dir, "ribosome_density_profiles.png")

  ggsave(output_path, grid.arrange(p1, p2, heights = c(3, 1)), width = 12, height = 10, dpi = 300)
  message(sprintf("Saved ribosome density plot to %s", output_path))

  return(list(p1 = p1, p2 = p2))
}

siLuc_df <- read.csv("/Users/rimi/Projects/binfo1_A/binfo1_export/siLuc_start_density.csv", row.names = 1, check.names = FALSE)
siLin28a_df <- read.csv("/Users/rimi/Projects/binfo1_A/binfo1_export/siLin28a_start_density.csv", row.names = 1, check.names = FALSE)

plots <- plot_ribosome_density(siLuc_df, siLin28a_df)
grid.arrange(plots$p1, plots$p2, heights = c(3, 1))
plot_ribosome_density(siLuc_df, siLin28a_df, output_dir = "analysis_results/figures")

# =========================== #
# Heatmap of Ribosome Density #
# =========================== #
library(pheatmap)
library(RColorBrewer)

luc <- read.csv("/Users/rimi/Projects/binfo1_A/binfo1_export/siLuc_start_density.csv", row.names = 1, check.names = FALSE)
lin <- read.csv("/Users/rimi/Projects/binfo1_A/binfo1_export/siLin28a_start_density.csv", row.names = 1, check.names = FALSE)

common_ids <- intersect(rownames(luc), rownames(lin))
luc <- luc[common_ids, ]
lin <- lin[common_ids, ]

luc <- luc[rowSums(luc) > 50, ]
lin <- lin[rownames(luc), ]

luc_z <- t(scale(t(luc)))
lin_z <- t(scale(t(lin)))

ordering <- order(apply(luc_z, 1, function(x) which.max(x)))
luc_sorted <- luc_z[ordering, ]
lin_sorted <- lin_z[ordering, ]

combined <- cbind(luc_sorted, lin_sorted)
colnames(combined) <- make.unique(rep(colnames(luc), 2))

heat <- pheatmap(combined,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Start Codon Ribosome Density: siLuc vs siLin28a",
         fontsize = 10)
ggsave("analysis_results/figures/heatmap_final.jpeg",
       plot = heat, width = 10, height = 8, dpi = 300)

# ================================ #
# Volcano: Start Codon Pause Score #
# ================================ #
common_ids <- intersect(rownames(luc), rownames(lin))
luc <- luc[common_ids, ]
lin <- lin[common_ids, ]

pause_window <- as.character(-3:9)
luc_pause <- rowSums(luc[, pause_window], na.rm = TRUE)
lin_pause <- rowSums(lin[, pause_window], na.rm = TRUE)

pause_df <- data.frame(
  gene = common_ids,
  luc = luc_pause,
  lin = lin_pause
) %>%
  mutate(
    log2FC = log2((lin + 1) / (luc + 1)),  # pseudocount
    diff = lin - luc
  )

vol <- ggplot(pause_df, aes(x = log2FC, y = abs(diff))) +
  geom_point(alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(title = "Ribosome Pause Score Change",
       x = "log2 Fold Change (siLin28a / siLuc)",
       y = "Absolute Pause Difference") +
  theme_minimal()
ggsave("analysis_results/figures/volcano_final.jpeg",
       plot = vol, width = 10, height = 8, dpi = 300)
