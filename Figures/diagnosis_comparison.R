rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})


work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)

gene_plot_data <- read.csv(
  "data/gene_plp_carrier_plot_data.csv",
  check.names = FALSE
)

trio_summary <- read.csv(
  "data/trio_plp_carrier_summary.csv",
  check.names = FALSE
)

sex_summary <- read.csv(
  "data/sex_plp_carrier_summary.csv",
  check.names = FALSE
)

trio_sex_summary <- read.csv(
  "data/trio_sex_plp_carrier_summary.csv",
  check.names = FALSE
)

trio_sex_annotations <- read.csv(
  "data/trio_sex_comparison_annotations.csv",
  check.names = FALSE
)

total_n <- 44029
reference_percent <- 7.7

gene_plot_data$Path_update5 <- factor(
  gene_plot_data$Path_update5,
  levels = c("LP", "P")
)

trio_summary$sfcate <- factor(
  trio_summary$sfcate,
  levels = c("Trio", "Non-trio")
)

sex_summary$sex <- factor(
  sex_summary$sex,
  levels = c("Female", "Male")
)

trio_sex_summary$sfcate <- factor(
  trio_sex_summary$sfcate,
  levels = c("Trio", "Non-trio")
)

trio_sex_summary$sex <- factor(
  trio_sex_summary$sex,
  levels = c("Female", "Male")
)

trio_sex_annotations$sfcate <- factor(
  trio_sex_annotations$sfcate,
  levels = c("Trio", "Non-trio")
)

trio_sex_annotations$x <- as.numeric(trio_sex_annotations$sfcate)

x_labels <- gene_plot_data %>%
  distinct(Gene, x) %>%
  arrange(x)


gene_carrier_plot <- ggplot(gene_plot_data) +

  geom_rect(
    aes(
      xmin = x - 0.32,
      xmax = x + 0.32,
      ymin = ymin,
      ymax = ymax,
      fill = Path_update5
    )
  ) +

  scale_y_continuous(
    name = "Proportion (‰)",
    limits = c(0, 2.05),
    sec.axis = sec_axis(
      ~ . * total_n / 1000,
      name = "Sample count",
      labels = function(x) sprintf("%.0f", x)
    )
  ) +

  scale_x_continuous(
    breaks = x_labels$x,
    labels = x_labels$Gene,
    expand = c(0.01, 0.01)
  ) +

  scale_fill_manual(
    values = c(
      "LP" = "#BBDEFB",
      "P" = "#1E88E5"
    ),
    breaks = c("LP", "P")
  ) +

  labs(
    x = NULL,
    title = NULL,
    fill = NULL
  ) +

  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.88, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.title.y.right = element_text(size = 11, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    )
  )


carrier_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(color = "black"),
  axis.text.x = element_text(size = 11, color = "black"),
  axis.text.y = element_text(size = 10, color = "black"),
  axis.title.y = element_text(size = 11, color = "black"),
  legend.position = "none"
)

make_carrier_comparison_plot <- function(
    data,
    group_col,
    fill_colors,
    y_top,
    y_label = NULL) {

  ggplot(
    data,
    aes(
      x = .data[[group_col]],
      y = prop * 100,
      fill = .data[[group_col]]
    )
  ) +

    geom_col(
      color = "black",
      width = 0.6
    ) +

    geom_errorbar(
      aes(ymin = lower * 100, ymax = upper * 100),
      width = 0.12,
      color = "black"
    ) +

    geom_hline(
      yintercept = reference_percent,
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +

    scale_fill_manual(values = fill_colors) +

    scale_y_continuous(
      limits = c(0, y_top),
      breaks = seq(0, 14, 2),
      expand = expansion(mult = c(0, 0))
    ) +

    labs(
      x = NULL,
      y = y_label,
      title = NULL
    ) +

    carrier_theme +

    annotate(
      "segment",
      x = 1,
      xend = 1,
      y = 10.5,
      yend = 11.1,
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 2,
      xend = 2,
      y = 10.5,
      yend = 11.1,
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 1,
      xend = 2,
      y = 11.1,
      yend = 11.1,
      linewidth = 0.6
    ) +

    annotate(
      "text",
      x = 1.5,
      y = 11.5,
      label = data$p_label[1],
      size = 4
    )
}

trio_carrier_plot <- make_carrier_comparison_plot(
  data = trio_summary,
  group_col = "sfcate",
  fill_colors = c(
    "Trio" = "#FFCC66",
    "Non-trio" = "#64B5F6"
  ),
  y_top = 14.8,
  y_label = "P/LP (%)"
)

sex_carrier_plot <- make_carrier_comparison_plot(
  data = sex_summary,
  group_col = "sex",
  fill_colors = c(
    "Female" = "#FFCC66",
    "Male" = "#64B5F6"
  ),
  y_top = 14.5
)


pdodge <- position_dodge(width = 0.65)

trio_sex_carrier_plot <- ggplot(
  trio_sex_summary,
  aes(x = sfcate, y = prop * 100, fill = sex)
) +

  geom_col(
    color = "black",
    width = 0.5,
    position = pdodge
  ) +

  geom_errorbar(
    aes(ymin = lower * 100, ymax = upper * 100),
    width = 0.12,
    color = "black",
    position = pdodge
  ) +

  geom_hline(
    yintercept = reference_percent,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +

  scale_fill_manual(
    values = c(
      "Female" = "#FFCC66",
      "Male" = "#64B5F6"
    ),
    name = NULL
  ) +

  scale_y_continuous(
    limits = c(0, 15.4),
    breaks = seq(0, 14, 2),
    expand = expansion(mult = c(0, 0))
  ) +

  geom_segment(
    data = trio_sex_annotations,
    aes(
      x = x - 0.18,
      xend = x - 0.18,
      y = y_lo,
      yend = y_hi
    ),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +

  geom_segment(
    data = trio_sex_annotations,
    aes(
      x = x + 0.18,
      xend = x + 0.18,
      y = y_lo,
      yend = y_hi
    ),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +

  geom_segment(
    data = trio_sex_annotations,
    aes(
      x = x - 0.18,
      xend = x + 0.18,
      y = y_hi,
      yend = y_hi
    ),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +

  geom_text(
    data = trio_sex_annotations,
    aes(x = x, y = y_text, label = p_label),
    inherit.aes = FALSE,
    size = 4
  ) +

  labs(x = NULL, y = NULL) +

  carrier_theme +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )


bottom_row <- plot_grid(
  trio_carrier_plot,
  sex_carrier_plot,
  trio_sex_carrier_plot,
  ncol = 3,
  rel_widths = c(1, 1, 2),
  labels = c("C", "D", "E"),
  label_size = 12
)

final_plot <- plot_grid(
  gene_carrier_plot,
  bottom_row,
  ncol = 1,
  labels = c("B", ""),
  label_size = 12
)

print(final_plot)

ggsave(
  filename = "group_comparisons.pdf",
  plot = final_plot,
  width = 7.5,
  height = 6.5,
  units = "in",
  device = "pdf",
  bg = "white",
  limitsize = FALSE
)