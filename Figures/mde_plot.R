rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(patchwork)
})


work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)

dat <- read_excel("data/circle.xlsx")

cate_order <- c(
  "Scales",
  "Developmental milestones",
  "Mental health and psychiatric conditions",
  "Neurodevelopmental and related conditions",
  "Perinatal complications",
  "Birth defects",
  "Growth abnormalities"
)

cate_colors <- c(
  "Scales"                                  = "#c1d3f4",
  "Developmental milestones"                = "#FEE796",
  "Mental health and psychiatric conditions" = "#F9CBAA",
  "Neurodevelopmental and related conditions" = "#d2deb8",
  "Perinatal complications"                 = "#9bc6dd",
  "Birth defects"                           = "#deccf2",
  "Growth abnormalities"                    = "#f3bcce"
)

dat_plot_base <- dat %>%
  select(cate2, count, `Effect ratio`, `Effect ratio (genpwr)`) %>%
  filter(
    !is.na(cate2),
    !is.na(count),
    !is.na(`Effect ratio`),
    !is.na(`Effect ratio (genpwr)`)
  ) %>%
  mutate(
    cate2 = factor(cate2, levels = cate_order),
    count = as.numeric(count),
    `Effect ratio` = as.numeric(`Effect ratio`),
    `Effect ratio (genpwr)` = as.numeric(`Effect ratio (genpwr)`)
  ) %>%
  arrange(cate2, count)

gap_size <- 3

dat_plot_base <- dat_plot_base %>%
  mutate(
    cate_id = as.integer(cate2),
    plot_order = row_number(),
    theta_pos = plot_order + (cate_id - 1) * gap_size
  )

total_theta <- nrow(dat_plot_base) + length(cate_order) * gap_size

breaks_real <- c(0, 1, 1.2, 1.4, 2, 5)
breaks_plot <- c(0, 1, 2, 3, 4, 5)

ref_lines <- data.frame(
  y = c(1, 2, 3, 4, 5),
  label = factor(
    c("1", "1.2", "1.4", "2", "5"),
    levels = c("1", "1.2", "1.4", "2", "5")
  ),
  color = c("black", "black", "black", "black", "grey75")
)


make_polar_plot <- function(data, value_col) {

  plot_data <- data %>%
    mutate(
      radial_pos = approx(
        x = breaks_real,
        y = breaks_plot,
        xout = .data[[value_col]],
        rule = 2
      )$y
    )

  ref_black <- ref_lines %>%
    filter(color == "black")

  ref_grey <- ref_lines %>%
    filter(color == "grey75")

  ggplot(plot_data, aes(x = theta_pos, y = radial_pos, color = cate2)) +

    geom_hline(
      data = ref_black,
      aes(yintercept = y, linetype = label),
      linewidth = 0.6,
      color = "black"
    ) +

    geom_hline(
      data = ref_grey,
      aes(yintercept = y, linetype = label),
      linewidth = 0.45,
      color = "grey75"
    ) +

    geom_point(
      shape = 16,
      size = 2.8,
      alpha = 0.95
    ) +

    coord_polar(
      theta = "x",
      start = 0,
      direction = 1,
      clip = "off"
    ) +

    scale_x_continuous(
      limits = c(0.5, total_theta + 0.5),
      breaks = NULL,
      expand = c(0, 0)
    ) +

    scale_y_continuous(
      limits = c(0, 5),
      breaks = breaks_plot,
      labels = c("0", "1", "1.2", "1.4", "2", "5"),
      expand = c(0, 0)
    ) +

    scale_color_manual(
      values = cate_colors,
      drop = FALSE
    ) +

    scale_linetype_manual(
      values = c(
        "1" = "solid",
        "1.2" = "dashed",
        "1.4" = "dotted",
        "2" = "dotdash",
        "5" = "solid"
      ),
      name = NULL
    ) +

    guides(
      linetype = guide_legend(
        override.aes = list(
          color = c("black", "black", "black", "black", "grey75")
        )
      )
    ) +

    labs(
      x = NULL,
      y = "Observed effect / MDE",
      color = NULL,
      linetype = NULL
    ) +

    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10, color = "grey20"),
      axis.title.y = element_text(size = 11),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.key.height = grid::unit(0.55, "cm"),
      plot.margin = margin(20, 35, 20, 20, "pt")
    )
}

p1 <- make_polar_plot(dat_plot_base, "Effect ratio")
p2 <- make_polar_plot(dat_plot_base, "Effect ratio (genpwr)")

combined <- p1 + p2 +
  plot_layout(guides = "collect")

print(combined)

ggsave(
  filename = "circle.pdf",
  plot = combined,
  width = 10,
  height = 5,
  units = "in",
  limitsize = FALSE
)