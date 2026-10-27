rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})


work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)

palette_main <- c(
  "Without" = "#eaeaea",
  "LP" = "#c2dff1",
  "P" = "#4386dd"
)

palette_variant <- c(
  "Stopgain" = "#64B5F6",
  "Nonsynonymous SNV" = "#80CBC4",
  "Frameshift deletion" = "#81C784",
  "Frameshift insertion" = "#FFF59D",
  "Splicing" = "#FFCC80",
  "Other" = "#A1887F"
)

palette_inheritance <- c(
  "AD" = "#f7ce76",
  "XD" = "#BBDEFB",
  "AR" = "#c3bef0",
  "XR" = "#A5D6A7"
)

df_main <- data.frame(
  category = c("Without", "LP", "P"),
  percentage = c(92.3, 5.7, 2.0)
) %>%
  mutate(
    category = factor(category, levels = names(palette_main))
  )

df_variant <- data.frame(
  category = c(
    "Stopgain",
    "Nonsynonymous SNV",
    "Frameshift deletion",
    "Frameshift insertion",
    "Splicing",
    "Other"
  ),
  percentage = c(29.0, 25.9, 19.9, 14.9, 9.8, 0.5)
) %>%
  mutate(
    category = factor(category, levels = names(palette_variant))
  )

df_inheritance <- data.frame(
  category = c("AD", "XD", "XR", "AR"),
  percentage = c(88.4, 5.3, 3.5, 2.8)
) %>%
  mutate(
    category = factor(category, levels = names(palette_inheritance))
  )


create_pie <- function(
    data,
    colors,
    title = NULL,
    show_legend = TRUE,
    highlight_cats = NULL) {

  plot_data <- data %>%
    mutate(
      border_width = ifelse(category %in% highlight_cats, 0.5, 0.25)
    )

  p <- ggplot(
    plot_data,
    aes(x = "", y = percentage, fill = category)
  ) +

    geom_col(
      aes(linewidth = border_width),
      width = 1,
      color = "black"
    ) +

    scale_linewidth_identity() +

    coord_polar(
      theta = "y",
      start = pi / 2.5
    ) +

    scale_fill_manual(values = colors) +

    labs(
      x = NULL,
      y = NULL,
      title = title,
      fill = NULL
    ) +

    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  p
}

pie_main <- create_pie(
  data = df_main,
  colors = palette_main,
  highlight_cats = c("LP", "P")
)

pie_variant <- create_pie(
  data = df_variant,
  colors = palette_variant
)

pie_inheritance <- create_pie(
  data = df_inheritance,
  colors = palette_inheritance
)

print(pie_main)
print(pie_variant)
print(pie_inheritance)

ggsave(
  filename = "variant_pathogenicity_distribution.pdf",
  plot = pie_main,
  width = 5,
  height = 4,
  units = "in",
  device = "pdf",
  limitsize = FALSE
)

ggsave(
  filename = "variant_functional_category_distribution.pdf",
  plot = pie_variant,
  width = 6.5,
  height = 4,
  units = "in",
  device = "pdf",
  limitsize = FALSE
)

ggsave(
  filename = "inheritance_mode_distribution.pdf",
  plot = pie_inheritance,
  width = 5,
  height = 4,
  units = "in",
  device = "pdf",
  limitsize = FALSE
)