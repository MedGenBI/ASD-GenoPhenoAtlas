# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)


work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)


############################# 1. Dual-axis plots ##############################

# plot.type <- 'fig4a'
plot.type <- 'fig4b'

data_path <- '../data/gene_level_evaluation_metrics.tsv'
data <- read.csv(data_path, check.names = F, sep = '\t')

if (plot.type == 'fig4a') {
  plot.data <- data[data$fdr_phenotypes > 0, ]
  pdf_out <- 'output/Fig4a.pdf'
} else if (plot.type == 'fig4b') {
  plot.data <- data[data$fdr_phenotypes == 0 & data$auc_median >= 0.8, ]
  pdf_out <- 'output/Fig4b.pdf'
}

cat("n row:", nrow(plot.data), "\n")

or_cap_value <- 50  

plot.data <- plot.data %>%
  mutate(
    is_capped = or > or_cap_value | is.infinite(or),
    or_plot = ifelse(or > or_cap_value | is.infinite(or), or_cap_value, or),
    auc_group = case_when(
      auc_median >= 0.8 ~ "AUC >= 0.8",
      auc_median >= 0.7 ~ "0.7 <= AUC < 0.8",
      TRUE              ~ "AUC < 0.7"
    ),
    auc_group = factor(auc_group, levels = c("AUC >= 0.8", "0.7 <= AUC < 0.8", "AUC < 0.7"))
  ) %>%
  mutate(gene = fct_reorder(gene, desc(auc_median)))

primary_max <- 1.0
secondary_max <- max(plot.data$p_phenotypes) * 1.2

scaling_factor <- secondary_max / primary_max

dual_axis_plot <- ggplot(plot.data, aes(x = gene)) +

  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = 0.8,
    linetype = "dashed",
    color = "#F09C2B",
    linewidth = 0.4
  ) +
  geom_col(
    aes(y = or_plot  / scaling_factor), 
    fill = "#a1c9f4",
    width = 0.7,
    alpha = 0.8
  ) +
  geom_text(
    data = subset(plot.data, is_capped == TRUE),
    aes(y = or_plot / scaling_factor, label = "^"),
    vjust = .3,
    size = 3,
    fontface = "bold",
    color = "#2a6ebb"
  ) +
  geom_errorbar(
    aes(ymin = auc_ci_lower, ymax = auc_ci_upper, color = auc_group),
    linewidth = 0.4,
    width = 0.3
  ) +
  geom_point(
    aes(y = auc_median, fill = auc_group, color = auc_group),
    shape = 23,
    size = 1.5,
    stroke = .6
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "AUC >= 0.8"        = "#dd666e",
      "0.7 <= AUC < 0.8"  = "#e79a36",
      "AUC < 0.7"         = "#64B5F6"
    )
  ) +
  scale_fill_manual(
    name = "",
    values = c(
      "AUC >= 0.8"        = "#FCA09D",
      "0.7 <= AUC < 0.8"  = "#F8C9A8",
      "AUC < 0.7"         = "#BBDEFB"
    )
  ) +
  scale_y_continuous(
    name = "Median AUC (95% CI)",
    limits = c(0, primary_max), expand = c(0, 0),
    sec.axis = sec_axis(
      trans = ~ . * scaling_factor, 
      name = "Odds Ratio",
      breaks = seq(0, or_cap_value, by = 10),
      labels = function(x) ifelse(x == or_cap_value, paste0(">=", x), x)
    )
  ) +
  scale_x_discrete(name = "Gene") +
  theme_classic(base_size = 7) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.line.x = element_line(linewidth = .2),
    axis.line.y.right = element_line(color = "black", linewidth = .2),
    axis.ticks.y.right = element_line(color = "black"),
    axis.text.y.right = element_text(color = "black", size = 7),
    axis.title.y.right = element_text(color = "black", size = 7),
    axis.line.y.left = element_line(color = "black", linewidth = .2),
    axis.ticks.y.left = element_line(color = "black"),
    axis.text.y.left = element_text(color = "black", size = 7),
    axis.title.y.left = element_text(color = "black", size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    legend.position = "top"
  ) +
  labs(
    title = "GPS Performance per Gene"
  )

# Save
ggsave(
  pdf_out, dual_axis_plot, 
  width = 6.69, height = 3, dpi = 300, units = "in"
)
################################################################################


############################### 2. Forest Plot #################################

plot.type <- 'wes'
# plot.type <- 'wgs'

if (plot.type == 'wes') {
  data_path <- file.path(work_dir, '../data/fig4.c.wes.forest.tsv')
  pdf_out <- 'output/Fig4c.forest.pdf'
} else if (plot.type == 'wgs') {
  data_path <- file.path(work_dir, '../data/fig4.d.wgs.forest.tsv')
  pdf_out <- 'output/Fig4d.forest.pdf'
}


data <- read.csv(data_path, check.names = F, sep = '\t')

top_p <- c("20.0%","10.0%","5.0%","1.0%", "0.1%")

df <- data %>%
  filter(top_percent %in% top_p) %>%
  mutate(
    top_percent = factor(top_percent, levels = top_p)
  ) %>%
  group_by(group) %>%
  mutate(
    delta = or[top_percent == "0.1%"] - or[top_percent == "20.0%"],
  ) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("High-confidence","Extended")))

df <- df %>%
  mutate(
    signif_group = case_when(
      p < 0.05   ~ "p<0.05",
      TRUE       ~ "NS"
    ),
    signif_group = factor(signif_group, levels = c("p<0.05", "NS")),
    fill_group = ifelse(p < 0.05, as.character(group), NA_character_),
    fill_group = factor(fill_group, levels = levels(group))
  )

df <- df %>%
  filter(or > 0, or_ci_low > 0, or_ci_upp > 0)

group_border_colors <- c(
  "High-confidence" = "#DD666E",      # Hihg-performance Set
  "Extended" = "#64B5F6"              # Extended Set
)

group_fill_colors <- c(
  "High-confidence" = "#FCA09D",      # Hihg-performance Set
  "Extended" = "#BBDEFB"              # Extended Set
)

neutral_line <- "#DDDDDD" 
n_col = 1

p <- ggplot(df, aes(x = top_percent, y = or, group = group )) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey70", linewidth = .3) +
  geom_line(aes(color = group), linewidth = .3) +
  geom_errorbar(aes(ymin = or_ci_low, ymax = or_ci_upp, color = group),
                width = 0.1, linewidth = 0.3) +
  geom_point(
    aes(color = group, fill  = fill_group),
    shape = 23, size = 1, stroke = .4
  ) +
  scale_color_manual(
    values = group_border_colors, 
    name = "group"
  ) +
  scale_fill_manual(
    values = group_fill_colors, 
    name = "group", 
    na.value = "white"
  ) + 
  facet_wrap(~ group, scales = "fixed", ncol = n_col) +
  scale_y_log10() +
  labs(
    x = "Percentile Cutoff",
    y = "Odds Ratio"
  ) +
  theme_minimal(base_size = 6) +
  theme(
    text = element_text(size = 6),
    strip.text =element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    panel.spacing = unit(0.6, "lines"),
    plot.background = element_rect(fill = "white", color = "#D1D1D1", linewidth = .5)
  ) 


# Save
ggsave(pdf_out, p, width = 1.2, height = 1.4, dpi = 300, units = "in", device = cairo_pdf)
################################################################################
