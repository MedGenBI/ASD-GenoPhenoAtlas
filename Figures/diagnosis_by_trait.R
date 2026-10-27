rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})


work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)

sig <- read.csv("data/phenotype_group_comparison_results.csv")
plp <- read.csv("data/PLP_results.csv")
raw <- read.csv("data/individual_phenotype_data.csv")

row.names(raw) <- raw$spid

raw$plp <- ifelse(
  raw$spid %in% unique(plp$spid),
  "T",
  "F"
)

sig <- sig %>%
  mutate(
    phenotype_clean = gsub("_", " ", phenotype),
    FDR = as.numeric(FDR)
  )

thresholds <- c(
  "Smiling" = 4,
  "Sitting upright" = 8,
  "Crawling" = 10,
  "Walking" = 15,
  "Spoon feeding self" = 18,
  "Speaking first word" = 18,
  "Speaking with combined word" = 27,
  "Speaking first phrase" = 36,
  "Attaining bladder control" = 42,
  "Attaining bowel control" = 48
)


prepare_panel_data <- function(data, trait_columns) {

  scale_df <- data[, c(trait_columns, "plp"), drop = FALSE]
  vars <- gsub("_", " ", names(scale_df)[-ncol(scale_df)])
  names(scale_df)[-ncol(scale_df)] <- vars

  plp_bool <- toupper(as.character(scale_df$plp)) %in%
    c("T", "TRUE", "1", "YES", "Y")

  list(
    data = scale_df,
    vars = vars,
    plp_bool = plp_bool
  )
}

calculate_bin_stats <- function(v, y, grp, bins, labels = bins) {

  out <- lapply(seq_along(bins), function(i) {

    z <- y[!is.na(grp) & grp == bins[i]]
    n <- length(z)
    carrier <- if (n > 0) sum(z) else 0
    prop <- if (n > 0) carrier / n else NA_real_
    se <- if (n > 0) sqrt(prop * (1 - prop) / n) else NA_real_

    data.frame(
      var = v,
      bin = labels[i],
      n = n,
      carrier = carrier,
      prop = prop,
      se = se,
      lower = if (!is.na(se)) pmax(0, prop - 1.96 * se) else NA_real_,
      upper = if (!is.na(se)) pmin(1, prop + 1.96 * se) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

bin_stats_scales <- function(v, idx, scale_df, plp_bool) {

  x <- scale_df[[v]]
  ok <- !is.na(x) & !is.na(plp_bool)
  x <- x[ok]
  y <- plp_bool[ok]

  if (length(x) == 0) {
    return(NULL)
  }

  qs <- quantile(
    x,
    probs = c(0.25, 0.75),
    type = 7,
    na.rm = TRUE
  )

  bins <- c("0_25", "25_75", "75_100")

  if (qs[1] < qs[2]) {
    grp <- cut(
      x,
      breaks = c(-Inf, qs[1], qs[2], Inf),
      labels = bins,
      include.lowest = TRUE,
      right = TRUE
    )
  } else {
    pr <- rank(x, ties.method = "average") / length(x)

    grp <- cut(
      pr,
      breaks = c(-Inf, 0.25, 0.75, Inf),
      labels = bins,
      include.lowest = TRUE,
      right = TRUE
    )
  }

  severity_labels <- if (idx <= 11) {
    c("mild", "moderate", "severe")
  } else {
    c("severe", "moderate", "mild")
  }

  calculate_bin_stats(
    v = v,
    y = y,
    grp = grp,
    bins = bins,
    labels = severity_labels
  )
}

bin_stats_milestone <- function(v, scale_df, plp_bool) {

  x <- scale_df[[v]]
  ok <- !is.na(x) & !is.na(plp_bool)
  x <- x[ok]
  y <- plp_bool[ok]

  if (length(x) == 0) {
    return(NULL)
  }

  threshold_value <- thresholds[v]
  bins <- c("0_95", "95_100")
  grp <- ifelse(x <= threshold_value, "0_95", "95_100")

  calculate_bin_stats(
    v = v,
    y = y,
    grp = grp,
    bins = bins
  )
}

bin_stats_binary <- function(v, scale_df, plp_bool) {

  x <- scale_df[[v]]
  ok <- !is.na(x) & !is.na(plp_bool)
  x <- x[ok]
  y <- plp_bool[ok]

  if (length(x) == 0) {
    return(NULL)
  }

  bins <- c("Absent (0)", "Present (1)")

  grp <- factor(
    x,
    levels = c(0, 1),
    labels = bins
  )

  calculate_bin_stats(
    v = v,
    y = y,
    grp = grp,
    bins = bins
  )
}

to_stars <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
}


make_bar_plot <- function(
    plot_df,
    fill_colors,
    y_top,
    errorbar_width,
    legend_title = NULL) {

  ggplot(plot_df, aes(x = var, y = prop * 100, fill = bin)) +

    geom_col(
      position = position_dodge(width = 0.78),
      width = 0.58,
      color = "black"
    ) +

    geom_errorbar(
      aes(ymin = lower * 100, ymax = upper * 100),
      position = position_dodge(width = 0.78),
      width = errorbar_width,
      linewidth = 0.4,
      color = "black"
    ) +

    geom_hline(
      yintercept = 8,
      linetype = "dashed",
      color = "black",
      linewidth = 0.4
    ) +

    scale_fill_manual(
      values = fill_colors,
      name = legend_title
    ) +

    scale_y_continuous(
      limits = c(0, y_top),
      breaks = seq(0, y_top, 2),
      expand = expansion(mult = c(0, 0))
    ) +

    labs(x = NULL, y = "", title = NULL) +

    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(
        angle = 30,
        hjust = 1,
        size = 9,
        color = "black"
      ),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 11, color = "black"),
      legend.position = "none",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

add_significance <- function(
    p,
    plot_df,
    sig_data,
    y_top,
    scales_panel = FALSE,
    bracket_width = 0.25,
    text_padding = 0.9) {

  top_df <- aggregate(
    upper ~ var,
    data = plot_df,
    FUN = max,
    na.rm = TRUE
  )

  top_df$y_base <- top_df$upper * 100

  sig_data <- sig_data %>%
    filter(
      phenotype_clean %in% levels(plot_df$var),
      !is.na(FDR),
      FDR < 0.05
    )

  ann <- merge(
    top_df,
    sig_data,
    by.x = "var",
    by.y = "phenotype_clean",
    all = FALSE
  )

  if (nrow(ann) == 0) {
    return(p)
  }

  ann$xi <- as.numeric(
    factor(ann$var, levels = levels(plot_df$var))
  )

  if (scales_panel) {

    ann <- ann %>%
      group_by(var, baseline, control) %>%
      mutate(
        comp_id = paste(baseline, control, sep = "_")
      ) %>%
      ungroup() %>%
      group_by(var) %>%
      mutate(
        comp_index = match(comp_id, unique(comp_id)),
        height_offset = 0.8 * comp_index
      ) %>%
      ungroup()

    ann$y_lo <- ann$y_base + 0.3 + ann$height_offset
    ann$y_hi <- ann$y_base + 0.6 + ann$height_offset
    ann$y_text <- ann$y_base + 0.7 + ann$height_offset

    ann$x_left <- ann$xi - 0.3
    ann$x_right <- ann$xi + ifelse(
      ann$control == "moderate",
      0.05,
      0.3
    )

  } else {

    ann$y_lo <- ann$y_base + 0.3
    ann$y_hi <- ann$y_base + 0.6
    ann$y_text <- ann$y_base + text_padding

    ann$x_left <- ann$xi - bracket_width
    ann$x_right <- ann$xi + bracket_width
  }

  ann$label <- to_stars(ann$FDR)

  seg_df <- rbind(
    data.frame(
      x = ann$x_left,
      xend = ann$x_left,
      y = ann$y_lo,
      yend = ann$y_hi
    ),
    data.frame(
      x = ann$x_right,
      xend = ann$x_right,
      y = ann$y_lo,
      yend = ann$y_hi
    ),
    data.frame(
      x = ann$x_left,
      xend = ann$x_right,
      y = ann$y_hi,
      yend = ann$y_hi
    )
  )

  new_top <- max(y_top, max(ann$y_text + 0.25, na.rm = TRUE))
  break_top <- ceiling(new_top / 2) * 2

  p +

    geom_segment(
      data = seg_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +

    geom_text(
      data = ann,
      aes(x = xi, y = y_text, label = label),
      inherit.aes = FALSE,
      size = 4
    ) +

    scale_y_continuous(
      limits = c(0, new_top),
      breaks = seq(0, break_top, by = 2),
      expand = expansion(mult = c(0, 0))
    )
}


panel_data <- prepare_panel_data(
  data = raw,
  trait_columns = names(raw)[13:32]
)

plot_df <- bind_rows(
  lapply(seq_along(panel_data$vars), function(i) {
    bin_stats_scales(
      v = panel_data$vars[i],
      idx = i,
      scale_df = panel_data$data,
      plp_bool = panel_data$plp_bool
    )
  })
)

plot_df$var <- factor(
  plot_df$var,
  levels = panel_data$vars
)

plot_df$bin <- factor(
  plot_df$bin,
  levels = c("mild", "moderate", "severe")
)

fill_colors <- c(
  "mild" = "#BBDEFB",
  "moderate" = "#64B5F6",
  "severe" = "#1E88E5"
)

p1 <- make_bar_plot(
  plot_df = plot_df,
  fill_colors = fill_colors,
  y_top = 16,
  errorbar_width = 0.45,
  legend_title = "Severity"
)

p1 <- add_significance(
  p = p1,
  plot_df = plot_df,
  sig_data = sig,
  y_top = 16,
  scales_panel = TRUE
)


panel_data <- prepare_panel_data(
  data = raw,
  trait_columns = names(raw)[33:42]
)

plot_df <- bind_rows(
  lapply(panel_data$vars, function(v) {
    bin_stats_milestone(
      v = v,
      scale_df = panel_data$data,
      plp_bool = panel_data$plp_bool
    )
  })
)

plot_df$var <- factor(
  plot_df$var,
  levels = panel_data$vars
)

plot_df$bin <- factor(
  plot_df$bin,
  levels = c("0_95", "95_100")
)

fill_colors <- c(
  "0_95" = "#FCE4D6",
  "95_100" = "#F8C9A8"
)

p2 <- make_bar_plot(
  plot_df = plot_df,
  fill_colors = fill_colors,
  y_top = 18,
  errorbar_width = 0.4
)

p2 <- add_significance(
  p = p2,
  plot_df = plot_df,
  sig_data = sig,
  y_top = 18,
  bracket_width = 0.22
)


make_binary_panel <- function(
    trait_columns,
    fill_colors,
    y_top,
    errorbar_width,
    text_padding = 0.9) {

  panel_data <- prepare_panel_data(
    data = raw,
    trait_columns = trait_columns
  )

  plot_df <- bind_rows(
    lapply(panel_data$vars, function(v) {
      bin_stats_binary(
        v = v,
        scale_df = panel_data$data,
        plp_bool = panel_data$plp_bool
      )
    })
  )

  plot_df$var <- factor(
    plot_df$var,
    levels = panel_data$vars
  )

  plot_df$bin <- factor(
    plot_df$bin,
    levels = c("Absent (0)", "Present (1)")
  )

  p <- make_bar_plot(
    plot_df = plot_df,
    fill_colors = fill_colors,
    y_top = y_top,
    errorbar_width = errorbar_width
  )

  add_significance(
    p = p,
    plot_df = plot_df,
    sig_data = sig,
    y_top = y_top,
    bracket_width = 0.25,
    text_padding = text_padding
  )
}

p3 <- make_binary_panel(
  trait_columns = names(raw)[43:61],
  fill_colors = c(
    "Absent (0)" = "#B2DFDB",
    "Present (1)" = "#80CBC4"
  ),
  y_top = 14,
  errorbar_width = 0.5
)

p4 <- make_binary_panel(
  trait_columns = names(raw)[62:71],
  fill_colors = c(
    "Absent (0)" = "#FFF9C4",
    "Present (1)" = "#FFF176"
  ),
  y_top = 18,
  errorbar_width = 0.45
)

p5 <- make_binary_panel(
  trait_columns = names(raw)[72:77],
  fill_colors = c(
    "Absent (0)" = "#d9e1f2",
    "Present (1)" = "#b3c4e8"
  ),
  y_top = 24,
  errorbar_width = 0.3,
  text_padding = 1.2
)

p6 <- make_binary_panel(
  trait_columns = names(raw)[78:84],
  fill_colors = c(
    "Absent (0)" = "#dcedc8",
    "Present (1)" = "#aed581"
  ),
  y_top = 24,
  errorbar_width = 0.3,
  text_padding = 1.2
)

p7 <- make_binary_panel(
  trait_columns = names(raw)[85:90],
  fill_colors = c(
    "Absent (0)" = "#f8bbd0",
    "Present (1)" = "#f48fb1"
  ),
  y_top = 24,
  errorbar_width = 0.3
)


all_aligned <- align_plots(
  p2, p4, p5, p6, p7,
  align = "v",
  axis = "l"
)

row1 <- plot_grid(
  p1,
  nrow = 1,
  labels = "A",
  label_size = 12,
  label_y = 0.98,
  vjust = 1
)

row2 <- plot_grid(
  all_aligned[[1]],
  all_aligned[[2]],
  nrow = 1,
  rel_widths = c(1, 1),
  labels = c("B", "D"),
  label_size = 12,
  label_y = 0.98,
  vjust = 1
)

row3 <- plot_grid(
  p3,
  nrow = 1,
  labels = "C",
  label_size = 12,
  label_y = 0.98,
  vjust = 1
)

row4 <- plot_grid(
  all_aligned[[3]],
  all_aligned[[4]],
  all_aligned[[5]],
  nrow = 1,
  rel_widths = c(1, 1, 1),
  labels = c("E", "F", "G"),
  label_size = 12,
  label_y = 0.98,
  vjust = 1
)

combined <- plot_grid(
  row1, row2, row3, row4,
  ncol = 1,
  rel_heights = c(1, 1, 1, 1)
)

final_plot <- ggdraw() +
  draw_plot(
    combined,
    x = 0.055,
    y = 0,
    width = 0.89,
    height = 1
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

print(final_plot)

ggsave(
  filename = "final.pdf",
  plot = final_plot,
  width = 210,
  height = 297,
  units = "mm",
  device = "pdf",
  limitsize = FALSE
)