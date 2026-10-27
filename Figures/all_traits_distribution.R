rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(ggbreak)
})


work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)

phenotype_data_path <- "data/clinical_phenotype_data.xlsx"

scale_increasing_severity <- c(
  "SCQ_total",
  "SCQ_communication",
  "SCQ_interaction",
  "SCQ_stereotype",
  "RBSR_total",
  "RBSR_stereotyped",
  "RBSR_injurious",
  "RBSR_compulsive",
  "RBSR_ritualistic",
  "RBSR_sameness",
  "RBSR_restricted"
)

scale_decreasing_severity <- c(
  "DCDQ_total",
  "DCDQ_movement",
  "DCDQ_handwriting",
  "DCDQ_coordination",
  "ABC_total",
  "ABC_communication",
  "ABC_daily_living_skills",
  "ABC_socialization",
  "ABC_motor_skills"
)

milestone_thresholds <- c(
  "Smiling" = 4,
  "Sitting_upright" = 8,
  "Crawling" = 10,
  "Walking" = 15,
  "Spoon_feeding_self" = 18,
  "Speaking_first_word" = 18,
  "Speaking_with_combined_word" = 27,
  "Speaking_first_phrase" = 36,
  "Attaining_bladder_control" = 42,
  "Attaining_bowel_control" = 48
)


read_phenotype_sheet <- function(sheet_name) {

  data <- read.xlsx(
    phenotype_data_path,
    sheet = sheet_name,
    check.names = FALSE
  )

  data %>%
    select(-any_of("spid"))
}

make_percentage_plot <- function(
    plot_data,
    variable_order,
    fill_colors,
    title,
    legend_title,
    border_color,
    use_axis_break = FALSE,
    right_margin_cm = 1) {

  p <- ggplot(
    plot_data,
    aes(
      x = factor(Variable, levels = variable_order),
      y = Percentage,
      fill = Group
    )
  ) +

    geom_col(
      position = "stack",
      width = 0.7,
      color = border_color,
      linewidth = 0.3
    ) +

    scale_fill_manual(
      values = fill_colors,
      name = legend_title
    ) +

    scale_y_continuous(
      limits = c(0, 100.001),
      breaks = c(0, 25, 50, 75, 100),
      labels = c("0", "25", "50", "75", "100"),
      expand = expansion(mult = c(0, 0.05))
    ) +

    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +

    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 50,
        hjust = 1,
        size = 12,
        color = "black"
      ),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(
        size = 12,
        face = "bold",
        margin = margin(t = 15)
      ),
      axis.title.y = element_text(size = 12, face = "bold"),
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, right_margin_cm, 1.5, 1.5, "cm"),
      legend.position = "none",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = grid::unit(0.8, "lines")
    )

  if (use_axis_break) {
    p <- p +

      scale_y_break(
        breaks = c(1, 50),
        scales = 17,
        space = 0.05,
        expand = TRUE
      ) +

      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.line.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
      )
  }

  p
}

save_percentage_plot <- function(plot, filename, n_variables) {

  plot_width <- 4 + n_variables * 0.25 * 1.75

  print(plot)

  ggsave(
    filename = filename,
    plot = plot,
    width = plot_width,
    height = 8,
    units = "in",
    device = "pdf",
    bg = "white",
    limitsize = FALSE
  )
}


scale_data <- read_phenotype_sheet("Scales")
scale_cols <- colnames(scale_data)

unexpected_scales <- setdiff(
  scale_cols,
  c(scale_increasing_severity, scale_decreasing_severity)
)

if (length(unexpected_scales) > 0) {
  stop(
    "Severity direction is undefined for: ",
    paste(unexpected_scales, collapse = ", ")
  )
}

scale_summary <- bind_rows(
  lapply(scale_cols, function(col_name) {

    col_data <- na.omit(scale_data[[col_name]])
    total_count <- length(col_data)

    if (total_count == 0) {
      return(
        data.frame(
          Variable = col_name,
          mild_pct = NA_real_,
          moderate_pct = NA_real_,
          severe_pct = NA_real_,
          total_count = 0
        )
      )
    }

    q <- quantile(
      col_data,
      probs = c(0.25, 0.75),
      na.rm = TRUE
    )

    low_pct <- mean(col_data <= q[1]) * 100
    medium_pct <- mean(col_data > q[1] & col_data <= q[2]) * 100
    high_pct <- mean(col_data > q[2]) * 100

    increasing_severity <- col_name %in% scale_increasing_severity

    data.frame(
      Variable = col_name,
      mild_pct = if (increasing_severity) low_pct else high_pct,
      moderate_pct = medium_pct,
      severe_pct = if (increasing_severity) high_pct else low_pct,
      total_count = total_count
    )
  })
)

scale_plot_data <- scale_summary %>%
  pivot_longer(
    cols = c(mild_pct, moderate_pct, severe_pct),
    names_to = "Group",
    values_to = "Percentage"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("severe_pct", "moderate_pct", "mild_pct"),
      labels = c("Severe", "Moderate", "Mild")
    )
  )

scale_plot <- make_percentage_plot(
  plot_data = scale_plot_data,
  variable_order = scale_cols,
  fill_colors = c(
    "Mild" = "#BBDEFB",
    "Moderate" = "#64B5F6",
    "Severe" = "#1E88E5"
  ),
  title = "Scales",
  legend_title = "Severity Group",
  border_color = "white"
)

save_percentage_plot(
  plot = scale_plot,
  filename = "scale_severity_distribution.pdf",
  n_variables = length(scale_cols)
)


milestone_data <- read_phenotype_sheet("Developmental milestones")

missing_milestones <- setdiff(
  names(milestone_thresholds),
  colnames(milestone_data)
)

if (length(missing_milestones) > 0) {
  stop(
    "Missing developmental milestone columns: ",
    paste(missing_milestones, collapse = ", ")
  )
}

milestone_summary <- bind_rows(
  lapply(names(milestone_thresholds), function(col_name) {

    col_data <- na.omit(milestone_data[[col_name]])
    total_count <- length(col_data)
    threshold <- milestone_thresholds[col_name]

    data.frame(
      Variable = col_name,
      below_pct = if (total_count > 0) {
        mean(col_data <= threshold) * 100
      } else {
        NA_real_
      },
      above_pct = if (total_count > 0) {
        mean(col_data > threshold) * 100
      } else {
        NA_real_
      },
      total_count = total_count
    )
  })
)

milestone_plot_data <- milestone_summary %>%
  pivot_longer(
    cols = c(below_pct, above_pct),
    names_to = "Group",
    values_to = "Percentage"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("above_pct", "below_pct"),
      labels = c("> Threshold", "<= Threshold")
    )
  )

milestone_plot <- make_percentage_plot(
  plot_data = milestone_plot_data,
  variable_order = names(milestone_thresholds),
  fill_colors = c(
    "<= Threshold" = "#FCE4D6",
    "> Threshold" = "#F8C9A8"
  ),
  title = "Developmental milestones",
  legend_title = "Threshold Status",
  border_color = "#F8C9A8"
)

save_percentage_plot(
  plot = milestone_plot,
  filename = "developmental_milestone_threshold_distribution.pdf",
  n_variables = length(milestone_thresholds)
)


make_binary_percentage_plot <- function(
    sheet_name,
    title,
    fill_colors,
    filename,
    right_margin_cm = 1) {

  data <- read_phenotype_sheet(sheet_name)
  variable_order <- colnames(data)

  binary_summary <- bind_rows(
    lapply(variable_order, function(col_name) {

      col_data <- na.omit(data[[col_name]])
      total_count <- length(col_data)

      if (any(!col_data %in% c(0, 1))) {
        stop("Non-binary values found in: ", sheet_name, " / ", col_name)
      }

      data.frame(
        Variable = col_name,
        percent_0 = if (total_count > 0) {
          mean(col_data == 0) * 100
        } else {
          NA_real_
        },
        percent_1 = if (total_count > 0) {
          mean(col_data == 1) * 100
        } else {
          NA_real_
        },
        Availability = if (nrow(data) > 0) {
          total_count / nrow(data) * 100
        } else {
          NA_real_
        }
      )
    })
  )

  plot_data <- binary_summary %>%
    pivot_longer(
      cols = c(percent_0, percent_1),
      names_to = "Group",
      values_to = "Percentage"
    ) %>%
    mutate(
      Group = factor(
        Group,
        levels = c("percent_1", "percent_0"),
        labels = c("Present (1)", "Absent (0)")
      )
    )

  p <- make_percentage_plot(
    plot_data = plot_data,
    variable_order = variable_order,
    fill_colors = fill_colors,
    title = title,
    legend_title = "Diagnosis Status",
    border_color = unname(fill_colors["Present (1)"]),
    use_axis_break = TRUE,
    right_margin_cm = right_margin_cm
  )

  save_percentage_plot(
    plot = p,
    filename = filename,
    n_variables = length(variable_order)
  )

  p
}

mental_health_plot <- make_binary_percentage_plot(
  sheet_name = "Mental health and psychiatric",
  title = "Mental health and psychiatric conditions",
  fill_colors = c(
    "Absent (0)" = "#B2DFDB",
    "Present (1)" = "#80CBC4"
  ),
  filename = "mental_health_psychiatric_condition_distribution.pdf"
)

neurodevelopmental_plot <- make_binary_percentage_plot(
  sheet_name = "Neurodevelopmental and related",
  title = "Neurodevelopmental and related conditions",
  fill_colors = c(
    "Absent (0)" = "#FFF9C4",
    "Present (1)" = "#FFF176"
  ),
  filename = "neurodevelopmental_condition_distribution.pdf",
  right_margin_cm = 1.5
)

perinatal_plot <- make_binary_percentage_plot(
  sheet_name = "Perinatal complications",
  title = "Perinatal complications",
  fill_colors = c(
    "Absent (0)" = "#d9e1f2",
    "Present (1)" = "#b3c4e8"
  ),
  filename = "perinatal_complication_distribution.pdf"
)

birth_defect_plot <- make_binary_percentage_plot(
  sheet_name = "Birth defects",
  title = "Birth defects",
  fill_colors = c(
    "Absent (0)" = "#dcedc8",
    "Present (1)" = "#aed581"
  ),
  filename = "birth_defect_distribution.pdf"
)

growth_plot <- make_binary_percentage_plot(
  sheet_name = "Growth abnormalities",
  title = "Growth abnormalities",
  fill_colors = c(
    "Absent (0)" = "#f8bbd0",
    "Present (1)" = "#f48fb1"
  ),
  filename = "growth_abnormality_distribution.pdf"
)