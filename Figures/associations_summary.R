library(dplyr)
library(tidyr)
library(ggplot2)


work_dir <- "/home/mgbi/projects/gpasd/plots"
setwd(work_dir)

association_results_path <- 'data/geno_pheno_association_results.tsv'
trait_metadata_path <- "data/trait_metadata.csv"

required_trait_metadata_columns <- c("trait_id", "trait_label")

trait_metadata <- read.csv(
  trait_metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

trait_label_by_id <- setNames(trait_metadata$trait_label, trait_metadata$trait_id)

association_results <- read.csv(
  association_results_path,
  check.names = FALSE,
  sep = "\t"
) %>%
  rename(
    trait_id = Trait,
    gene = Gene,
    beta = Beta,
    p_value = P,
    efdr_significant = Empirical_FDR_Sig
  ) %>%
  mutate(efdr_significant = as.logical(as.character(efdr_significant))) %>%
  left_join(
    trait_metadata %>% select(trait_id, trait_label),
    by = "trait_id"
  )

unmapped_association_trait_ids <- association_results %>%
  filter(is.na(trait_label)) %>%
  distinct(trait_id) %>%
  pull(trait_id)
if (length(unmapped_association_trait_ids) > 0) {
  stop(
    "The following Trait IDs are missing from trait_metadata: ",
    paste(unmapped_association_trait_ids, collapse = ", ")
  )
}

significance_type_levels <- c("efdr_sig_count", "nominal_p_sig_count")

############################# 1. Bar Chart: Genes ##############################

gene_association_status_summary <- data.frame(
  association_status = factor(
    c("eFDR < 0.05", "P < 0.05", "Not Significant"),
    levels = c("eFDR < 0.05", "P < 0.05", "Not Significant")
  ),
  gene_count = c(51, 179, 6)
)

association_status_colors <- c("eFDR < 0.05" = "#1E88E5",
                               "P < 0.05" = "#F8C9A8",
                               "Not Significant" = "grey80")

association_status_border_colors <- c(
  "eFDR < 0.05" = "#1565C0",
  "P < 0.05" = "#C97845",
  "Not Significant" = "#969696"
)
gene_status_bar_width <- 0.62

gene_association_status_plot <- ggplot(
  gene_association_status_summary,
  aes(x = association_status, y = gene_count, fill = association_status)
) +
  geom_col(
    aes(colour = association_status),
    width = gene_status_bar_width,
    linewidth = 0.35
  ) +
  geom_text(aes(label = gene_count), vjust = -0.7, size = 2, fontface = "plain", color = "black") +
  scale_fill_manual(values = association_status_colors) +
  scale_colour_manual(values = association_status_border_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13)), breaks = scales::pretty_breaks(n = 6)) +
  ggprism::theme_prism(base_size = 5.5, border = TRUE) +
  labs(x = NULL, y = "gene counts") +
  theme(
    
    axis.text = element_text(face = "plain", size = 5, color = "black"),
    axis.title.y = element_text(face = "plain", size = 5.5, margin = margin(r = 1)),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 5, l = 5),
    plot.title = element_text(face = "plain", size = 7, hjust = 0.5),

    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.15, colour = "black"),
    axis.ticks.length = unit(0.05, "cm")
  )

ggsave('association_summary_a.pdf', gene_association_status_plot, width = 2, height = 2, dpi = 300, units = "in")


######################### 2. Stacked Bar Chart: Genes ##########################
gene_association_results <- association_results %>%
  group_by(gene) %>%
  filter(any(efdr_significant, na.rm = TRUE)) %>%
  ungroup()

gene_significance_summary <- gene_association_results %>%
  group_by(gene) %>%
  summarise(
    nominal_p_sig_count = sum(p_value < 0.05, na.rm = TRUE),
    efdr_sig_count = sum(efdr_significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%

  arrange(desc(efdr_sig_count), desc(nominal_p_sig_count)) %>%
  mutate(gene = factor(gene, levels = gene))

gene_significance_long <- gene_significance_summary %>%
  pivot_longer(
    c(nominal_p_sig_count, efdr_sig_count),
    names_to = "significance_type",
    values_to = "significant_count"
  ) %>%
  mutate(
    significance_type = factor(
      significance_type,
      levels = significance_type_levels
    )
  )

gene_significance_bar_plot <- ggplot(
  gene_significance_long,
  aes(x = gene, y = significant_count, fill = significance_type)
) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  scale_fill_manual(
    breaks = significance_type_levels,
    values = c(nominal_p_sig_count = "#BBDEFB", efdr_sig_count = "#1E88E5"),
    labels = c(nominal_p_sig_count = "P < 0.05", efdr_sig_count = "eFDR < 0.05"),
    name = NULL
  ) +
  labs(x = NULL, y = "Significant counts") +
  scale_y_continuous(expand = c(0, 1)) +
  theme_bw(base_size = 5.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5, margin = margin(r=1)),
    legend.text  = element_text(size = 6),
    legend.title = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.position = "inside",
    legend.position.inside = c(.85, .85),
    panel.grid = element_blank()
  )

ggsave('association_summary_b.pdf', gene_significance_bar_plot, width = 4.69, height = 2, dpi = 300, units = "in")


######################## 3. Stacked Bar Chart: Traits ##########################
trait_significance_summary <- association_results %>%
  group_by(trait_id, trait_label) %>%
  summarise(
    nominal_p_sig_count = sum(p_value < 0.05, na.rm = TRUE),
    efdr_sig_count = sum(efdr_significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(efdr_sig_count), desc(nominal_p_sig_count)) %>%
  mutate(trait_id = factor(trait_id, levels = trait_id))

trait_significance_long <- trait_significance_summary %>%
  pivot_longer(
    c(nominal_p_sig_count, efdr_sig_count),
    names_to = "significance_type",
    values_to = "significant_count"
  ) %>%
  mutate(
    significance_type = factor(
      significance_type,
      levels = significance_type_levels
    )
  )

trait_significance_bar_plot <- ggplot(
  trait_significance_long,
  aes(x = trait_id, y = significant_count, fill = significance_type)
) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  scale_fill_manual(
    breaks = significance_type_levels,
    values = c(nominal_p_sig_count = "#BBDEFB", efdr_sig_count = "#1E88E5"),
    labels = c(nominal_p_sig_count = "P < 0.05", efdr_sig_count = "eFDR < 0.05"),
    name = NULL
  ) +
  scale_x_discrete(labels = trait_label_by_id) +
  labs(x = NULL, y = "Significant counts") +
  scale_y_continuous(expand = c(0, 1)) +
  theme_bw(base_size = 5.5) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    legend.text  = element_text(size = 5),
    legend.title = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.position = "inside",
    legend.position.inside = c(.85, .85),
    panel.grid = element_blank()
  )

ggsave('association_summary_c.pdf', trait_significance_bar_plot, width = 6.69, height = 2.8, dpi = 300, units = "in")
