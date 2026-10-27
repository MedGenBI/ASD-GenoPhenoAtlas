# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

work_dir <- file.path(base_dir, 'plots')
setwd(work_dir)


association_results_path <- 'data/geno_pheno_association_results.tsv'

trait_matrix_path <- "data/phenotype.csv"
trait_metadata_path <- "data/trait_metadata.csv"

required_trait_metadata_columns <- c(
  "trait_id", "trait_label", "trait_type", "clinical_domain"
)

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
  select(-any_of(c("Trait category", "Trait domain"))) %>%
  mutate(efdr_significant = as.logical(as.character(efdr_significant))) %>%
  left_join(
    trait_metadata %>%
      select(trait_id, trait_label, trait_type, clinical_domain),
    by = "trait_id"
  )

unmapped_association_trait_ids <- association_results %>%
  filter(is.na(trait_label)) %>%
  distinct(trait_id) %>%
  pull(trait_id)
if (length(unmapped_association_trait_ids) > 0) {
  stop(
    "Trait IDs missing from trait_metadata: ",
    paste(unmapped_association_trait_ids, collapse = ", ")
  )
}

association_results <- association_results %>%
  group_by(gene) %>%
  filter(any(efdr_significant, na.rm = TRUE)) %>%
  ungroup()

trait_matrix <- read.csv(trait_matrix_path, check.names = FALSE)


figure_width_in   <- 7.2
figure_height_in  <- 10.0

base_font_size_pt      <- 6.5
axis_x_text_size_pt    <- 6.0
axis_y_text_size_pt    <- 6.0
bar_axis_text_size_pt  <- 6.0

bubble_size_range   <- c(1, 3)
stacked_bar_width   <- 0.7

panel_gap_pt <- 1
panel_bottom_margin_pt <- 18

panel_margin_internal <- margin(
  panel_gap_pt, panel_gap_pt, panel_gap_pt, panel_gap_pt,
  unit = "pt"
)
panel_margin_bottom_row <- margin(
  panel_gap_pt, panel_gap_pt, panel_bottom_margin_pt, panel_gap_pt,
  unit = "pt"
)

# Output path
pdf_output_path <- "bubble.pdf"


clinical_domain_order <- c(
  "Scales", "Developmental Milestones", "Mental Health and Psychiatric Conditions",
  "Neurodevelopmental and Related Conditions", "Perinatal Complications",
  "Birth Defects", "Growth Abnormalities"
)

unexpected_clinical_domains <- setdiff(
  unique(trait_metadata$clinical_domain),
  clinical_domain_order
)
missing_clinical_domains <- setdiff(
  clinical_domain_order,
  unique(trait_metadata$clinical_domain)
)
if (length(unexpected_clinical_domains) > 0 || length(missing_clinical_domains) > 0) {
  stop(
    "clinical_domain values in trait_metadata do not match clinical_domain_order."
  )
}

trait_id_order <- colnames(trait_matrix)[-1]
unmapped_trait_matrix_ids <- setdiff(trait_id_order, trait_metadata$trait_id)
if (length(unmapped_trait_matrix_ids) > 0) {
  stop(
    "Trait IDs in trait_matrix missing from trait_metadata: ",
    paste(unmapped_trait_matrix_ids, collapse = ", ")
  )
}
trait_id_levels <- rev(trait_id_order)

association_results$trait_id <- factor(
  association_results$trait_id,
  levels = trait_id_levels
)
association_results$efdr_significant <- as.logical(as.character(association_results$efdr_significant))

clinical_domain_heights <- trait_metadata %>%
  filter(trait_id %in% trait_id_order) %>%
  distinct(clinical_domain, trait_id) %>%
  count(clinical_domain, name = "trait_count") %>%
  filter(clinical_domain %in% clinical_domain_order) %>%
  arrange(match(clinical_domain, clinical_domain_order)) %>%
  pull(trait_count)
if (sum(clinical_domain_heights) != length(trait_id_order)) {
  stop("Clinical-domain trait counts do not match the trait_matrix column count.")
}


gene_significance_summary <- association_results %>%
  group_by(gene) %>%
  summarise(
    nominal_p_sig_count = sum(p_value < 0.05, na.rm = TRUE),
    efdr_sig_count = sum(efdr_significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(efdr_sig_count), desc(nominal_p_sig_count))

genes_efdr_gt6  <- gene_significance_summary %>% filter(efdr_sig_count > 6)               %>% pull(gene)
genes_efdr_4to6 <- gene_significance_summary %>% filter(efdr_sig_count >= 4, efdr_sig_count <= 6) %>% pull(gene)
genes_efdr_2to3 <- gene_significance_summary %>% filter(efdr_sig_count >= 2, efdr_sig_count <= 3) %>% pull(gene)
genes_efdr_eq1  <- gene_significance_summary %>% filter(efdr_sig_count == 1)               %>% pull(gene)

gene_groups_by_efdr_count <- list(
  genes_efdr_gt6,
  genes_efdr_4to6,
  genes_efdr_2to3,
  genes_efdr_eq1
)
gene_group_widths <- lengths(gene_groups_by_efdr_count)

ordered_gene_levels <- unlist(gene_groups_by_efdr_count, use.names = FALSE)
association_results$gene <- factor(association_results$gene, levels = ordered_gene_levels)

significance_type_levels <- c("efdr_sig_count", "nominal_p_sig_count")

gene_significance_long <- gene_significance_summary %>%
  select(gene, nominal_p_sig_count, efdr_sig_count) %>%
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


trait_significance_long <- association_results %>%
  group_by(trait_id) %>%
  summarise(
    nominal_p_sig_count = sum(p_value < 0.05, na.rm = TRUE),
    efdr_sig_count = sum(efdr_significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
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


n_genes_total  <- length(ordered_gene_levels)
n_traits_total <- sum(clinical_domain_heights)
effect_size_limits <- range(abs(association_results$beta[association_results$p_value < 0.05]), na.rm = TRUE)

# Panel themes
bubble_theme_left_bottom <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x  = element_text(size = axis_x_text_size_pt),
  axis.text.y  = element_text(size = axis_y_text_size_pt),
  legend.position = "none",
  plot.margin = panel_margin_bottom_row
)
bubble_theme_bottom_row <- theme(
  axis.text.y = element_blank(), axis.ticks.y = element_blank(),
  axis.title.y = element_blank(), axis.title.x = element_blank(),
  axis.text.x  = element_text(size = axis_x_text_size_pt),
  legend.position = "none", plot.margin = panel_margin_bottom_row
)
bubble_theme_left_column <- theme(
  axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  axis.title.x = element_blank(), axis.title.y = element_blank(),
  axis.text.y  = element_text(size = axis_y_text_size_pt),
  legend.position = "none", plot.margin = panel_margin_internal
)
bubble_theme_internal <- theme(
  axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  axis.title.x = element_blank(), axis.text.y = element_blank(),
  axis.ticks.y = element_blank(), axis.title.y = element_blank(),
  legend.position = "none", plot.margin = panel_margin_internal
)

bubble_panel_list <- list()
for (col_idx in seq_along(gene_groups_by_efdr_count)) {
  genes_in_group <- gene_groups_by_efdr_count[[col_idx]]
  gene_group_data <- association_results %>% filter(gene %in% genes_in_group)
  
  for (row_idx in seq_along(clinical_domain_order)) {
    current_clinical_domain <- clinical_domain_order[row_idx]
    clinical_domain_panel_data <- gene_group_data %>%
      filter(clinical_domain == current_clinical_domain)
    
    bubble_panel <- ggplot(
      clinical_domain_panel_data,
      aes(x = gene, y = trait_id,
          size = ifelse(p_value < 0.05, abs(beta), 0))
    ) +
      geom_point(
        aes(
          colour = ifelse(beta > 0, "Positive", "Negative"),
          alpha  = factor(ifelse(efdr_significant, "EFDR", "Nominal"),
                          levels = c("Nominal","EFDR"))
        ),
        shape = 16
      ) +
      scale_alpha_manual(values = c("Nominal"=0.3,"EFDR"=1)) +
      scale_size_continuous(
        name   = "|Effect size|",
        breaks = c(0.1,0.5,1,1.5,2,3),
        limits = effect_size_limits,
        range = bubble_size_range,
      ) +
      scale_color_manual(
        values = c("Positive"="firebrick3","Negative"="#325fa5")
      ) +
      scale_y_discrete(labels = trait_label_by_id) +
      scale_x_discrete(limits = genes_in_group) +
      theme_bw(base_size = base_font_size_pt) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    is_left_column <- (col_idx == 1)
    is_bottom_row <- (row_idx == length(clinical_domain_order))
    bubble_panel <- bubble_panel + (
      if (is_left_column && is_bottom_row) bubble_theme_left_bottom
      else if (is_left_column)            bubble_theme_left_column
      else if (is_bottom_row)              bubble_theme_bottom_row
      else                                 bubble_theme_internal
    )
    
    bubble_panel_list[[length(bubble_panel_list) + 1]] <- bubble_panel
  }
}

bubble_matrix_plot <- wrap_plots(
  bubble_panel_list,
  ncol = length(gene_groups_by_efdr_count),
  nrow = length(clinical_domain_order),
  byrow = FALSE,
  heights = clinical_domain_heights,
  widths  = gene_group_widths
)


gene_significance_bar_limit <- max(gene_significance_long %>%
                        group_by(gene) %>%
                        summarise(total = sum(significant_count), .groups="drop") %>%
                        pull(total),
                      na.rm = TRUE)

gene_bar_panel_list <- list()
for (col_idx in seq_along(gene_groups_by_efdr_count)) {
  genes_in_group <- gene_groups_by_efdr_count[[col_idx]]
  
  gene_bar_panel <- ggplot(
    gene_significance_long,
    aes(x = gene, y = significant_count, fill = significance_type)
  ) +
    geom_col(position = position_stack(reverse = TRUE), width = stacked_bar_width) +
    scale_fill_manual(
      values = alpha(c(nominal_p_sig_count="#325fa5", efdr_sig_count="firebrick3"), 0.5),
      breaks = significance_type_levels,
      labels = c(nominal_p_sig_count = "P < 0.05", efdr_sig_count = "EFDR < 0.05"),
      name   = NULL
    ) +
    scale_x_discrete(limits = genes_in_group) +
    scale_y_continuous(limits = c(0, gene_significance_bar_limit)) +
    labs(x = NULL, y = "") +
    theme_bw(base_size = base_font_size_pt)
  
  gene_bar_theme_left <- theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line.y.left = element_line(),
    axis.title.y = element_text(size = base_font_size_pt),
    legend.position = "none",
    plot.margin = panel_margin_internal
  )
  gene_bar_theme_internal <- theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = panel_margin_internal
  )
  
  gene_bar_panel_list[[col_idx]] <- gene_bar_panel +
    if (col_idx == 1) gene_bar_theme_left else gene_bar_theme_internal
}

gene_bar_plot <- wrap_plots(gene_bar_panel_list, nrow = 1, widths = gene_group_widths)


trait_bar_limit <- max(trait_significance_long %>%
                          group_by(trait_id) %>%
                          summarise(total = sum(significant_count), .groups="drop") %>%
                         pull(total),
                       na.rm = TRUE)

trait_bar_panel_list <- list()
for (row_idx in seq_along(clinical_domain_order)) {
  current_clinical_domain <- clinical_domain_order[row_idx]
  trait_ids_in_clinical_domain <- trait_metadata %>%
    filter(clinical_domain == current_clinical_domain) %>%
    mutate(trait_id = factor(trait_id, levels = trait_id_levels)) %>%
    arrange(trait_id) %>%
    pull(trait_id) %>%
    as.character()
  
  trait_bar_panel <- ggplot(
    trait_significance_long,
    aes(x = significant_count, y = trait_id, fill = significance_type)
  ) +
    geom_col(position = position_stack(reverse = TRUE), width = stacked_bar_width) +
    scale_fill_manual(
      values = alpha(c(nominal_p_sig_count="#325fa5", efdr_sig_count="firebrick3"), 0.5),
      breaks = significance_type_levels,
      labels = c(nominal_p_sig_count = "P < 0.05", efdr_sig_count = "EFDR < 0.05"),
      name   = NULL
    ) +
    scale_y_discrete(
      limits = trait_ids_in_clinical_domain,
      labels = trait_label_by_id
    ) +
    scale_x_continuous(limits = c(0, trait_bar_limit)) +
    labs(x = if (row_idx == length(clinical_domain_order)) "" else NULL, y = NULL) +
    theme_bw(base_size = base_font_size_pt)
  
  trait_bar_theme_internal <- theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid   = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = panel_margin_internal
  )
  trait_bar_theme_bottom <- theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = bar_axis_text_size_pt),
    panel.grid   = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(),
    legend.position = "none",
    plot.margin = panel_margin_bottom_row
  )
  
  trait_bar_panel_list[[row_idx]] <- trait_bar_panel +
    if (row_idx == length(clinical_domain_order)) trait_bar_theme_bottom else trait_bar_theme_internal
}

trait_bar_plot <- wrap_plots(
  trait_bar_panel_list,
  ncol = 1,
  heights = clinical_domain_heights
)


final_plot <- wrap_plots(
  A = gene_bar_plot,
  B = trait_bar_plot,
  C = bubble_matrix_plot,
  design = "
  A#
  CB
  ",
  heights = c(5, n_traits_total),
  widths  = c(n_genes_total, 5)
)

# Save as PDF
ggsave(
  filename = pdf_output_path,
  plot     = final_plot,
  width    = figure_width_in,
  height   = figure_height_in,
  units    = "in",
  device   = cairo_pdf
)
