# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)


# Load data
assoc_path <- "../data/geno_pheno_associations.fdr_48genes.tsv"
phenotypes_path <- "../data/phenotype/spark_clinical_infor_44962.phenotype.csv"

assoc      <- read.csv(assoc_path, check.names = FALSE, sep = "\t")
phenotypes <- read.csv(phenotypes_path, check.names = FALSE)

fig_width_in   <- 7.2
fig_height_in  <- 10.0

base_pt        <- 7.5
axis_x_pt      <- 7.0
axis_y_pt      <- 7.0
bar_text_pt    <- 7.0
title_pt       <- 8.0
family_font    <- NULL

bubble_range   <- c(1, 4)
bar_width      <- 0.7

# Output path
pdf_out <- "output/Fig3.pdf"


# Clinical domains
cate_order <- c(
  "Scales", "Developmental Milestones", "Mental Health and Psychiatric Conditions",
  "Neurodevelopmental and Related Conditions", "Perinatal Complications",
  "Birth Defects", "Growth Abnormalities"
)

phenotype_order <- colnames(phenotypes)[2:79]
phenotype_levels <- rev(phenotype_order)

assoc$phenotype <- factor(assoc$phenotype, levels = phenotype_levels)

cate_heights <- assoc %>%
  distinct(phe_cate, phenotype) %>%
  count(phe_cate, name = "n_pheno") %>%
  filter(phe_cate %in% cate_order) %>%
  arrange(match(phe_cate, cate_order)) %>%
  pull(n_pheno)

gene_sig_summary <- assoc %>%
  group_by(gene) %>%
  summarise(
    sig_p = sum(p < 0.05, na.rm = TRUE),
    sig_q = sum(p_fdr < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sig_q), desc(sig_p))

genes_qgt5  <- gene_sig_summary %>% filter(sig_q > 5)                %>% pull(gene)
genes_q3to5 <- gene_sig_summary %>% filter(sig_q >= 3, sig_q <= 5)   %>% pull(gene)
genes_q2    <- gene_sig_summary %>% filter(sig_q == 2)               %>% pull(gene)
genes_q1    <- gene_sig_summary %>% filter(sig_q == 1)               %>% pull(gene)

gene_groups <- list(genes_qgt5, genes_q3to5, genes_q2, genes_q1)
gene_group_widths <- lengths(gene_groups)

gene_levels <- unlist(gene_groups, use.names = FALSE)
assoc$gene  <- factor(assoc$gene, levels = gene_levels)

gene_sig_long <- gene_sig_summary %>%
  transmute(gene, p_sig = sig_p, q_sig = sig_q) %>%
  pivot_longer(c(p_sig, q_sig), names_to = "sig_type", values_to = "count")

pheno_sig_long <- assoc %>%
  group_by(phenotype) %>%
  summarise(
    p_sig = sum(p < 0.05, na.rm = TRUE),
    q_sig = sum(p_fdr < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(p_sig, q_sig), names_to = "sig_type", values_to = "count")


# Bubble panels
size_limits <- range(abs(assoc$beta[assoc$p < 0.05]), na.rm = TRUE)

t_both <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x  = element_text(size = axis_x_pt),
  axis.text.y  = element_text(size = axis_y_pt),
  legend.position = "none",
  plot.margin = margin(0,0,20,0)
)
t_x_only <- theme(
  axis.text.y = element_blank(), axis.ticks.y = element_blank(),
  axis.title.y = element_blank(), axis.title.x = element_blank(),
  axis.text.x  = element_text(size = axis_x_pt),
  legend.position = "none", plot.margin = margin(0,0,20,0)
)
t_y_only <- theme(
  axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  axis.title.x = element_blank(), axis.title.y = element_blank(),
  axis.text.y  = element_text(size = axis_y_pt),
  legend.position = "none", plot.margin = margin(0,0,0,0)
)
t_none <- theme(
  axis.text.x = element_blank(), axis.ticks.x = element_blank(),
  axis.title.x = element_blank(), axis.text.y = element_blank(),
  axis.ticks.y = element_blank(), axis.title.y = element_blank(),
  legend.position = "none", plot.margin = margin(0,0,0,0)
)

bubble_plots <- list()
for (col_idx in seq_along(gene_groups)) {
  genes_this <- gene_groups[[col_idx]]
  df_col <- assoc %>% filter(gene %in% genes_this)
  
  for (row_idx in seq_along(cate_order)) {
    cate <- cate_order[row_idx]
    df_pan <- df_col %>% filter(phe_cate == cate)
    
    p <- ggplot(
      df_pan,
      aes(x = gene, y = phenotype,
          size = ifelse(p < 0.05, abs(beta), 0))
    ) +
      geom_point(
        aes(
          colour = ifelse(beta > 0, "Positive", "Negative"),
          alpha  = factor(ifelse(p_fdr < 0.05, "FDR", "Nominal"),
                          levels = c("Nominal","FDR"))
        ),
        shape = 16
      ) +
      scale_alpha_manual(values = c("Nominal"=0.3,"FDR"=1)) +
      scale_size_continuous(
        name   = "|Effect size|",
        breaks = c(0.1,0.5,1,1.5,2,3),
        limits = size_limits,
        range = bubble_range,
      ) +
      scale_color_manual(
        values = c("Positive"="firebrick3","Negative"="#325fa5")
      ) +
      scale_x_discrete(limits = genes_this) +
      theme_bw(base_size = base_pt) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    is_left <- (col_idx == 1)
    is_bottom <- (row_idx == length(cate_order))
    p <- p + (
      if (is_left && is_bottom) t_both
      else if (is_left)         t_y_only
      else if (is_bottom)       t_x_only
      else                      t_none
    )
    
    bubble_plots[[length(bubble_plots) + 1]] <- p
  }
}

main_bubble_plot <- wrap_plots(
  bubble_plots,
  ncol = length(gene_groups),
  nrow = length(cate_order),
  byrow = FALSE,
  heights = cate_heights,
  widths  = gene_group_widths
)


# TOP stacked gene bars
gene_bar_limit <- max(gene_sig_long %>%
                        group_by(gene) %>%
                        summarise(total = sum(count), .groups="drop") %>%
                        pull(total),
                      na.rm = TRUE)

top_bar_list <- list()
for (col_idx in seq_along(gene_groups)) {
  genes_this <- gene_groups[[col_idx]]
  
  bar <- ggplot(gene_sig_long, aes(x = gene, y = count, fill = sig_type)) +
    geom_col(position = "stack", width = bar_width) +
    scale_fill_manual(
      values = alpha(c(p_sig="#325fa5", q_sig="firebrick3"), 0.5),
      labels = c("p < 0.05", "q < 0.05"),
      name   = NULL
    ) +
    scale_x_discrete(limits = genes_this) +
    scale_y_continuous(limits = c(0, gene_bar_limit)) +
    labs(x = NULL, y = "") +
    theme_bw(base_size = base_pt)
  
  theme_left <- theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line.y.left = element_line(),
    axis.title.y = element_text(size = base_pt),
    legend.position = "none",
    plot.margin = margin(0,0,0,0)
  )
  theme_inner <- theme(
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
    plot.margin = margin(0,0,0,0)
  )
  
  top_bar_list[[col_idx]] <- bar + if (col_idx == 1) theme_left else theme_inner
}

top_bar_plot <- wrap_plots(top_bar_list, nrow = 1, widths = gene_group_widths)


# RIGHT stacked phenotype bars
pheno_bar_limit <- max(pheno_sig_long %>%
                         group_by(phenotype) %>%
                         summarise(total = sum(count), .groups="drop") %>%
                         pull(total),
                       na.rm = TRUE)

right_bar_list <- list()
for (row_idx in seq_along(cate_order)) {
  cate <- cate_order[row_idx]
  phen_this <- assoc %>%
    filter(phe_cate == cate) %>%
    distinct(phenotype) %>%
    arrange(phenotype) %>%
    pull(phenotype)
  
  bar <- ggplot(pheno_sig_long, aes(x = count, y = phenotype, fill = sig_type)) +
    geom_col(position = "stack", width = bar_width) +
    scale_fill_manual(
      values = alpha(c(p_sig="#325fa5", q_sig="firebrick3"), 0.5),
      labels = c("p < 0.05", "q < 0.05"),
      name   = NULL
    ) +
    scale_y_discrete(limits = phen_this) +
    scale_x_continuous(limits = c(0, pheno_bar_limit)) +
    labs(x = if (row_idx == length(cate_order)) "" else NULL, y = NULL) +
    theme_bw(base_size = base_pt)
  
  theme_mid <- theme(
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
    plot.margin = margin(0,0,0,0)
  )
  theme_bottom <- theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = bar_text_pt),
    panel.grid   = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(),
    legend.position = "none",
    plot.margin = margin(0,0,0,0)
  )
  
  right_bar_list[[row_idx]] <- bar + if (row_idx == length(cate_order)) theme_bottom else theme_mid
}

right_bar_plot <- wrap_plots(right_bar_list, ncol = 1, heights = cate_heights)


# Final layout
final_plot <- wrap_plots(
  A = top_bar_plot,
  B = right_bar_plot,
  C = main_bubble_plot,
  design = "
  A#
  CB
  ",
  heights = c(5, 78),
  widths  = c(48, 5)
)

# Save
ggsave(
  filename = pdf_out,
  plot     = final_plot,
  width    = fig_width_in,
  height   = fig_height_in,
  units    = "in",
  device   = cairo_pdf
)
