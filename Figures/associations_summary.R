library(dplyr)
library(tidyr)
library(ggplot2)


work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)


############################# 1. Bar Chart: Genes ##############################

bar_df <- data.frame(
  category = factor(
    c("FDR < 0.05", "p < 0.05", "Not Significant"),
    levels = c("FDR < 0.05", "p < 0.05", "Not Significant")
  ),
  count = c(48, 173, 13)
)

bar_colors <- c("FDR < 0.05" = "#1E88E5",
                "p < 0.05" = "#F8C9A8",
                "Not Significant" = "grey80")

p <- ggplot(bar_df, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.62, color = "black", linewidth = 0) + 
  geom_text(aes(label = count), vjust = -0.7, size = 2, fontface = "plain", color = "black") + 
  scale_fill_manual(values = bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.13)), breaks = scales::pretty_breaks(n = 6)) +
  ggprism::theme_prism(base_size = 5.5, border = TRUE) +   
  labs(x = NULL, y = "gene counts") +
  theme(
    
    axis.text = element_text(face = "plain", size = 5, color = "black"),
    axis.title.y = element_text(face = "plain", size = 5.5, margin = margin(r = 1)),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 5, b = 5, l = 5),
    plot.title = element_text(face = "plain", size = 7, hjust = 0.5),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.15, color = "black"),
    axis.ticks = element_line(linewidth = 0.15),
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

## Save
ggsave('output/ExtendFig4a.pdf', p, width = 2, height = 2, dpi = 300, units = "in")
################################################################################


######################### 2. Stacked Bar Chart: Genes ##########################
data_fdr <- read.csv('../data/geno_pheno_associations.fdr_48genes.tsv', check.names = F, sep = '\t')

sig_counts <- data_fdr %>%
  group_by(gene) %>%
  summarise(
    p_sig = sum(p < 0.05, na.rm = TRUE),
    q_sig = sum(p_fdr < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(q_sig), desc(p_sig)) %>%
  mutate(gene = factor(gene, levels = gene)) %>%
  pivot_longer(c(p_sig, q_sig), names_to = "sig_type", values_to = "count")

p <- ggplot(sig_counts, aes(x = gene, y = count, fill = sig_type)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(
    values = c(p_sig = "#BBDEFB", q_sig = "#1E88E5"),
    labels = c("p_sig" = "p < 0.05", "q_sig" = "FDR < 0.05"),
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

## Save
ggsave('output/ExtendFig4b.pdf', p, width = 4.69, height = 2, dpi = 300, units = "in")
################################################################################



######################## 3. Stacked Bar Chart: Traits ##########################
data_p <- read.csv('../data/geno_pheno_associations.p_173genes.tsv' , check.names = F, sep = '\t')

sig_counts <- data_p %>%
  group_by(phenotype) %>%
  summarise(
    p_sig = sum(p < 0.05, na.rm = TRUE),
    q_sig = sum(p_fdr < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sort_key_q = q_sig, sort_key_p = p_sig) %>%
  arrange(desc(sort_key_q), desc(sort_key_p)) %>%
  mutate(phenotype = factor(phenotype, levels = phenotype)) %>%
  pivot_longer(c(p_sig, q_sig), names_to = "sig_type", values_to = "count")

p <- ggplot(sig_counts, aes(x = phenotype, y = count, fill = sig_type)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(
    values = c(p_sig = "#BBDEFB", q_sig = "#1E88E5"),
    labels = c("p_sig" = "p < 0.05", "q_sig" = "FDR < 0.05"),
    name = NULL
  ) +
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

## Save
ggsave('output/ExtendFig4c.pdf', p, width = 6.69, height = 2.8, dpi = 300, units = "in")
################################################################################
