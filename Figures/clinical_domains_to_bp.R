library(dplyr)
library(tidyr)
library(ggplot2)


work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)

###################### Stacked Bar Chart: Clinical Domain ######################

# Load data
data.phe_pathway <- read.csv(
  '../data/phenotype.pathway.weight.tsv',
  check.names = FALSE, sep = '\t'
)

# Clinical domains
p_cates <- c(
  'Scales', 'Developmental Milestones', 'Mental Health and Psychiatric Conditions',
  'Neurodevelopmental and Related Conditions', 'Perinatal Complications', 'Birth Defects',
  'Growth Abnormalities'
)

stacked <- data.phe_pathway %>%
  filter(!is.na(gene_count), gene_count > 0) %>%
  group_by(phenotype) %>%
  arrange(desc(gene_count), .by_group = TRUE) %>%
  mutate(
    xmin = cumsum(lag(gene_count, default = 0)),
    xmax = xmin + gene_count
  ) %>%
  ungroup() %>%
  mutate(phenotype = factor(phenotype, levels = p_cates))

bar_height <- 0.7

p <- ggplot(stacked) +
  geom_rect(
    aes(
      xmin = xmin, xmax = xmax,
      ymin = as.numeric(phenotype) - bar_height/2,
      ymax = as.numeric(phenotype) + bar_height/2,
      fill = pathway
    ),
    color = NA
  ) +
  scale_fill_manual(
    values = c(
      cognition = "#FBB4AE",
      `Rett syndrome` = "#B3CDE3",
      `chromatin remodeling` = "#CCEBC5",
      `brain development` = "#DECBE4"
    ),
    labels = c(
      cognition = "cognition",
      `Rett syndrome` = "Rett syndrome",
      `chromatin remodeling` = "chromatin remodeling",
      `brain development` = "brain development"
    ),
    name = NULL
  ) +
  scale_x_continuous(expand = c(0, 0.3)) +
  scale_y_continuous(
    breaks = seq_along(levels(stacked$phenotype)),
    labels = levels(stacked$phenotype)
  ) +
  labs(x = "the number of FDR-significant genes", y = NULL) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(hjust = 1),
    legend.position = "top",
    legend.key.size = unit(4, "mm")
  )

## save
ggsave('output/ExtendFig6b.pdf', p, width = 4, height = 2, dpi = 300, units = "in")
################################################################################
