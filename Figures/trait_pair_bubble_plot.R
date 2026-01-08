rm(list = ls())
library(ggplot2);library(dplyr);library(tidyr);library(patchwork);library(openxlsx)

assoc_data <- read.csv("../PLP.csv")

assoc_data <- assoc_data %>%
  rename(Exposed_proportion = "Exposed.P/LP.proportion.(%)") %>%
  mutate(
    Exposed_proportion = ifelse(is.na(Exposed_proportion), 0, Exposed_proportion),
    Exposed_proportion = ifelse(Exposed_proportion > 25, 25, Exposed_proportion),
    logFDR = -log10(FDR),
    logFDR = ifelse(logFDR > 10, 10, logFDR),
    significance = ifelse(FDR < 0.05, "significant", "non-significant")
  )

pheno_data <- read.csv("../category.csv")
sig_pheno <- unique(pheno_data$Trait[pheno_data$FDR < 0.05])
pheno_order <- sig_pheno
all_phenos <- unique(c(assoc_data$Trait1, assoc_data$Trait2))
pheno_order <- c(pheno_order, setdiff(all_phenos, pheno_order))

assoc_data <- assoc_data %>%
  rowwise() %>%
  mutate(pheno_pair = paste(sort(c(Trait1, Trait2)), collapse = "|")) %>%
  distinct(pheno_pair, .keep_all = TRUE) %>%
  ungroup() %>%
  select(-pheno_pair)

assoc_data$Trait1 <- factor(assoc_data$Trait1, levels = pheno_order)
assoc_data$Trait2 <- factor(assoc_data$Trait2, levels = rev(pheno_order))

assoc_data <- assoc_data %>% filter(!is.na(Trait1) & !is.na(Trait2))

ggplot(assoc_data, aes(x = Trait1, y = Trait2)) +
  geom_point(data = filter(assoc_data, significance == "non-significant"),
             aes(size = Exposed_proportion),
             shape = 21, color = "white", fill = "lightgrey", alpha = 0.7) +
  geom_point(data = filter(assoc_data, significance == "significant"),
             aes(size = Exposed_proportion, fill = logFDR),
             shape = 21, color = "white", alpha = 1) +
  scale_size_continuous(name = "Exposed P/LP proportion (%)",
                        range = c(0, 5),
                        breaks = seq(0, 25, by = 5),
                        limits = c(0, 25)) +
  scale_fill_gradientn(name = "-log10(FDR)",
                       colours = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15"),
                       values = scales::rescale(c(0, 2.5, 5, 7.5, 10)),
                       limits = c(0, 10),
                       na.value = "lightgrey") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90", size = 0.2),
        legend.position = "right",
        panel.border = element_rect(fill = NA, color = "grey80")) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(override.aes = list(color = "black", fill = "black")))

ggsave("bubble.pdf",plot = main_plot,width = 297,height = 210,units = "mm",device = "pdf")
