rm(list = ls())
library(tidyr); library(dplyr); library(ggplot2); library(openxlsx); library(cowplot)

sig <- read.csv("../FDR.csv");plp <- read.csv("../PLP.csv");raw <- read.csv("../44962.csv")
row.names(raw) <- raw$spid;raw <- cbind(raw, plp = ifelse(raw[, "spid"] %in% unique(plp[, "spid"]), "T", "F"));raw2 <- raw

# scales: a
scales <- raw2[, c(13:32, 91)];scale_df <- as.data.frame(scales, stringsAsFactors = FALSE)
vars <- colnames(scale_df)[1:20];new_vars <- gsub("_", " ", vars);colnames(scale_df)[1:20] <- new_vars;vars <- new_vars
plp_char <- as.character(scale_df[[21]]);plp_bool <- toupper(plp_char) %in% c("T", "TRUE", "1", "YES", "Y")

bin_stats_one_var_scales <- function(v, idx) {
  x <- scale_df[[v]]
  ok <- !is.na(x) & !is.na(plp_bool)
  x <- x[ok]; y <- plp_bool[ok]
  if (length(x) == 0) return(NULL)
  qs <- quantile(x, probs = c(0.25, 0.75), type = 7, na.rm = TRUE)
  if (qs[1] < qs[2]) {
    grp <- cut(x, breaks = c(-Inf, qs[1], qs[2], Inf), labels = c("0_25", "25_75", "75_100"), include.lowest = TRUE, right = TRUE)
  } else {
    pr <- rank(x, ties.method = "average") / length(x)
    grp <- cut(pr, breaks = c(-Inf, 0.25, 0.75, Inf), labels = c("0_25", "25_75", "75_100"), include.lowest = TRUE, right = TRUE)
  }
  if (idx <= 11) {
    severity_map <- c("0_25" = "mild", "25_75" = "moderate", "75_100" = "severe")
  } else {
    severity_map <- c("0_25" = "severe", "25_75" = "moderate", "75_100" = "mild")
  }
  bins <- c("0_25", "25_75", "75_100")
  out <- lapply(bins, function(b) {
    z <- y[grp == b]; n <- length(z); carrier <- if (n > 0) sum(z) else 0
    prop <- if (n > 0) carrier / n else NA_real_
    se <- if (n > 0) sqrt(prop * (1 - prop) / n) else NA_real_
    data.frame(var = v, bin = severity_map[b], n = n, carrier = carrier, prop = prop, se = se,
               lower = if (!is.na(se)) pmax(0, prop - 1.96 * se) else NA_real_,
               upper = if (!is.na(se)) pmin(1, prop + 1.96 * se) else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

plot_df_list <- lapply(seq_along(vars), function(i) bin_stats_one_var_scales(vars[i], i))
plot_df <- do.call(rbind, plot_df_list)
plot_df$var <- factor(plot_df$var, levels = vars)
plot_df$bin <- factor(plot_df$bin, levels = c("mild", "moderate", "severe"))
fill_colors <- c("mild" = "#BBDEFB", "moderate" = "#64B5F6", "severe" = "#1E88E5")
y_top <- 16

p <- ggplot(plot_df, aes(x = var, y = prop * 100, fill = bin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.78), width = 0.58, color = "black") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),
                position = position_dodge(width = 0.78), width = 0.45, size = 0.4, color = "black") +
  geom_hline(yintercept = 8, linetype = "dashed", color = "black", size = 0.4) +
  scale_fill_manual(values = fill_colors, name = "Severity") +
  scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, 2), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "", title = NULL) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 11, color = "black"),
        legend.position = "none", legend.title = element_text(size = 10), legend.text = element_text(size = 9))

sig_clean <- sig %>% mutate(phenotype_clean = gsub("_", " ", phenotype), FDR = as.numeric(FDR), significant = FDR < 0.05) %>% filter(significant)
top_df <- aggregate(upper ~ var, data = plot_df, FUN = max, na.rm = TRUE)
top_df$y_base <- top_df$upper * 100
ann <- merge(top_df, sig_clean, by.x = "var", by.y = "phenotype_clean", all.x = FALSE)

if (nrow(ann) > 0) {
  ann$xi <- as.numeric(factor(ann$var, levels = levels(plot_df$var)))
  ann <- ann %>% group_by(var, baseline, control) %>% mutate(comp_id = paste(baseline, control, sep = "_")) %>%
    ungroup() %>% group_by(var) %>% mutate(comp_index = match(comp_id, unique(comp_id)), height_offset = 0.8 * comp_index) %>% ungroup()
  ann$y_lo <- ann$y_base + 0.3 + ann$height_offset
  ann$y_hi <- ann$y_base + 0.6 + ann$height_offset
  ann$y_text <- ann$y_base + 0.7 + ann$height_offset
  seg_df <- data.frame(); text_df <- data.frame()
  for (i in 1:nrow(ann)) {
    xi <- ann$xi[i]; y_lo <- ann$y_lo[i]; y_hi <- ann$y_hi[i]; y_text <- ann$y_text[i]
    if (ann$control[i] == "moderate") { w <- 0.05 } else if (ann$control[i] == "severe") { w <- 0.30 } else { w <- 0.30 }
    seg_df <- rbind(seg_df,
                    data.frame(x = xi - 0.3, xend = xi - 0.3, y = y_lo, yend = y_hi, group = i),
                    data.frame(x = xi + w, xend = xi + w, y = y_lo, yend = y_hi, group = i),
                    data.frame(x = xi - 0.3, xend = xi + w, y = y_hi, yend = y_hi, group = i))
    to_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
    label <- to_stars(ann$FDR[i])
    text_df <- rbind(text_df, data.frame(x = xi, y = y_text, label = label, group = i))
  }
  new_top <- max(y_top, max(ann$y_text + 0.25, na.rm = TRUE))
  break_top <- ceiling(new_top / 2) * 2
  p1 <- p +
    geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = yend, group = group),
                 inherit.aes = FALSE, size = 0.4) +
    geom_text(data = text_df, aes(x = x, y = y, label = label, group = group),
              inherit.aes = FALSE, size = 4) +
    scale_y_continuous(limits = c(0, new_top), breaks = seq(0, break_top, by = 2), expand = expansion(mult = c(0, 0)))
} else {
  p1 <- p
}

# developmental milestones: b
scales <- raw2[, c(33:42, 91)];scale_df <- as.data.frame(scales, stringsAsFactors = FALSE)
vars <- colnames(scale_df)[1:10];new_vars <- gsub("_", " ", vars)
colnames(scale_df)[1:10] <- new_vars;vars <- new_vars
plp_char <- as.character(scale_df[[11]])
plp_bool <- toupper(plp_char) %in% c("T", "TRUE", "1", "YES", "Y")

thresholds <- c("Smiling" = 4, "Sitting upright" = 8, "Crawling" = 10, "Walking" = 15,
                "Spoon feeding self" = 18, "Speaking first word" = 18, "Speaking with combined word" = 27,
                "Speaking first phrase" = 36, "Attaining bladder control" = 42, "Attaining bowel control" = 48)

bin_stats_one_var_milestone <- function(v) {
  x <- scale_df[[v]]
  ok <- !is.na(x) & !is.na(plp_bool)
  x <- x[ok]; y <- plp_bool[ok]
  if (length(x) == 0) return(NULL)
  threshold_value <- thresholds[v]
  grp <- ifelse(x <= threshold_value, "0_95", "95_100")
  bins <- c("0_95", "95_100")
  out <- lapply(bins, function(b) {
    z <- y[grp == b]; n <- length(z); carrier <- if (n > 0) sum(z) else 0
    prop <- if (n > 0) carrier / n else NA_real_
    se <- if (n > 0) sqrt(prop * (1 - prop) / n) else NA_real_
    data.frame(var = v, bin = b, n = n, carrier = carrier, prop = prop, se = se,
               lower = if (!is.na(se)) pmax(0, prop - 1.96 * se) else NA_real_,
               upper = if (!is.na(se)) pmin(1, prop + 1.96 * se) else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

plot_df <- do.call(rbind, lapply(vars, bin_stats_one_var_milestone))
plot_df$var <- factor(plot_df$var, levels = vars)
plot_df$bin <- factor(plot_df$bin, levels = c("0_95", "95_100"))
fill_colors <- c("0_95" = "#FCE4D6", "95_100" = "#F8C9A8")
y_top <- 18

p <- ggplot(plot_df, aes(x = var, y = prop * 100, fill = bin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.78), width = 0.58, color = "black") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),
                position = position_dodge(width = 0.78), width = 0.4, size = 0.4, color = "black") +
  geom_hline(yintercept = 8, linetype = "dashed", color = "black", size = 0.4) +
  scale_fill_manual(values = fill_colors, name = NULL) +
  scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, 2), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "", title = NULL) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 11, color = "black"),
        legend.position = "none")

sig_df <- as.data.frame(sig, stringsAsFactors = FALSE)
sig_df$phenotype <- as.character(sig_df$phenotype)
sig_df$phenotype_clean <- gsub("_", " ", sig_df$phenotype)
sig_df$FDR <- as.numeric(sig_df$FDR)
top_df <- aggregate(upper ~ var, data = plot_df, FUN = max, na.rm = TRUE)
top_df$y_base <- top_df$upper * 100
ann <- merge(top_df, sig_df[, c("phenotype_clean", "FDR")], by.x = "var", by.y = "phenotype_clean", all.x = FALSE, all.y = FALSE)
ann <- ann[!is.na(ann$FDR) & ann$FDR < 0.05, , drop = FALSE]

if (nrow(ann) > 0) {
  ann$xi <- as.numeric(factor(ann$var, levels = levels(plot_df$var)))
  pad1 <- 0.30; pad2 <- 0.60; pad3 <- 0.90
  ann$y_lo <- ann$y_base + pad1
  ann$y_hi <- ann$y_base + pad2
  ann$y_text <- ann$y_base + pad3
  w <- 0.22
  seg_df <- rbind(
    data.frame(x = ann$xi - w, xend = ann$xi - w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi + w, xend = ann$xi + w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi - w, xend = ann$xi + w, y = ann$y_hi, yend = ann$y_hi)
  )
  to_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  ann$label <- to_stars(ann$FDR)
  new_top <- max(y_top, max(ann$y_text + 0.25, na.rm = TRUE))
  break_top <- ceiling(new_top / 2) * 2
  p2 <- p +
    geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, size = 0.4) +
    geom_text(data = ann, aes(x = xi, y = y_text, label = label),
              inherit.aes = FALSE, size = 4) +
    scale_y_continuous(limits = c(0, new_top), breaks = seq(0, break_top, by = 2), expand = expansion(mult = c(0, 0)))
} else {
  p2 <- p
}

# mental health and psychiatric: c
scales <- raw2[, c(43:61, 91)]
scale_df <- as.data.frame(scales, stringsAsFactors = FALSE)
vars <- colnames(scale_df)[1:19];new_vars <- gsub("_", " ", vars);colnames(scale_df)[1:19] <- new_vars;vars <- new_vars
plp_char <- as.character(scale_df[[20]]);plp_bool <- toupper(plp_char) %in% c("T", "TRUE", "1", "YES", "Y")

bin_stats_one_var_binary <- function(v) {
  x <- scale_df[[v]]
  ok <- !is.na(x) & !is.na(plp_bool)
  x <- x[ok]; y <- plp_bool[ok]
  if (length(x) == 0) return(NULL)
  grp <- factor(x, levels = c(0, 1), labels = c("Absent (0)", "Present (1)"))
  bins <- c("Absent (0)", "Present (1)")
  out <- lapply(bins, function(b) {
    z <- y[grp == b]; n <- length(z); carrier <- if (n > 0) sum(z) else 0
    prop <- if (n > 0) carrier / n else NA_real_
    se <- if (n > 0) sqrt(prop * (1 - prop) / n) else NA_real_
    data.frame(var = v, bin = b, n = n, carrier = carrier, prop = prop, se = se,
               lower = if (!is.na(se)) pmax(0, prop - 1.96 * se) else NA_real_,
               upper = if (!is.na(se)) pmin(1, prop + 1.96 * se) else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

plot_df <- do.call(rbind, lapply(vars, bin_stats_one_var_binary))
plot_df$var <- factor(plot_df$var, levels = vars)
plot_df$bin <- factor(plot_df$bin, levels = c("Absent (0)", "Present (1)"))
fill_colors <- c("Absent (0)" = "#B2DFDB", "Present (1)" = "#80CBC4")
y_top <- 14

p <- ggplot(plot_df, aes(x = var, y = prop * 100, fill = bin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.78), width = 0.58, color = "black") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),
                position = position_dodge(width = 0.78), width = 0.5, size = 0.4, color = "black") +
  geom_hline(yintercept = 8, linetype = "dashed", color = "black", size = 0.4) +
  scale_fill_manual(values = fill_colors, name = NULL) +
  scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, 2), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "", title = NULL) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 11, color = "black"),
        legend.position = "none")

sig_df <- as.data.frame(sig, stringsAsFactors = FALSE);sig_df$phenotype <- as.character(sig_df$phenotype)
sig_df$phenotype_clean <- gsub("_", " ", sig_df$phenotype);sig_df$FDR <- as.numeric(sig_df$FDR)
plot_vars <- levels(plot_df$var);sig_df <- sig_df[sig_df$phenotype_clean %in% plot_vars, ]
top_df <- aggregate(upper ~ var, data = plot_df, FUN = max, na.rm = TRUE);top_df$y_base <- top_df$upper * 100
ann <- merge(top_df, sig_df[, c("phenotype_clean", "FDR")], by.x = "var", by.y = "phenotype_clean", all.x = FALSE, all.y = FALSE)
ann <- ann[!is.na(ann$FDR) & ann$FDR < 0.05, , drop = FALSE]

if (nrow(ann) > 0) {
  ann$xi <- as.numeric(factor(ann$var, levels = levels(plot_df$var)))
  pad1 <- 0.30; pad2 <- 0.60; pad3 <- 0.90
  ann$y_lo <- ann$y_base + pad1
  ann$y_hi <- ann$y_base + pad2
  ann$y_text <- ann$y_base + pad3
  w <- 0.25
  seg_df <- rbind(
    data.frame(x = ann$xi - w, xend = ann$xi - w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi + w, xend = ann$xi + w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi - w, xend = ann$xi + w, y = ann$y_hi, yend = ann$y_hi)
  )
  to_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  ann$label <- to_stars(ann$FDR)
  new_top <- max(y_top, max(ann$y_text + 0.25, na.rm = TRUE))
  break_top <- ceiling(new_top / 2) * 2
  p3 <- p +
    geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, size = 0.4) +
    geom_text(data = ann, aes(x = xi, y = y_text, label = label),
              inherit.aes = FALSE, size = 4) +
    scale_y_continuous(limits = c(0, new_top), breaks = seq(0, break_top, by = 2), expand = expansion(mult = c(0, 0)))
} else {
  p3 <- p
}

# neurodevelopmental and related: d
scales <- raw2[, c(62:71, 91)];scale_df <- as.data.frame(scales, stringsAsFactors = FALSE)
vars <- colnames(scale_df)[1:10];new_vars <- gsub("_", " ", vars)
colnames(scale_df)[1:10] <- new_vars;vars <- new_vars
plp_char <- as.character(scale_df[[11]])
plp_bool <- toupper(plp_char) %in% c("T", "TRUE", "1", "YES", "Y")
plot_df <- do.call(rbind, lapply(vars, bin_stats_one_var_binary))
plot_df$var <- factor(plot_df$var, levels = vars)
plot_df$bin <- factor(plot_df$bin, levels = c("Absent (0)", "Present (1)"))
fill_colors <- c("Absent (0)" = "#FFF9C4", "Present (1)" = "#FFF176")
y_top <- 18

p <- ggplot(plot_df, aes(x = var, y = prop * 100, fill = bin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.78), width = 0.58, color = "black") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),
                position = position_dodge(width = 0.78), width = 0.45, size = 0.4, color = "black") +
  geom_hline(yintercept = 8, linetype = "dashed", color = "black", size = 0.4) +
  scale_fill_manual(values = fill_colors, name = NULL) +
  scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, 2), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "", title = NULL) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 11, color = "black"),
        legend.position = "none")

sig_df <- as.data.frame(sig, stringsAsFactors = FALSE)
sig_df$phenotype <- as.character(sig_df$phenotype)
sig_df$phenotype_clean <- gsub("_", " ", sig_df$phenotype)
sig_df$FDR <- as.numeric(sig_df$FDR)
plot_vars <- levels(plot_df$var)
sig_df <- sig_df[sig_df$phenotype_clean %in% plot_vars, ]
top_df <- aggregate(upper ~ var, data = plot_df, FUN = max, na.rm = TRUE)
top_df$y_base <- top_df$upper * 100
ann <- merge(top_df, sig_df[, c("phenotype_clean", "FDR")], by.x = "var", by.y = "phenotype_clean", all.x = FALSE, all.y = FALSE)
ann <- ann[!is.na(ann$FDR) & ann$FDR < 0.05, , drop = FALSE]

if (nrow(ann) > 0) {
  ann$xi <- as.numeric(factor(ann$var, levels = levels(plot_df$var)))
  pad1 <- 0.30; pad2 <- 0.60; pad3 <- 0.90
  ann$y_lo <- ann$y_base + pad1
  ann$y_hi <- ann$y_base + pad2
  ann$y_text <- ann$y_base + pad3
  w <- 0.25
  seg_df <- rbind(
    data.frame(x = ann$xi - w, xend = ann$xi - w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi + w, xend = ann$xi + w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi - w, xend = ann$xi + w, y = ann$y_hi, yend = ann$y_hi)
  )
  to_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  ann$label <- to_stars(ann$FDR)
  new_top <- max(y_top, max(ann$y_text + 0.25, na.rm = TRUE))
  break_top <- ceiling(new_top / 2) * 2
  p4 <- p +
    geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, size = 0.4) +
    geom_text(data = ann, aes(x = xi, y = y_text, label = label),
              inherit.aes = FALSE, size = 4) +
    scale_y_continuous(limits = c(0, new_top), breaks = seq(0, break_top, by = 2), expand = expansion(mult = c(0, 0)))
} else {
  p4 <- p
}

# perinatal complications: e
scales <- raw2[, c(72:77, 91)];scale_df <- as.data.frame(scales, stringsAsFactors = FALSE)
vars <- colnames(scale_df)[1:6];new_vars <- gsub("_", " ", vars)
colnames(scale_df)[1:6] <- new_vars;vars <- new_vars
plp_char <- as.character(scale_df[[7]]);plp_bool <- toupper(plp_char) %in% c("T", "TRUE", "1", "YES", "Y")
plot_df <- do.call(rbind, lapply(vars, bin_stats_one_var_binary))
plot_df$var <- factor(plot_df$var, levels = vars)
plot_df$bin <- factor(plot_df$bin, levels = c("Absent (0)", "Present (1)"))
fill_colors <- c("Absent (0)" = "#d9e1f2", "Present (1)" = "#b3c4e8")
y_top <- 24

p <- ggplot(plot_df, aes(x = var, y = prop * 100, fill = bin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.78), width = 0.58, color = "black") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),
                position = position_dodge(width = 0.78), width = 0.3, size = 0.4, color = "black") +
  geom_hline(yintercept = 8, linetype = "dashed", color = "black", size = 0.4) +
  scale_fill_manual(values = fill_colors, name = NULL) +
  scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, 2), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "", title = NULL) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 11, color = "black"),
        legend.position = "none")

sig_df <- as.data.frame(sig, stringsAsFactors = FALSE)
sig_df$phenotype <- as.character(sig_df$phenotype)
sig_df$phenotype_clean <- gsub("_", " ", sig_df$phenotype)
sig_df$FDR <- as.numeric(sig_df$FDR)
plot_vars <- levels(plot_df$var)
sig_df <- sig_df[sig_df$phenotype_clean %in% plot_vars, ]
top_df <- aggregate(upper ~ var, data = plot_df, FUN = max, na.rm = TRUE)
top_df$y_base <- top_df$upper * 100
ann <- merge(top_df, sig_df[, c("phenotype_clean", "FDR")], by.x = "var", by.y = "phenotype_clean", all.x = FALSE, all.y = FALSE)
ann <- ann[!is.na(ann$FDR) & ann$FDR < 0.05, , drop = FALSE]

if (nrow(ann) > 0) {
  ann$xi <- as.numeric(factor(ann$var, levels = levels(plot_df$var)))
  pad1 <- 0.30; pad2 <- 0.60; pad3 <- 1.2
  ann$y_lo <- ann$y_base + pad1
  ann$y_hi <- ann$y_base + pad2
  ann$y_text <- ann$y_base + pad3
  w <- 0.25
  seg_df <- rbind(
    data.frame(x = ann$xi - w, xend = ann$xi - w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi + w, xend = ann$xi + w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi - w, xend = ann$xi + w, y = ann$y_hi, yend = ann$y_hi)
  )
  to_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  ann$label <- to_stars(ann$FDR)
  new_top <- max(y_top, max(ann$y_text + 0.25, na.rm = TRUE))
  break_top <- ceiling(new_top / 2) * 2
  p5 <- p +
    geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, size = 0.4) +
    geom_text(data = ann, aes(x = xi, y = y_text, label = label),
              inherit.aes = FALSE, size = 4) +
    scale_y_continuous(limits = c(0, new_top), breaks = seq(0, break_top, by = 2), expand = expansion(mult = c(0, 0)))
} else {
  p5 <- p
}

# birth defects: f
scales <- raw2[, c(78:84, 91)]
scale_df <- as.data.frame(scales, stringsAsFactors = FALSE)
vars <- colnames(scale_df)[1:7]
new_vars <- gsub("_", " ", vars)
colnames(scale_df)[1:7] <- new_vars
vars <- new_vars
plp_char <- as.character(scale_df[[8]])
plp_bool <- toupper(plp_char) %in% c("T", "TRUE", "1", "YES", "Y")
plot_df <- do.call(rbind, lapply(vars, bin_stats_one_var_binary))
plot_df$var <- factor(plot_df$var, levels = vars)
plot_df$bin <- factor(plot_df$bin, levels = c("Absent (0)", "Present (1)"))
fill_colors <- c("Absent (0)" = "#dcedc8", "Present (1)" = "#aed581")
y_top <- 24

p <- ggplot(plot_df, aes(x = var, y = prop * 100, fill = bin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.78), width = 0.58, color = "black") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),
                position = position_dodge(width = 0.78), width = 0.3, size = 0.4, color = "black") +
  geom_hline(yintercept = 8, linetype = "dashed", color = "black", size = 0.4) +
  scale_fill_manual(values = fill_colors, name = NULL) +
  scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, 2), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "", title = NULL) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 11, color = "black"),
        legend.position = "none")

sig_df <- as.data.frame(sig, stringsAsFactors = FALSE)
sig_df$phenotype <- as.character(sig_df$phenotype)
sig_df$phenotype_clean <- gsub("_", " ", sig_df$phenotype)
sig_df$FDR <- as.numeric(sig_df$FDR)
plot_vars <- levels(plot_df$var)
sig_df <- sig_df[sig_df$phenotype_clean %in% plot_vars, ]
top_df <- aggregate(upper ~ var, data = plot_df, FUN = max, na.rm = TRUE)
top_df$y_base <- top_df$upper * 100
ann <- merge(top_df, sig_df[, c("phenotype_clean", "FDR")], by.x = "var", by.y = "phenotype_clean", all.x = FALSE, all.y = FALSE)
ann <- ann[!is.na(ann$FDR) & ann$FDR < 0.05, , drop = FALSE]

if (nrow(ann) > 0) {
  ann$xi <- as.numeric(factor(ann$var, levels = levels(plot_df$var)))
  pad1 <- 0.30; pad2 <- 0.60; pad3 <- 1.2
  ann$y_lo <- ann$y_base + pad1
  ann$y_hi <- ann$y_base + pad2
  ann$y_text <- ann$y_base + pad3
  w <- 0.25
  seg_df <- rbind(
    data.frame(x = ann$xi - w, xend = ann$xi - w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi + w, xend = ann$xi + w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi - w, xend = ann$xi + w, y = ann$y_hi, yend = ann$y_hi)
  )
  to_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  ann$label <- to_stars(ann$FDR)
  new_top <- max(y_top, max(ann$y_text + 0.25, na.rm = TRUE))
  break_top <- ceiling(new_top / 2) * 2
  p6 <- p +
    geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, size = 0.4) +
    geom_text(data = ann, aes(x = xi, y = y_text, label = label),
              inherit.aes = FALSE, size = 4) +
    scale_y_continuous(limits = c(0, new_top), breaks = seq(0, break_top, by = 2), expand = expansion(mult = c(0, 0)))
} else {
  p6 <- p
}

# growth abnormalities: g
scales <- raw2[, c(85:90, 91)]
scale_df <- as.data.frame(scales, stringsAsFactors = FALSE)
vars <- colnames(scale_df)[1:6]
new_vars <- gsub("_", " ", vars)
colnames(scale_df)[1:6] <- new_vars
vars <- new_vars
plp_char <- as.character(scale_df[[7]])
plp_bool <- toupper(plp_char) %in% c("T", "TRUE", "1", "YES", "Y")
plot_df <- do.call(rbind, lapply(vars, bin_stats_one_var_binary))
plot_df$var <- factor(plot_df$var, levels = vars)
plot_df$bin <- factor(plot_df$bin, levels = c("Absent (0)", "Present (1)"))
fill_colors <- c("Absent (0)" = "#f8bbd0", "Present (1)" = "#f48fb1")
y_top <- 24

p <- ggplot(plot_df, aes(x = var, y = prop * 100, fill = bin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.78), width = 0.58, color = "black") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),
                position = position_dodge(width = 0.78), width = 0.3, size = 0.4, color = "black") +
  geom_hline(yintercept = 8, linetype = "dashed", color = "black", size = 0.4) +
  scale_fill_manual(values = fill_colors, name = NULL) +
  scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, 2), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "", title = NULL) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 11, color = "black"),
        legend.position = "none")

sig_df <- as.data.frame(sig, stringsAsFactors = FALSE)
sig_df$phenotype <- as.character(sig_df$phenotype)
sig_df$phenotype_clean <- gsub("_", " ", sig_df$phenotype)
sig_df$FDR <- as.numeric(sig_df$FDR)
plot_vars <- levels(plot_df$var)
sig_df <- sig_df[sig_df$phenotype_clean %in% plot_vars, ]
top_df <- aggregate(upper ~ var, data = plot_df, FUN = max, na.rm = TRUE)
top_df$y_base <- top_df$upper * 100
ann <- merge(top_df, sig_df[, c("phenotype_clean", "FDR")], by.x = "var", by.y = "phenotype_clean", all.x = FALSE, all.y = FALSE)
ann <- ann[!is.na(ann$FDR) & ann$FDR < 0.05, , drop = FALSE]

if (nrow(ann) > 0) {
  ann$xi <- as.numeric(factor(ann$var, levels = levels(plot_df$var)))
  pad1 <- 0.30; pad2 <- 0.60; pad3 <- 0.90
  ann$y_lo <- ann$y_base + pad1
  ann$y_hi <- ann$y_base + pad2
  ann$y_text <- ann$y_base + pad3
  w <- 0.25
  seg_df <- rbind(
    data.frame(x = ann$xi - w, xend = ann$xi - w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi + w, xend = ann$xi + w, y = ann$y_lo, yend = ann$y_hi),
    data.frame(x = ann$xi - w, xend = ann$xi + w, y = ann$y_hi, yend = ann$y_hi)
  )
  to_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  ann$label <- to_stars(ann$FDR)
  new_top <- max(y_top, max(ann$y_text + 0.25, na.rm = TRUE))
  break_top <- ceiling(new_top / 2) * 2
  p7 <- p +
    geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, size = 0.4) +
    geom_text(data = ann, aes(x = xi, y = y_text, label = label),
              inherit.aes = FALSE, size = 4) +
    scale_y_continuous(limits = c(0, new_top), breaks = seq(0, break_top, by = 2), expand = expansion(mult = c(0, 0)))
} else {
  p7 <- p
}

all_aligned <- align_plots(p2, p4, p5, p6, p7, align = "v", axis = "l")
aligned_p2 <- all_aligned[[1]]
aligned_p4 <- all_aligned[[2]]
aligned_p5 <- all_aligned[[3]]
aligned_p6 <- all_aligned[[4]]
aligned_p7 <- all_aligned[[5]]

row1 <- plot_grid(p1, nrow = 1, labels = "A", label_size = 12, label_y = 0.98, vjust = 1)
row2 <- plot_grid(aligned_p2, aligned_p4, nrow = 1, rel_widths = c(1, 1), labels = c("B", "D"), label_size = 12, label_y = 0.98, vjust = 1)
row3 <- plot_grid(p3, nrow = 1, labels = "C", label_size = 12, label_y = 0.98, vjust = 1)
row4 <- plot_grid(aligned_p5, aligned_p6, aligned_p7, nrow = 1, rel_widths = c(1, 1, 1), labels = c("E", "F", "G"), label_size = 12, label_y = 0.98, vjust = 1)

final_plot <- plot_grid(row1, row2, row3, row4, ncol = 1, rel_heights = c(1, 1, 1, 1))
final_plot <- ggdraw() + draw_plot(final_plot, x = 0.055, y = 0, width = 0.89, height = 1) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("final.pdf", plot = final_plot, width = 210, height = 297, units = "mm", device = "pdf")
