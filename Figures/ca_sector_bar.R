rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
})


work_dir <- "/home/mgbi/projects/gpasd/plots"

setwd(work_dir)

management <- data.frame(
  category = c("Lifestyle changes", "Clinical trial", "Procedure", "Surveillance", "Medication"),
  percent = c(19.2, 27.7, 43.2, 77.9, 87.9),
  n = c(123, 177, 276, 498, 562),
  stringsAsFactors = FALSE
)

management$category <- factor(
  management$category, 
  levels = c("Lifestyle changes", "Clinical trial", "Procedure", "Surveillance", "Medication")
)

management$n_label <- sprintf("n = %s", management$n)


positive_n <- 639

p <- ggplot(management) +
  
  geom_col(
    aes(x = percent, y = category), 
    fill = "#333333",       
    width = 0.75,            
    color = "white", linewidth = 0.5
  ) +
  
  geom_text(
    aes(x = percent, y = category, label = n_label),
    hjust = -0.2,           
    size = 3.3,
    fontface = "italic",
    color = "#333333"
  ) +
  
  annotate(
    "text", x = 0, y = 0.2, 
    label = "Positive\ndiagnosis", 
    size = 3.6, fontface = "bold", color = "#333333"
  ) +
  annotate(
    "text", x = 0, y = -0.4, 
    label = sprintf("n = %s", positive_n), 
    size = 3.3, fontface = "italic", color = "#333333"
  ) +
  
  scale_x_continuous(
    limits = c(0, 100),                
    breaks = c(0, 25, 50, 75, 100),    
    labels = c("0", "25", "50", "75", "100")
  ) +

  scale_y_discrete(expand = expansion(add = c(0.8, 0.5))) +
  
  coord_radial(
    theta = "x", 
    start=3.2*pi/2, end=0.8*pi/2,    
    inner.radius = 0.2,  
    expand = FALSE,
    clip = "off",          
    r.axis.inside = FALSE, 
    rotate.angle = FALSE   
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", color = "black", hjust = 1),
    axis.text.x = element_text(size = 9, color = "grey40"),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "#F0F0F0", color = NA),
    panel.grid.major.x = element_line(color = "white", linewidth = 0.8), 
    panel.grid.major.y = element_line(color = "white", linewidth = 0.8), 
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )

print(p)

ggsave(
  filename = "sector bar.pdf",
  plot = p,
  width = 7.5,
  height = 4.5,
  units = "in",
  limitsize = FALSE
)