
library(ggplot2)
library(tidyr)

# Creating the data
data <- data.frame(
  Analysis = c("DEG", "DET", "DET-G", "DUT", "DUT-G"),
  ROSMAP = c(995, 1398, 992, 293, 185),
  MSBB = c(320, 2450, 1979, 541, 425)
)

# Convert the data to long format
data_long <- pivot_longer(data, cols = c("ROSMAP", "MSBB"), names_to = "Dataset", values_to = "Count")

data_long$Dataset <- factor(data_long$Dataset, levels = c("ROSMAP", "MSBB"))


# Color palette
color_palette <- c(
  "DEG" = "#FF7F7F",    # Red
  "DET" = "#ADD8E6",    # Light blue
  "DET-G" = "#20B2AA",  # Turquoise
  "DUT" = "#9370DB",    # Purple
  "DUT-G" = "#647CB0"   # Blue
)


ggplot(data_long, aes(x = Dataset, y = Count, fill = Analysis)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.8), 
            vjust = -0.3, size = 10) +  # Adding labels to the bars
  scale_fill_manual(values = color_palette) +
  labs(
    title = "Comparison of Analyses Across ROSMAP and MSBB",
    fill = "Analysis"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

path <- "~/Desktop/Figürler"

ggsave(filename = file.path(path, "comparison_plot.tiff"),
       plot = last_plot(),
       width = 12, height = 10, device='tiff', dpi=600)

