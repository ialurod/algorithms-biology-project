results <- read.csv("output_file.csv", header=TRUE)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# Store your exact thresholds so we don't have to type them twice
my_thresholds <- c(2, 3, 4, 5)

# --- PLOT 1: The Accuracy Bar Chart ---
df_long <- results %>%
  # Filter out any lingering data (like 10) just in case
  filter(test_threshold %in% my_thresholds) %>%
  pivot_longer(cols = c(globin_hits, outlier_hits),
               names_to = "Hit_Type",
               values_to = "Count") %>%
  mutate(
    Hit_Type = ifelse(Hit_Type == "globin_hits", "True Alpha Globins", "False Positives (Beta/Delta)"),
    # Update to the 4 levels!
    Threshold_Label = factor(paste("Drop-off Threshold:", test_threshold), 
                             levels = paste("Drop-off Threshold:", my_thresholds))
  )

p1 <- ggplot(df_long, aes(x = factor(test_k), y = Count, fill = Hit_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 3.5, fontface = "bold", color = "black") +
  
  # Force exactly 4 columns so it locks in with the bottom graph
  facet_wrap(~ Threshold_Label, ncol = 4) + 
  
  scale_fill_manual(values = c("True Alpha Globins" = "#2C3E50", 
                               "False Positives (Beta/Delta)" = "#E74C3C")) +
  scale_y_continuous(limits = c(0, 110)) + 
  labs(
    title = "Impact of K-mer Size & Threshold on BLAST Performance",
    subtitle = "Top: Accuracy (Sensitivity vs. Specificity) | Bottom: Computational Efficiency",
    y = "Sequences Retrieved",
    fill = "Classification"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "top", 
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.title.x = element_blank(), 
    axis.text.x = element_blank()   
  )

# --- PLOT 2: The Runtime Line Chart ---
# We must apply the exact same filter and factor levels to the original results!
results_filtered <- results %>%
  filter(test_threshold %in% my_thresholds) %>%
  mutate(Threshold_Label = factor(paste("Drop-off Threshold:", test_threshold), 
                                  levels = paste("Drop-off Threshold:", my_thresholds)))

p2 <- ggplot(results_filtered, aes(x = factor(test_k), y = runtime_seconds, group = 1)) +
  geom_line(color = "#27AE60", linewidth = 1.2) + 
  geom_point(color = "#27AE60", size = 3) +
  
  # Force exactly 4 columns here too!
  facet_wrap(~ Threshold_Label, ncol = 4) + 
  
  labs(
    x = "K-mer Size (Initial Seed Length)",
    y = "Runtime (Seconds)",
    caption = "Figure 1. Trade-offs in Seed-and-Extend Alignment Parameters"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(), 
    strip.text = element_blank()        
  )

# --- COMBINE THEM WITH PATCHWORK ---
final_plot <- p1 / p2 + plot_layout(heights = c(2, 1))

# Display the plot
print(final_plot)