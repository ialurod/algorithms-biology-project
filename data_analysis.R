results <- read.csv("output_file.csv", header=TRUE)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(ggrepel)

my_thresholds <- c(10, 13, 16, 19)

# --- PLOT 1: The Accuracy Bar Chart ---
df_long <- results %>%
  filter(test_threshold %in% my_thresholds) %>%
  filter(grepl("P69905", query_id)) %>% 
  pivot_longer(cols = c(globin_hits, outlier_hits),
               names_to = "Hit_Type",
               values_to = "Count") %>%
  mutate(
    Hit_Type = ifelse(Hit_Type == "globin_hits", "True Alpha Globins", "False Positives (Beta/Delta)"),
    Threshold_Label = factor(paste("Seed Threshold:", test_threshold), 
                             levels = paste("Seed Threshold:", my_thresholds))
  )

p1 <- ggplot(df_long, aes(x = factor(test_k), y = Count, fill = Hit_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 3.5, fontface = "bold", color = "black") +
  facet_wrap(~ Threshold_Label, ncol = 4) + 
  scale_fill_manual(values = c("True Alpha Globins" = "#2C3E50", 
                               "False Positives (Beta/Delta)" = "#E74C3C")) +
  scale_y_continuous(limits = c(0, 110)) + 
  labs(
    title = "Impact of K-mer Size & BLOSUM Seed Threshold on BLAST Performance",
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
results_filtered <- results %>%
  filter(test_threshold %in% my_thresholds) %>%
  filter(grepl("P69905", query_id)) %>% 
  mutate(Threshold_Label = factor(paste("Seed Threshold:", test_threshold), 
                                  levels = paste("Seed Threshold:", my_thresholds)))

p2 <- ggplot(results_filtered, aes(x = factor(test_k), y = runtime_seconds, group = 1)) +
  geom_line(color = "#27AE60", linewidth = 1.2) + 
  geom_point(color = "#27AE60", size = 3) +
  facet_wrap(~ Threshold_Label, ncol = 4) + 
  labs(
    x = "K-mer Size (Initial Seed Length)",
    y = "Runtime (Seconds)",
    caption = "Figure 1. Trade-offs in BLOSUM62 Seed Alignment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_blank(), 
    strip.text = element_blank()        
  )

final_plot <- p1 / p2 + plot_layout(heights = c(2, 1))
print(final_plot)

# --- PLOTS DIVERGENTES ---
df_divergent <- results %>%
  filter(!grepl("P69905", query_id)) %>% 
  filter(test_threshold %in% my_thresholds) %>%
  mutate(
    Threshold_Label = factor(paste("Seed Threshold:", test_threshold), levels = paste("Seed Threshold:", my_thresholds))
  )

p3 <- ggplot(df_divergent, aes(x = factor(test_k), y = globin_hits, fill = query_id)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Threshold_Label, ncol = 4) + 
  scale_fill_viridis_d(option = "viridis", begin = 0.2, end = 0.9) + 
  labs(title = "Divergence Analysis: BLOSUM Seed Accuracy", y = "True Alpha Globin Hits", fill = "Query Divergence") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16), axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "top", strip.background = element_rect(fill = "gray90", color = NA), strip.text = element_text(face = "bold", size = 12))

p4 <- ggplot(df_divergent, aes(x = factor(test_k), y = runtime_seconds, color = query_id, group = query_id)) +
  geom_line(linewidth = 1.2) + 
  geom_point(size = 3) +
  facet_wrap(~ Threshold_Label, ncol = 4) + 
  scale_color_viridis_d(option = "viridis", begin = 0.2, end = 0.9) + 
  labs(x = "K-mer Size", y = "Runtime (Seconds)", color = "Query Divergence") +
  theme_minimal(base_size = 14) +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = "none")

divergent_plot <- p3 / p4 + plot_layout(heights = c(2, 1))
print(divergent_plot)

# --- PLOT GLOBAL TRADE-OFF ---
df_global <- results %>%
  mutate(
    Divergence = factor(case_when(
      grepl("P69905", query_id) ~ "0% (Baseline)",
      grepl("10", query_id) ~ "10% Mutation",
      grepl("20", query_id) ~ "20% Mutation",
      grepl("30", query_id) ~ "30% Mutation"
    ), levels = c("0% (Baseline)", "10% Mutation", "20% Mutation", "30% Mutation")),
    Is_Optimal = ifelse(test_k == 4 & test_threshold == 16, "Optimal Balance (K=4, Seed=16)", "Other Configurations")
  )

p_global <- ggplot(df_global, aes(x = runtime_seconds, y = globin_hits, color = Divergence)) +
  geom_point(aes(shape = Is_Optimal, size = Is_Optimal), alpha = 0.85) +
  geom_text_repel(data = filter(df_global, Is_Optimal == "Optimal Balance (K=4, Seed=16)"),
                  aes(label = paste(globin_hits, "hits\n(", round(runtime_seconds, 2), "s)", sep="")),
                  show.legend = FALSE, size = 3.5, fontface = "bold", 
                  box.padding = 0.8, segment.color = "gray50") +
  scale_color_viridis_d(option = "viridis", begin = 0.1, end = 0.9) +
  scale_shape_manual(values = c("Optimal Balance (K=4, Seed=16)" = 17, "Other Configurations" = 16)) + 
  scale_size_manual(values = c("Optimal Balance (K=4, Seed=16)" = 5, "Other Configurations" = 3)) +
  labs(
    title = "Global Performance & BLOSUM Parameter Balance",
    x = "Runtime (Seconds)",
    y = "Accuracy (Recovered Globin Hits %)",
    color = "Query Divergence",
    shape = "Parameter Profile",
    size = "Parameter Profile"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face="bold", size=16, color = "#2C3E50"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90")
  )

print(p_global)
ggsave("plot_global_balance.png", plot = p_global, width = 12, height = 7, dpi = 300)