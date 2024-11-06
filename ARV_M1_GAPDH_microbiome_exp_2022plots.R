# Load the required libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(rstatix)  # For Wilcoxon test and pairwise comparisons
library(showtext)
library(ggsignif)  # For adding significance annotations to the plot)

# Add Times New Roman font
font_add(family = "Times New Roman", regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf")
showtext_auto()

# Set a global theme
theme_set(
  theme_linedraw(base_size = 24, base_family = "Times New Roman") +
    theme(
      axis.title = element_text(size = 24),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 24, face = "bold"),
      strip.text = element_text(size = 20),
      legend.key.size = unit(2, "lines"),
      panel.spacing = unit(1, "lines")
    )
)

# Load the data from Excel
file_path <- "ARV_micriobiome_Exp_2022_M1_GAPDH.xlsx"
jej_data <- read_excel(file_path, sheet = "Jejunum")
hrt_data <- read_excel(file_path, sheet = "Heart")

# Combine the datasets and add a new column to distinguish between them
jej_data <- jej_data %>% mutate(Tissue = "Jejunum")
hrt_data <- hrt_data %>% mutate(Tissue = "Heart")

# Filter out the negative control and select only relevant groups and numeric days
selected_groups <- c("Low AL", "Low S1133", "High AL", "High S1133")
selected_days <- c("7dpi", "21dpi")

jej_data <- jej_data %>% filter(Group %in% selected_groups, dpi %in% selected_days)
nrow(jej_data)
# filter samples where ARV M1 Ct is 46
jej_data <- jej_data %>% filter(M1_Ct < 46)
nrow(jej_data)

hrt_data <- hrt_data %>% filter(Group %in% selected_groups, dpi %in% selected_days)
nrow(hrt_data)
# filter samples where ARV M1 Ct is 46
hrt_data <- hrt_data %>% filter(M1_Ct < 46)
nrow(hrt_data)

# Set the order of the Group, Tissue, and dpi factors
combined_data <- bind_rows(jej_data, hrt_data) %>%
  mutate(Group = factor(Group, levels = c("Low AL", "Low S1133", "High AL", "High S1133")),
         dpi = factor(dpi, levels = c("7dpi", "21dpi")),
         Tissue = factor(Tissue, levels = c("Jejunum", "Heart")))

# Check if the combined data is still empty
if (nrow(combined_data) == 0) {
  stop("No data available after filtering. Check the group names and dpi values.")
}

# Perform Kruskal-Wallis Test and Wilcoxon post-hoc comparisons

# Redirect output to a text file
sink("kruskal_wallis_results.txt")

kruskal_results <- list()
wilcoxon_results <- list()

# Loop over each tissue and day combination
for(tissue in unique(combined_data$Tissue)) {
  for(day in unique(combined_data$dpi)) {
    
    # Subset the data for the specific tissue and day
    subset_data <- combined_data %>% filter(Tissue == tissue, dpi == day)
    
    # Run Kruskal-Wallis Test
    kruskal_result <- kruskal_test(subset_data, `log2^dCt` ~ Group)
    kruskal_results[[paste(tissue, day, sep = "_")]] <- kruskal_result
    # Print the Kruskal-Wallis summary
    cat("Kruskal-Wallis test for", tissue, "at dpi", day, "\n")
    print(kruskal_result)
    
    
    # If Kruskal-Wallis test is significant, run pairwise Wilcoxon test
    if (kruskal_result$p < 0.05) {
      wilcoxon_result <- pairwise_wilcox_test(subset_data, `log2^dCt` ~ Group, p.adjust.method = "none")
      wilcoxon_results[[paste(tissue, day, sep = "_")]] <- wilcoxon_result
      # Print the Dunn's test summary
      cat("Wilcoxon's test for", tissue, "at dpi", day, "\n")
      print(wilcoxon_result)
    }
    
    cat("\n")
    
  }
}
sink()

# Proceed with plotting the results
p <- ggplot(combined_data, aes(x = Group, y = `log2^dCt`, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  facet_grid(Tissue ~ dpi, scales = "free") +
  scale_fill_manual(values = c("Low AL" = "darkgoldenrod1", "Low S1133" = "#3F00FF", "High AL" = "red", "High S1133" = "purple")) +
  labs(
    y = "Relative Viral RNA (log2^dCt)",
    x = "Group",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18, family = "Times New Roman"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 18, family = "Times New Roman"),
    axis.title.y = element_text(size = 22, family = "Times New Roman"),
    legend.title = element_text(size = 22, family = "Times New Roman"),
    legend.text = element_text(size = 20, family = "Times New Roman"),
    strip.text = element_text(size = 24, face = "bold", family = "Times New Roman"),
    legend.key.size = unit(1.5, "cm"),
    panel.spacing = unit(2, "lines")
  ) +
  ylim(min(combined_data$`log2^dCt`) - 2, max(combined_data$`log2^dCt`) + 2)  # Adjust the y-axis limits

p

# Add Wilcoxon significance annotations
for (result_name in names(wilcoxon_results)) {
  wilcoxon_result <- wilcoxon_results[[result_name]]
  
  # Extract significant comparisons
  significant_pairs <- wilcoxon_result %>% filter(p.adj < 0.05) %>%
    select(group1, group2)
  
  # Add the significant pairs to the plot
  if (nrow(significant_pairs) > 0) {
    for (i in 1:nrow(significant_pairs)) {
      p <- p + stat_signif(
        comparisons = list(c(significant_pairs$group1[i], significant_pairs$group2[i])),
        map_signif_level = TRUE,
        textsize = 6,
        tip_length = 0.05,
        y_position = max(combined_data$`log2^dCt`) + 1,
        family = "Times New Roman",
        test = "wilcox.test",
        test.args = list(p.adjust.method = "none")
        
      )
    }
  }
}

p


# Save the plot
ggsave("viral_loads_comparison_boxplot_KW_with_wilcoxon_significance.svg", plot = p, width = 12, height = 8, dpi = 600)
