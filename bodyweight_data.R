# Load the necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(multcomp)  # For Tukey HSD
library(ggsignif)
library(tibble)
library(tidyr)

# Load the data
file_path <- "bw_for_R.txt"
data <- read_delim(file_path, delim = "\t")

# Remove 35 dpi data
data <- data[data$day != "35dpi", ]
unique(data$day)

# rm rows containing NAs
data <- data[complete.cases(data), ]
unique(data$Treatment)

# Order dpi to be plotted
data$day <- factor(data$day, levels = c("8 dpi", "14 dpi", "21 dpi", "28 dpi"))

# order Treatments
data$Treatment <- factor(data$Treatment, levels = c("NC", "Low AL", "High AL", 
                          "Low S1133", "High S1133"))



# Prepare to store ANOVA and Tukey results
anova_results <- list()
tukey_results <- list()

# Open a file to save the results
sink("anova_tukey_results.txt")

# Perform ANOVA and Tukey's post-hoc test for each day (dpi)
for(dpi_value in unique(data$day)) {
  # Subset the data for the current dpi
  subset_data <- data %>% filter(day == dpi_value)
  
  # Perform ANOVA
  aov_result <- aov(bw ~ Treatment, data = subset_data)
  
  # Print ANOVA results to the file
  cat("ANOVA Results for", dpi_value, "dpi:\n")
  print(summary(aov_result))
  cat("\n")
  
  # Store the ANOVA result
  anova_results[[dpi_value]] <- summary(aov_result)
  
  # If ANOVA is significant, perform Tukey's HSD post-hoc test
  if(summary(aov_result)[[1]]["Pr(>F)"][1,] < 0.05) {
    tukey_result <- TukeyHSD(aov_result)
    
    # Print Tukey HSD results to the file
    cat("Tukey HSD Post-Hoc Results for", dpi_value, "dpi:\n")
    print(tukey_result)
    cat("\n")
    
    # Store the Tukey HSD result
    tukey_results[[dpi_value]] <- tukey_result
  }
}

# Close the file connection
sink()

# Create the boxplot with significance annotations
p <- ggplot(data, aes(x = Treatment, y = bw, fill = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ day, scales = "free") +  # Facet by day (dpi) with free y-scales
  scale_fill_manual(values = c("NC" = "seagreen", "Low AL" = "darkgoldenrod1", 
                               "High AL" = "red", "Low S1133" = "#3F00FF", 
                               "High S1133" = "purple")) +  # Custom color palette
  labs(
    y = "Body Weight (g)",
    fill = "Treatment"
  ) +
  theme(
    axis.text.y = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    strip.placement = "outside",
    legend.position = "right"
  )

# 
# 
# # Add significance annotations (based on Tukey results)
# for(dpi_value in names(tukey_results)) {
#   tukey_test <- tukey_results[[dpi_value]]$Treatment
#   
#   # Subset data for this dpi to calculate y-positions
#   dpi_data <- data %>% filter(day == dpi_value)
#   
#   # Calculate max y (body weight) for this dpi
#   max_bw <- max(dpi_data$bw, na.rm = TRUE)
#   
#   # Extract significant pairs
#   significant_pairs <- as.data.frame(tukey_test) %>%
#     filter(`p adj` < 0.05) %>%
#     rownames_to_column(var = "comparison") %>%
#     separate(comparison, into = c("group1", "group2"), sep = "-")
#   
#   # Add significant annotations within each facet (dpi)
#   for (i in 1:nrow(significant_pairs)) {
#     p <- p + geom_signif(
#       comparisons = list(c(significant_pairs$group1[i], significant_pairs$group2[i])),
#       map_signif_level = TRUE,  # Adjust the y-position for each significance marker dynamically
#       tip_length = 0.01,
#       y_position = ylim()[2] - 0.1 * max_bw,  # Adjust the y-position (0.1 is a scaling factor
#       textsize = 4,
#       na.rm = TRUE
#     )
#   }
# }

# Display the plot
print(p)

# Save the plot
ggsave("body_weight_boxplot_with_significance.svg", plot = p, width = 12, height = 8, dpi = 300)
