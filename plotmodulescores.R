library(ggplot2)
library(dplyr)
library(tidyr)

# Load your table
scores <- read.csv("sex_specific_tumor_signatures_sc_by_cell_type.csv", sep="\t")


# Assuming your data frame is called scores
# Pivot longer for plotting
scores_long <- scores %>%
  pivot_longer(cols = c(FemSpecTumorUp, MaleSpecTumorUp, SharedTumorUp),
               names_to = "module", values_to = "score") %>%
  mutate(group = paste(sex, condition, sep="_"))

# Generate one heatmap per cell type
plots <- scores_long %>%
  group_split(cell_type) %>%
  lapply(function(df) {
    ggplot(df, aes(x = module, y = group, fill = score)) +
      geom_tile(color = "grey70") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(title = paste("Module scores:", unique(df$cell_type)),
           x = "Module", y = "Sex_Condition") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

# Example: show first plot
print(plots[[13]])

# Save all plots to files
for (i in seq_along(plots)) {
  ggsave(filename = paste0("heatmap_", unique(scores_long$cell_type)[i], ".png"),
         plot = plots[[i]], width = 6, height = 4)
}