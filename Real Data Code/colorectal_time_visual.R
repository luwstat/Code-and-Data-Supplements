# Load required libraries
library(ggplot2)
library(dplyr)

# Read the data
df <- read.table(file = "F:/project1/revise_realdata/C.csv", sep = "," ,header = TRUE)

# Create a censoring type column
df <- df %>%
  mutate(censoring_type = case_when(
    L == 0 & R == 0 ~ "exact",
    L == 0 & R != 0 ~ "left-censored",
    L != 0 & R == 0 ~ "right-censored",
    TRUE ~ "interval-censored"
  ))

# Assign each subject a plotting row
df$row_id <- as.integer(factor(df$SUBJID))

# Prepare a plotting dataframe
df_plot <- df %>%
  mutate(start = ifelse(censoring_type == "left-censored", 0, L),
         end = ifelse(censoring_type == "right-censored", max(t,L, R, na.rm = TRUE), 
                      ifelse(censoring_type == "exact", t, R)))

# Plot
ggplot(df_plot, aes(x = start, xend = end, y = row_id, yend = row_id, color = censoring_type)) +
  geom_segment(linewidth = 0.7) +
  geom_point(data = subset(df_plot, censoring_type == "exact"),
             aes(x = t, y = row_id), shape = 16, size = 2) +
  scale_color_manual(values = c("exact" = "black", "left-censored" = "lightblue",
                                "right-censored" = "pink", "interval-censored" = "purple")) +
  labs(x = "Time", y = "Subject", title = "Survival Time and Censoring Type Visualization") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")
