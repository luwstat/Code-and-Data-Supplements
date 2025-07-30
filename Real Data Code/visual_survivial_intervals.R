####################Read data set ####################
#read data
load("C:/Users/LW HomeOffice/Dropbox/projectLeftTruncation/realdata/finalDataInfection3LT.rdata")

#number of observations per set
n <- nrow(dataInfection3LTbirth)
# the number of lifetime partners reported at the time of enrollment 
x_1 = dataInfection3LTbirth[,"nptrlife0"]
# 1: African American 0: otherwise
x_2 = as.numeric(dataInfection3LTbirth[,"race"] == 1)
C = dataInfection3LTbirth[,"agefrstsx"]

L <- dataInfection3LTbirth[,c("L_GC")]
R <- replace(dataInfection3LTbirth[,c("R_GC")],is.na(dataInfection3LTbirth[,c("R_GC")]), Inf)

d_1 <- rep(0, n)
d_2 <- rep(0, n)
d_3 <- rep(0, n)
d_0 <- as.numeric(L == R)
for (i in which(d_0 == 0)) {
  if (L[i] == C[i]) {
    d_1[i] <- 1
  }
  else if (R[i] == Inf) {
    d_3[i] <- 1
  }
  else {
    d_2[i] <- 1
  }
}

df <- cbind(L,R,C,d_0,d_1,d_2,d_3,x_1,x_2)
df <- as.data.frame(df)

# Load required libraries
library(ggplot2)
library(dplyr)

# Load your dataset
# df <- read.csv("your_data.csv")  # Uncomment if not already loaded

# Define censoring type
df <- df %>%
  mutate(
    censoring_type = case_when(
      L == R ~ "exact",
      is.infinite(R) ~ "right-censored",
      L == 0 ~ "left-censored",
      TRUE ~ "interval-censored"
    ),
    row_id = row_number()  # For y-axis plotting
  )

# Replace Inf with a large number for plotting
max_r <- max(df$R[is.finite(df$R)], na.rm = TRUE)
df$R_plot <- ifelse(is.infinite(df$R), max_r + 5, df$R)

# Plot
ggplot(df, aes(x = L, xend = R_plot, y = row_id, yend = row_id, color = censoring_type)) +
  geom_segment(size = 0.7) +
  # Black dots for exact events
  geom_point(data = filter(df, censoring_type == "exact"),
             aes(x = L, y = row_id), color = "black", size = 2) +
  # Black dots for truncation time
  geom_point(aes(x = C, y = row_id), color = "green", shape = 21, fill = "lightgreen", size = 1.8) +
  scale_color_manual(
    values = c(
      "right-censored" = "pink",
      "interval-censored" = "orchid",  # light purple
      "left-censored" = "lightblue",
      "exact" = "black"
    )
  ) +
  labs(
    x = "Time", y = "Subject",
    title = "Survival Intervals with Censoring Type and Truncation Time",
    color = "Censoring Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )
