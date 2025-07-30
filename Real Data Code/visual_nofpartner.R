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
# Basic summary statistics
summary(df$x_1)
sd(df$x_1, na.rm = TRUE)  # Standard deviation
length(unique(df$x_1))    # Number of unique values

# Table of frequencies (if it's discrete/small range)
table(df$x_1)

# Save to PDF or run in enlarged plotting window if needed
# pdf("x1_hist_density_boxplot_tight.pdf", width = 7, height = 5)

layout(matrix(1:2, ncol = 1), heights = c(4, 1))  # Make boxplot panel smaller

# Reduce outer and inner margins
par(mar = c(0, 4, 4, 2))  # top plot: no bottom margin, keep labels

# --- Histogram + density ---
hist(df$x_1,
     breaks = 20,
     col = "lightblue",
     border = "white",
     probability = TRUE,
     main = "",
     xlab = "",
     ylab = "Density")

lines(density(df$x_1, na.rm = TRUE),
      col = "darkblue", lwd = 2)

# Now boxplot
par(mar = c(4, 4, 0, 2))  # bottom plot: no top margin

boxplot(df$x_1,
        horizontal = TRUE,
        col = "orange",
        axes = FALSE,
        ylim = range(df$x_1, na.rm = TRUE))

# dev.off()  # if saving


# Reset layout (e.g., mfrow or mfcol) back to 1 plot
par(mfrow = c(1, 1))

# Reset margins to default (bottom, left, top, right)
par(mar = c(5.1, 4.1, 4.1, 2.1))
# View basic information
table(df$x_2)                            # Frequency table
prop.table(table(df$x_2))               # Proportion table
summary(as.factor(df$x_2))              # Summary as factor (for labeled display)
length(unique(df$x_2))                  # Number of unique race categories

# Optional: convert to labeled factor if needed
# df$x_2 <- factor(df$x_2, labels = c("White", "Black", "Asian", "Other"))

# Bar plot
barplot(table(df$x_2), 
        col = "skyblue", 
        main = "Barplot of Race (x_2)", 
        xlab = "Race Categories", 
        ylab = "Count")

# Pie chart (if appropriate)
pie(table(df$x_2), 
    col = rainbow(length(unique(df$x_2))), 
    main = "Race Distribution (x_2)")
