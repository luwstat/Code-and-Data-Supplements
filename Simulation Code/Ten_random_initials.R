#install.packages("devtools") #install this pacakge if you haven’t installed it yet.
library(devtools)
# Install the proposed R pacakge regPOspline
#remotes::install_github("luwstat/regPOspline", dependencies = TRUE)
library(regPOspline)
################################################################################
#  Run po_fit() with 10 randomly generated initial_reg vectors,               #
#  extract coefficient estimates + AIC/BIC,                                   #
#  pull the (x , y) data out of each Baseline_Surv ggplot object,             #
#  build a single summary table, and                                          #
#  super-impose all baseline–survival curves in one figure.                   #
################################################################################
library(ggplot2)      # plotting
library(dplyr)        # bind_rows() convenience (optional)

## 2. Containers for results ---------------------------------------------------
results <- vector("list", 10)
summary_table  <- data.frame(
  Model = character(),
  x1_bt = numeric(),
  x2_bt = numeric(),
  AIC   = numeric(),
  BIC   = numeric(),
  stringsAsFactors = FALSE
)
baseline_list  <- list()   # will hold one data frame per model

## 3. Fit the models -----------------------------------------------------------
for (i in seq(1,10)) {
  set.seed(i) 
  initial_beta <- runif(2,-2,2)
  initial_gama <- runif(10,0,2)
  fit <- tryCatch(
    po_fit(
      L             = generated_data[, "L"],
      R             = generated_data[, "R"],
      truncation    = TRUE,
      C             = generated_data[, "C"],
      x             = generated_data[, c("x_1", "x_2")],
      initial_reg   = initial_beta,
      order         = 3,
      nknot         = 9,
      initial_spline = initial_gama
    ),
    error = function(e) {
      message("  -> ERROR: ", e$message)
      return(NULL)
    }
  )
  
  results[[i]] <- fit                          # keep the entire object
  if (is.null(fit)) {
    next                                      # skip remainder if fit failed
  } else {
    cat("Initial values:",
        "\n  Beta:", initial_beta,
        "\n  Gamma:", initial_gama, "\n")
  }                       
  
  ## ---- 3a.  Summary table row
  coef_est <- fit$coefficient_est
  summary_table <- rbind(
    summary_table,
    data.frame(
      Model = paste0("Initial_", i),
      x1_bt = coef_est["x_1", "bt"],
      x2_bt = coef_est["x_2", "bt"],
      AIC   = fit$AIC,
      BIC   = fit$BIC,
      stringsAsFactors = FALSE
    )
  )
  
  ## ---- 3b.  Pull baseline–survival data out of the ggplot
  if (inherits(fit$Baseline_Surv, "gg")) {
    build_obj   <- ggplot_build(fit$Baseline_Surv)
    
    # In case the plot has multiple layers, combine them all
    layer_dfs <- lapply(build_obj$data, function(df) {
      data.frame(Time = df$x, Survival = df$y)
    })
    baseline_df <- bind_rows(layer_dfs)
    
    baseline_df$Model <- paste0("Model_", i)   # tag curve with model id
    baseline_list[[i]] <- baseline_df
  }
}

## 4. Show / save the coefficient–AIC/BIC table --------------------------------
print(summary_table)                           # view in console
# write.csv(summary_table, "summary_table.csv", row.names = FALSE)  # if desired

## 5. Combine baseline–survival curves into one big data frame -----------------
baseline_surv_df <- bind_rows(baseline_list)   # dplyr::bind_rows keeps factors tidy

## 6. Plot all baseline–survival curves ----------------------------------------
ggplot(baseline_surv_df, aes(x = Time, y = Survival, colour = Model)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  labs(
    title = "Baseline Survival Functions Across 10 Initialisations",
    x     = "Time",
    y     = "Survival Probability"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
