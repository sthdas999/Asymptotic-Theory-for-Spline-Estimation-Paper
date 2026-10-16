# ================================================================
# Real-data analysis for the multivariate partially linear model
# Dataset: UCI Appliances Energy Prediction
#
# Model:
#   Y_i = A Z_i + m(X_i) + error_i
#
# Y_i = (log1p(Appliances_i), log1p(lights_i))^T
# Z_i = (sin(2*pi*hour/24), cos(2*pi*hour/24), weekend)^T
# X_i = T_out,i
#
# Observations are serially ordered at 10-minute intervals.
# HAC standard errors are practical sensitivity estimates.
# The i.i.d. asymptotic theory does not by itself justify inference
# under serial dependence.
# ================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ================================================================
# 1. READ CSV FILE
# ================================================================

csv_file <- "energydata_complete.csv"

if (!file.exists(csv_file)) {
  message("CSV file not found in the current working directory.")
  message("Current working directory: ", getwd())
  message("Please select energydata_complete.csv manually.")
  csv_file <- file.choose()
}

raw <- read.csv(csv_file, check.names = FALSE)

required <- c("date", "Appliances", "lights", "T_out")
missing_required <- setdiff(required, names(raw))

if (length(missing_required) > 0L) {
  stop(
    "Required columns not found: ",
    paste(missing_required, collapse = ", ")
  )
}

cat("\nCSV file loaded successfully.\n")
cat("Number of observations:", nrow(raw), "\n")
cat("Number of variables:", ncol(raw), "\n")
print(names(raw))

# ================================================================
# 2. CREATE DATE-TIME AND TIME-BASED PREDICTORS
# ================================================================

raw$date_time <- as.POSIXct(
  as.character(raw$date),
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)

if (anyNA(raw$date_time)) {
  raw$date_time <- as.POSIXct(
    as.character(raw$date),
    tz = "UTC"
  )
}

if (anyNA(raw$date_time)) {
  stop(
    "Unable to parse some values in the date column. ",
    "Please inspect the date format."
  )
}

hour <- as.integer(format(raw$date_time, "%H"))
weekday <- as.POSIXlt(raw$date_time)$wday

raw$hour_sin <- sin(2 * pi * hour / 24)
raw$hour_cos <- cos(2 * pi * hour / 24)
raw$weekend <- as.integer(weekday %in% c(0, 6))

cat("\nTime-based predictors created successfully.\n")

# ================================================================
# 3. DESCRIPTIVE STATISTICS
# ================================================================

summary_vars <- c(
  "Appliances", "lights", "T_out",
  "T6", "RH_out", "Windspeed"
)

summary_vars <- summary_vars[summary_vars %in% names(raw)]

descriptive <- do.call(
  rbind,
  lapply(summary_vars, function(v) {
    x <- raw[[v]]
    
    data.frame(
      variable = v,
      n = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      q1 = unname(quantile(x, 0.25, na.rm = TRUE)),
      q3 = unname(quantile(x, 0.75, na.rm = TRUE)),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  })
)

write.csv(
  descriptive,
  "descriptive_statistics.csv",
  row.names = FALSE
)

print(descriptive)

# ================================================================
# 4. CHECK WHETHER rv1 AND rv2 ARE IDENTICAL
# ================================================================

random_variable_check <- data.frame(
  item = "rv1 and rv2 identical in every row",
  result = if (all(c("rv1", "rv2") %in% names(raw))) {
    identical(raw$rv1, raw$rv2)
  } else {
    NA
  }
)

write.csv(
  random_variable_check,
  "random_variable_check.csv",
  row.names = FALSE
)

print(random_variable_check)

# ================================================================
# 5. PREPARE ANALYSIS DATA
# ================================================================

analysis <- data.frame(
  date_time = raw$date_time,
  Y_app = log1p(raw$Appliances),
  Y_lights = log1p(raw$lights),
  X = raw$T_out,
  hour_sin = raw$hour_sin,
  hour_cos = raw$hour_cos,
  weekend = raw$weekend
)

analysis <- analysis[complete.cases(analysis), ]
analysis <- analysis[order(analysis$date_time), ]
rownames(analysis) <- NULL

Y <- as.matrix(
  analysis[, c("Y_app", "Y_lights"), drop = FALSE]
)

colnames(Y) <- c(
  "log1p(Appliances)",
  "log1p(lights)"
)

Z <- as.matrix(
  analysis[, c("hour_sin", "hour_cos", "weekend"), drop = FALSE]
)

colnames(Z) <- c(
  "hour_sin",
  "hour_cos",
  "weekend"
)

X <- analysis$X
n <- nrow(analysis)

if (n < 100L) {
  stop("Too few complete observations for this analysis.")
}

cat("\nAnalysis data prepared.\n")
cat("Complete observations:", n, "\n")

# ================================================================
# 6. CHRONOLOGICAL TRAINING AND TESTING SPLIT
# ================================================================

n_train <- floor(0.80 * n)

train_id <- seq_len(n_train)
test_id <- seq.int(n_train + 1L, n)

Y_train <- Y[train_id, , drop = FALSE]
Y_test <- Y[test_id, , drop = FALSE]

Z_train <- Z[train_id, , drop = FALSE]
Z_test <- Z[test_id, , drop = FALSE]

X_train <- X[train_id]
X_test <- X[test_id]

cat("Training observations:", n_train, "\n")
cat("Holdout observations:", n - n_train, "\n")

# ================================================================
# 7. HAC LONG-RUN COVARIANCE ESTIMATOR
# ================================================================

hac_long_run <- function(G, L = 72L) {
  
  G <- as.matrix(G)
  nn <- nrow(G)
  
  if (nn < 2L) {
    stop("At least two observations are required for HAC estimation.")
  }
  
  Gc <- sweep(G, 2, colMeans(G), FUN = "-")
  S <- crossprod(Gc) / nn
  
  L <- min(as.integer(L), nn - 1L)
  
  if (L >= 1L) {
    for (ell in seq_len(L)) {
      
      Gamma_ell <- crossprod(
        Gc[(ell + 1L):nn, , drop = FALSE],
        Gc[seq_len(nn - ell), , drop = FALSE]
      ) / nn
      
      weight <- 1 - ell / (L + 1)
      
      S <- S + weight * (
        Gamma_ell + t(Gamma_ell)
      )
    }
  }
  
  (S + t(S)) / 2
}

# ================================================================
# 8. FIT THE MULTIVARIATE PARTIALLY LINEAR MODEL
# ================================================================

fit_profile <- function(spline_df = 8L, hac_lag = 72L) {
  
  spline_df <- as.integer(spline_df)
  hac_lag <- as.integer(hac_lag)
  
  n_internal <- spline_df - 4L
  
  if (n_internal < 1L) {
    stop("spline_df must be at least 5.")
  }
  
  # Use the training covariate distribution to determine knots.
  # Boundary knots cover the observed covariate support.
  boundary_knots <- range(X)
  
  probs <- seq(
    0, 1,
    length.out = n_internal + 2L
  )[-c(1L, n_internal + 2L)]
  
  internal_knots <- as.numeric(
    quantile(
      X_train,
      probs = probs,
      names = FALSE,
      type = 7
    )
  )
  
  internal_knots <- sort(unique(internal_knots))
  
  if (length(internal_knots) != n_internal) {
    stop("Repeated internal knots occurred. Try a smaller spline_df.")
  }
  
  # Construct a common cubic B-spline basis.
  B_train <- splines::bs(
    X_train,
    knots = internal_knots,
    degree = 3,
    intercept = TRUE,
    Boundary.knots = boundary_knots
  )
  
  B_test <- splines::bs(
    X_test,
    knots = internal_knots,
    degree = 3,
    intercept = TRUE,
    Boundary.knots = boundary_knots
  )
  
  # Estimate E(Y | X) and E(Z | X) on the training sample.
  coef_muY <- qr.solve(B_train, Y_train)
  coef_muZ <- qr.solve(B_train, Z_train)
  
  muY_train <- B_train %*% coef_muY
  muZ_train <- B_train %*% coef_muZ
  
  Y_tilde <- Y_train - muY_train
  Z_tilde <- Z_train - muZ_train
  
  # Estimate the coefficient matrix A.
  Omega_hat <- crossprod(Z_tilde) / nrow(Z_tilde)
  
  if (rcond(Omega_hat) < 1e-10) {
    warning(
      "Residualized Z Gram matrix may be ill-conditioned."
    )
  }
  
  A_hat <- t(
    solve(
      crossprod(Z_tilde),
      crossprod(Z_tilde, Y_tilde)
    )
  )
  
  colnames(A_hat) <- colnames(Z_train)
  rownames(A_hat) <- colnames(Y_train)
  
  # Residuals from the partialled-out regression.
  eps_hat <- Y_tilde - Z_tilde %*% t(A_hat)
  
  q <- ncol(Y_train)
  d <- ncol(Z_train)
  
  # Construct score vectors vec(epsilon_i z_tilde_i^T).
  G <- t(
    vapply(
      seq_len(nrow(eps_hat)),
      function(i) {
        as.vector(
          outer(eps_hat[i, ], Z_tilde[i, ])
        )
      },
      numeric(q * d)
    )
  )
  
  # HAC covariance estimate.
  S_hac <- hac_long_run(G, L = hac_lag)
  
  Omega_inv <- solve(Omega_hat)
  bread <- kronecker(Omega_inv, diag(q))
  
  vcov_vecA <- bread %*% S_hac %*% bread /
    nrow(Y_train)
  
  vcov_vecA <- (vcov_vecA + t(vcov_vecA)) / 2
  
  se_vecA <- sqrt(
    pmax(diag(vcov_vecA), 0)
  )
  
  # Coefficient table: column-major ordering of vec(A).
  coef_rows <- list()
  ii <- 1L
  
  for (j in seq_len(d)) {
    for (k in seq_len(q)) {
      
      idx <- (j - 1L) * q + k
      
      est <- A_hat[k, j]
      se <- se_vecA[idx]
      
      coef_rows[[ii]] <- data.frame(
        predictor = colnames(Z_train)[j],
        outcome = colnames(Y_train)[k],
        estimate = est,
        HAC_SE = se,
        CI_lower = est - qnorm(0.975) * se,
        CI_upper = est + qnorm(0.975) * se,
        z_value = if (se > 0) est / se else NA_real_,
        p_value = if (se > 0) {
          2 * pnorm(-abs(est / se))
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
      
      ii <- ii + 1L
    }
  }
  
  coef_table <- do.call(rbind, coef_rows)
  
  # Estimate the nonparametric component:
  # m(x) = E(Y | X=x) - A E(Z | X=x).
  muY_test <- B_test %*% coef_muY
  muZ_test <- B_test %*% coef_muZ
  
  m_test <- muY_test - muZ_test %*% t(A_hat)
  
  Yhat_test <- Z_test %*% t(A_hat) + m_test
  
  # Training fitted values and residuals.
  Yhat_train <- Z_train %*% t(A_hat) +
    muY_train - muZ_train %*% t(A_hat)
  
  resid_train <- Y_train - Yhat_train
  
  # Holdout prediction metrics on the log1p scale.
  metrics <- do.call(
    rbind,
    lapply(seq_len(ncol(Y_test)), function(k) {
      
      err <- Y_test[, k] - Yhat_test[, k]
      
      data.frame(
        outcome = colnames(Y_test)[k],
        n_test = length(err),
        RMSE = sqrt(mean(err^2)),
        MAE = mean(abs(err)),
        R_squared = if (
          sum((Y_test[, k] - mean(Y_test[, k]))^2) > 0
        ) {
          1 - sum(err^2) /
            sum((Y_test[, k] - mean(Y_test[, k]))^2)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    })
  )
  
  # Ljung-Box residual serial-dependence diagnostics.
  diagnostic_rows <- list()
  di <- 1L
  
  diag_lags <- c(24L, 72L, 144L)
  
  for (k in seq_len(ncol(resid_train))) {
    for (ell in diag_lags) {
      
      bt <- Box.test(
        resid_train[, k],
        lag = ell,
        type = "Ljung-Box"
      )
      
      diagnostic_rows[[di]] <- data.frame(
        outcome = colnames(Y_train)[k],
        lag = ell,
        Ljung_Box_statistic = unname(bt$statistic),
        p_value = bt$p.value,
        stringsAsFactors = FALSE
      )
      
      di <- di + 1L
    }
  }
  
  diagnostics <- do.call(rbind, diagnostic_rows)
  
  # Smooth components over the observed temperature range.
  x_grid <- seq(
    min(X),
    max(X),
    length.out = 250L
  )
  
  B_grid <- splines::bs(
    x_grid,
    knots = internal_knots,
    degree = 3,
    intercept = TRUE,
    Boundary.knots = boundary_knots
  )
  
  muY_grid <- B_grid %*% coef_muY
  muZ_grid <- B_grid %*% coef_muZ
  
  m_grid <- muY_grid - muZ_grid %*% t(A_hat)
  
  smooth_grid <- data.frame(
    T_out = x_grid,
    m_log1p_Appliances = m_grid[, 1],
    m_log1p_lights = m_grid[, 2]
  )
  
  list(
    A = A_hat,
    vcov = vcov_vecA,
    coef_table = coef_table,
    metrics = metrics,
    diagnostics = diagnostics,
    residuals = resid_train,
    smooth_grid = smooth_grid,
    spline_df = spline_df,
    hac_lag = hac_lag,
    train_n = length(train_id),
    test_n = length(test_id),
    train_end = analysis$date_time[n_train],
    test_start = analysis$date_time[n_train + 1L],
    Yhat_test = Yhat_test,
    Yhat_train = Yhat_train,
    internal_knots = internal_knots,
    boundary_knots = boundary_knots
  )
}

# ================================================================
# 9. FIT THE MAIN MODEL
# ================================================================

main_fit <- fit_profile(
  spline_df = 8L,
  hac_lag = 72L
)

write.csv(
  main_fit$coef_table,
  "profile_coefficients_HAC.csv",
  row.names = FALSE
)

write.csv(
  main_fit$metrics,
  "holdout_prediction_metrics.csv",
  row.names = FALSE
)

write.csv(
  main_fit$diagnostics,
  "residual_serial_dependence_diagnostics.csv",
  row.names = FALSE
)

write.csv(
  main_fit$smooth_grid,
  "estimated_smooth_component.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(main_fit$A),
  "estimated_A_matrix.csv"
)

cat("\nEstimated coefficient matrix A:\n")
print(main_fit$A)

cat("\nHoldout prediction metrics:\n")
print(main_fit$metrics)

cat("\nHAC coefficient estimates:\n")
print(main_fit$coef_table)

# ================================================================
# 10. SPLINE-DIMENSION SENSITIVITY
# ================================================================

spline_dfs <- c(6L, 8L, 10L)

spline_sensitivity <- do.call(
  rbind,
  lapply(spline_dfs, function(df) {
    
    obj <- fit_profile(
      spline_df = df,
      hac_lag = 72L
    )
    
    cbind(
      data.frame(spline_df = df),
      obj$metrics
    )
  })
)

write.csv(
  spline_sensitivity,
  "spline_dimension_sensitivity.csv",
  row.names = FALSE
)

print(spline_sensitivity)

# ================================================================
# 11. HAC BANDWIDTH SENSITIVITY
# ================================================================

hac_lags <- c(12L, 36L, 72L, 144L)

hac_sensitivity <- do.call(
  rbind,
  lapply(hac_lags, function(L) {
    
    obj <- fit_profile(
      spline_df = 8L,
      hac_lag = L
    )
    
    cbind(
      data.frame(hac_lag = L),
      obj$coef_table
    )
  })
)

write.csv(
  hac_sensitivity,
  "HAC_bandwidth_sensitivity.csv",
  row.names = FALSE
)

# ================================================================
# 12. RESPONSE AND RESIDUAL CORRELATIONS
# ================================================================

correlations <- data.frame(
  quantity = c(
    "Correlation of log1p outcomes",
    "Correlation of training residuals"
  ),
  correlation = c(
    cor(Y[, 1], Y[, 2]),
    cor(
      main_fit$residuals[, 1],
      main_fit$residuals[, 2]
    )
  )
)

write.csv(
  correlations,
  "response_and_residual_correlations.csv",
  row.names = FALSE
)

print(correlations)

# ================================================================
# 13. PLOT THE FIRST OBSERVED WEEK
# ================================================================

png(
  "appliances_energy_time_series_example.png",
  width = 1200,
  height = 700,
  res = 120
)

week_n <- min(
  7L * 24L * 6L,
  nrow(analysis)
)

plot(
  analysis$date_time[seq_len(week_n)],
  Y[seq_len(week_n), 1],
  type = "l",
  xlab = "Time",
  ylab = "log(1 + Appliances energy, Wh)",
  main = "Appliances energy use: first observed week"
)

dev.off()

# ================================================================
# 14. PLOT ESTIMATED SMOOTH TEMPERATURE COMPONENTS
# ================================================================

png(
  "estimated_smooth_temperature_effects.png",
  width = 1200,
  height = 700,
  res = 120
)

plot(
  main_fit$smooth_grid$T_out,
  main_fit$smooth_grid$m_log1p_Appliances,
  type = "l",
  lwd = 2,
  xlab = "Outdoor temperature T_out (degrees C)",
  ylab = "Estimated smooth component",
  main = "Estimated smooth temperature components"
)

lines(
  main_fit$smooth_grid$T_out,
  main_fit$smooth_grid$m_log1p_lights,
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = c(
    "log1p(Appliances)",
    "log1p(lights)"
  ),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

dev.off()

# ================================================================
# 15. PLOT HOLDOUT ACTUAL VERSUS PREDICTED VALUES
# ================================================================

png(
  "holdout_actual_vs_predicted.png",
  width = 1200,
  height = 700,
  res = 120
)

plot(
  Y_test[, 1],
  main_fit$Yhat_test[, 1],
  pch = 16,
  cex = 0.45,
  xlab = "Observed log1p(Appliances)",
  ylab = "Predicted log1p(Appliances)",
  main = "Chronological holdout: observed versus predicted"
)

abline(
  a = 0,
  b = 1,
  lty = 2,
  lwd = 2
)

dev.off()

# ================================================================
# 16. FINAL SUMMARY
# ================================================================

cat("\n================================================\n")
cat("ANALYSIS COMPLETED\n")
cat("================================================\n")

cat(
  "Output CSV and PNG files saved in:\n",
  getwd(), "\n"
)

cat("Rows used:", n, "\n")
cat("Training rows:", n_train, "\n")
cat("Holdout rows:", n - n_train, "\n")

cat(
  "Training period ends:",
  format(main_fit$train_end),
  "\n"
)

cat(
  "Holdout period starts:",
  format(main_fit$test_start),
  "\n"
)

cat(
  "Main model: cubic B-spline df =",
  main_fit$spline_df,
  "; HAC lag =",
  main_fit$hac_lag,
  "\n"
)

cat("\nGenerated output files:\n")

print(c(
  "descriptive_statistics.csv",
  "random_variable_check.csv",
  "profile_coefficients_HAC.csv",
  "holdout_prediction_metrics.csv",
  "residual_serial_dependence_diagnostics.csv",
  "estimated_smooth_component.csv",
  "estimated_A_matrix.csv",
  "spline_dimension_sensitivity.csv",
  "HAC_bandwidth_sensitivity.csv",
  "response_and_residual_correlations.csv",
  "appliances_energy_time_series_example.png",
  "estimated_smooth_temperature_effects.png",
  "holdout_actual_vs_predicted.png"
))