# ================================================================
# Real-data analysis for the multivariate partially linear model
# Dataset: UCI Appliances Energy Prediction (energydata_complete.csv)
# Model: log1p(Y_i) = A Z_i + m(T_out,i) + error_i
# Y_i = (Appliances_i, lights_i)^T
# Z_i = (sin(2*pi*hour/24), cos(2*pi*hour/24), weekend)^T
# X_i = T_out,i (smooth, univariate covariate)
#
# Important: observations are serially ordered at 10-minute intervals.
# HAC standard errors below are practical sensitivity estimates; the
# manuscript's i.i.d. asymptotic theory does not by itself justify them.
# Run this script from the directory containing the supplied ZIP file.
# ================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

zip_file <- "appliances+energy+prediction.zip"
if (!file.exists(zip_file)) {
  stop("Cannot find 'appliances+energy+prediction.zip'. Set zip_file to its full path.")
}

csv_name <- "energydata_complete.csv"
zip_contents <- unzip(zip_file, list = TRUE)$Name
if (!(csv_name %in% zip_contents)) {
  stop("The ZIP archive does not contain energydata_complete.csv.")
}
raw <- read.csv(unz(zip_file, csv_name), check.names = FALSE)

required <- c("date", "Appliances", "lights", "T_out")
missing_required <- setdiff(required, names(raw))
if (length(missing_required) > 0) {
  stop("Required columns not found: ", paste(missing_required, collapse = ", "))
}

# Parse timestamps and construct a chronological dataset.
raw$date_time <- as.POSIXct(raw$date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
if (anyNA(raw$date_time)) stop("Some timestamps could not be parsed; inspect raw$date.")
raw <- raw[order(raw$date_time), ]
rownames(raw) <- NULL

# This application focuses on two energy outcomes and one smooth temperature effect.
# Log1p is used to reduce right-skew and to accommodate zero lighting readings.
# Coefficients and prediction errors are therefore on the log(1 + Wh) scale.
if (any(raw$Appliances < 0, na.rm = TRUE) || any(raw$lights < 0, na.rm = TRUE)) {
  stop("Energy outcomes contain negative values; log1p transformation is not appropriate.")
}

raw$hour <- as.POSIXlt(raw$date_time, tz = "UTC")$hour
raw$wday <- as.POSIXlt(raw$date_time, tz = "UTC")$wday  # Sunday=0
raw$hour_sin <- sin(2 * pi * raw$hour / 24)
raw$hour_cos <- cos(2 * pi * raw$hour / 24)
raw$weekend <- as.numeric(raw$wday %in% c(0, 6))

# Check key data integrity features and report them to disk.
num_cols <- names(raw)[vapply(raw, is.numeric, logical(1))]
missing_by_column <- data.frame(
  variable = names(raw),
  missing_n = vapply(raw, function(z) sum(is.na(z)), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(missing_by_column, "missing_values_by_column.csv", row.names = FALSE)

profile <- data.frame(
  item = c("Number of rows", "Number of original columns", "First timestamp",
           "Last timestamp", "Number of missing cells in original data",
           "Number of duplicated full rows", "Number of duplicated timestamps",
           "Number of unique timestamps", "Median interval in minutes",
           "Number of distinct timestamp intervals in minutes"),
  value = c(
    nrow(raw), ncol(raw) - 6,
    format(min(raw$date_time), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    format(max(raw$date_time), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    sum(is.na(raw[, setdiff(names(raw), c("date_time", "hour", "wday", "hour_sin", "hour_cos", "weekend")), drop = FALSE])),
    sum(duplicated(raw[, setdiff(names(raw), c("date_time", "hour", "wday", "hour_sin", "hour_cos", "weekend")), drop = FALSE])),
    sum(duplicated(raw$date_time)), length(unique(raw$date_time)),
    if (nrow(raw) > 1) median(as.numeric(diff(raw$date_time), units = "mins")) else NA,
    if (nrow(raw) > 1) length(unique(as.numeric(diff(raw$date_time), units = "mins"))) else NA
  ), stringsAsFactors = FALSE
)
write.csv(profile, "dataset_profile.csv", row.names = FALSE)

# Descriptive statistics on the original Wh scale; all values are computed from the ZIP.
summary_vars <- c("Appliances", "lights", "T_out", "T6", "RH_out", "Windspeed")
summary_vars <- summary_vars[summary_vars %in% names(raw)]
descriptive <- do.call(rbind, lapply(summary_vars, function(v) {
  x <- raw[[v]]
  data.frame(variable = v, n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE),
             sd = sd(x, na.rm = TRUE), median = median(x, na.rm = TRUE),
             q1 = unname(quantile(x, .25, na.rm = TRUE)),
             q3 = unname(quantile(x, .75, na.rm = TRUE)),
             min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
}))
write.csv(descriptive, "descriptive_statistics.csv", row.names = FALSE)

# Exact duplicate check for the two synthetic/random columns, if present.
random_variable_check <- data.frame(
  item = "rv1 and rv2 identical in every row",
  result = if (all(c("rv1", "rv2") %in% names(raw))) identical(raw$rv1, raw$rv2) else NA
)
write.csv(random_variable_check, "random_variable_check.csv", row.names = FALSE)

# Analysis variables. Rows with missing analysis variables are excluded explicitly.
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

Y <- as.matrix(analysis[, c("Y_app", "Y_lights")])
colnames(Y) <- c("log1p(Appliances)", "log1p(lights)")
Z <- as.matrix(analysis[, c("hour_sin", "hour_cos", "weekend")])
colnames(Z) <- c("hour_sin", "hour_cos", "weekend")
X <- analysis$X
n <- nrow(analysis)
if (n < 100) stop("Too few complete observations for this analysis.")

# Chronological split: the earliest 80% is used for estimation; the latest 20%
# is held out for a future-period evaluation. No random shuffling is performed.
n_train <- floor(0.80 * n)
train_id <- seq_len(n_train)
test_id <- seq.int(n_train + 1L, n)
Y_train <- Y[train_id, , drop = FALSE]
Y_test <- Y[test_id, , drop = FALSE]
Z_train <- Z[train_id, , drop = FALSE]
Z_test <- Z[test_id, , drop = FALSE]
X_train <- X[train_id]
X_test <- X[test_id]

# HAC estimator for the long-run covariance of the score sequence.
# Bartlett/Newey-West weights: 1 - lag/(L+1).
hac_long_run <- function(G, L = 72L) {
  G <- as.matrix(G)
  nn <- nrow(G)
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
      S <- S + weight * (Gamma_ell + t(Gamma_ell))
    }
  }
  (S + t(S)) / 2
}

fit_profile <- function(spline_df = 8L, hac_lag = 72L) {
  # Use a fixed cubic B-spline basis; df is prespecified, not selected adaptively.
  # For a cubic basis, the number of interior knots is df - 4.
  spline_df <- as.integer(spline_df)
  n_internal <- spline_df - 4L
  if (n_internal < 1L) stop("spline_df must be at least 5.")
  boundary_knots <- range(X)  # covariate support; response values are not used
  probs <- seq(0, 1, length.out = n_internal + 2L)[-c(1L, n_internal + 2L)]
  internal_knots <- as.numeric(quantile(X_train, probs = probs, names = FALSE, type = 7))
  internal_knots <- sort(unique(internal_knots))
  if (length(internal_knots) != n_internal) {
    stop("Repeated internal knots occurred. Try a smaller spline_df.")
  }

  B_train <- splines::bs(X_train, knots = internal_knots, degree = 3,
                         intercept = TRUE, Boundary.knots = boundary_knots)
  B_test <- splines::bs(X_test, knots = internal_knots, degree = 3,
                        intercept = TRUE, Boundary.knots = boundary_knots)
  B_all <- splines::bs(X, knots = internal_knots, degree = 3,
                       intercept = TRUE, Boundary.knots = boundary_knots)

  # Spline least-squares estimates of E(Y | X) and E(Z | X).
  coef_muY <- qr.solve(B_train, Y_train)
  coef_muZ <- qr.solve(B_train, Z_train)
  muY_train <- B_train %*% coef_muY
  muZ_train <- B_train %*% coef_muZ
  Y_tilde <- Y_train - muY_train
  Z_tilde <- Z_train - muZ_train

  Omega_hat <- crossprod(Z_tilde) / nrow(Z_tilde)
  if (rcond(Omega_hat) < 1e-10) warning("Residualized Z Gram matrix may be ill-conditioned.")
  A_hat <- t(solve(crossprod(Z_tilde), crossprod(Z_tilde, Y_tilde)))
  colnames(A_hat) <- colnames(Z_train)
  rownames(A_hat) <- colnames(Y_train)

  eps_hat <- Y_tilde - Z_tilde %*% t(A_hat)
  q <- ncol(Y_train)
  d <- ncol(Z_train)
  # vec(epsilon_i ztilde_i^T), column-major, has dimension q*d.
  G <- t(vapply(seq_len(nrow(eps_hat)), function(i) {
    as.vector(outer(eps_hat[i, ], Z_tilde[i, ]))
  }, numeric(q * d)))
  S_hac <- hac_long_run(G, L = hac_lag)
  Omega_inv <- solve(Omega_hat)
  bread <- kronecker(Omega_inv, diag(q))
  vcov_vecA <- bread %*% S_hac %*% bread / nrow(Y_train)
  vcov_vecA <- (vcov_vecA + t(vcov_vecA)) / 2
  se_vecA <- sqrt(pmax(diag(vcov_vecA), 0))

  # Coefficient table follows R's column-major vec(A) ordering.
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
        estimate = est, HAC_SE = se,
        CI_lower = est - qnorm(.975) * se,
        CI_upper = est + qnorm(.975) * se,
        z_value = ifelse(se > 0, est / se, NA_real_),
        p_value = ifelse(se > 0, 2 * pnorm(-abs(est / se)), NA_real_),
        stringsAsFactors = FALSE
      )
      ii <- ii + 1L
    }
  }
  coef_table <- do.call(rbind, coef_rows)

  # The estimated nonparametric component is m(x) = mu_Y(x) - A mu_Z(x).
  muY_test <- B_test %*% coef_muY
  muZ_test <- B_test %*% coef_muZ
  m_test <- muY_test - muZ_test %*% t(A_hat)
  Yhat_test <- Z_test %*% t(A_hat) + m_test

  # Training residuals and time-ordered holdout metrics on the log1p scale.
  Yhat_train <- Z_train %*% t(A_hat) + muY_train - muZ_train %*% t(A_hat)
  resid_train <- Y_train - Yhat_train
  metrics <- do.call(rbind, lapply(seq_len(ncol(Y_test)), function(k) {
    err <- Y_test[, k] - Yhat_test[, k]
    data.frame(
      outcome = colnames(Y_test)[k], n_test = length(err),
      RMSE = sqrt(mean(err^2)), MAE = mean(abs(err)),
      R_squared = 1 - sum(err^2) / sum((Y_test[, k] - mean(Y_test[, k]))^2),
      stringsAsFactors = FALSE
    )
  }))

  # Residual serial-dependence diagnostics (diagnostic, not a proof of validity).
  diagnostic_rows <- list()
  di <- 1L
  diag_lags <- c(24L, 72L, 144L)
  for (k in seq_len(ncol(resid_train))) {
    for (ell in diag_lags) {
      bt <- Box.test(resid_train[, k], lag = ell, type = "Ljung-Box")
      diagnostic_rows[[di]] <- data.frame(
        outcome = colnames(Y_train)[k], lag = ell,
        Ljung_Box_statistic = unname(bt$statistic), p_value = bt$p.value,
        stringsAsFactors = FALSE
      )
      di <- di + 1L
    }
  }
  diagnostics <- do.call(rbind, diagnostic_rows)

  # Smooth component evaluated over a grid, for plots and inspection.
  x_grid <- seq(min(X), max(X), length.out = 250L)
  B_grid <- splines::bs(x_grid, knots = internal_knots, degree = 3,
                        intercept = TRUE, Boundary.knots = boundary_knots)
  muY_grid <- B_grid %*% coef_muY
  muZ_grid <- B_grid %*% coef_muZ
  m_grid <- muY_grid - muZ_grid %*% t(A_hat)
  smooth_grid <- data.frame(T_out = x_grid,
                            m_log1p_Appliances = m_grid[, 1],
                            m_log1p_lights = m_grid[, 2])

  list(A = A_hat, vcov = vcov_vecA, coef_table = coef_table,
       metrics = metrics, diagnostics = diagnostics,
       residuals = resid_train, smooth_grid = smooth_grid,
       spline_df = spline_df, hac_lag = hac_lag,
       train_n = length(train_id), test_n = length(test_id),
       train_end = analysis$date_time[n_train],
       test_start = analysis$date_time[n_train + 1L])
}

# Main specification.
main_fit <- fit_profile(spline_df = 8L, hac_lag = 72L)
write.csv(main_fit$coef_table, "profile_coefficients_HAC.csv", row.names = FALSE)
write.csv(main_fit$metrics, "holdout_prediction_metrics.csv", row.names = FALSE)
write.csv(main_fit$diagnostics, "residual_serial_dependence_diagnostics.csv", row.names = FALSE)
write.csv(main_fit$smooth_grid, "estimated_smooth_component.csv", row.names = FALSE)
write.csv(as.data.frame(main_fit$A), "estimated_A_matrix.csv")

# Spline-dimension sensitivity: all estimates are evaluated on the same temporal split.
spline_dfs <- c(6L, 8L, 10L)
spline_sensitivity <- do.call(rbind, lapply(spline_dfs, function(df) {
  obj <- fit_profile(spline_df = df, hac_lag = 72L)
  cbind(data.frame(spline_df = df), obj$metrics)
}))
write.csv(spline_sensitivity, "spline_dimension_sensitivity.csv", row.names = FALSE)

# HAC bandwidth sensitivity for the main spline dimension.
hac_lags <- c(12L, 36L, 72L, 144L)
hac_sensitivity <- do.call(rbind, lapply(hac_lags, function(L) {
  obj <- fit_profile(spline_df = 8L, hac_lag = L)
  cbind(data.frame(hac_lag = L), obj$coef_table)
}))
write.csv(hac_sensitivity, "HAC_bandwidth_sensitivity.csv", row.names = FALSE)

# Correlation of the two transformed outcomes and of the fitted training residuals.
correlations <- data.frame(
  quantity = c("Correlation of log1p outcomes", "Correlation of training residuals"),
  correlation = c(cor(Y[, 1], Y[, 2]), cor(main_fit$residuals[, 1], main_fit$residuals[, 2]))
)
write.csv(correlations, "response_and_residual_correlations.csv", row.names = FALSE)

# Plots saved as PNG files.
png("appliances_energy_time_series_example.png", width = 1200, height = 700, res = 120)
week_n <- min(7L * 24L * 6L, nrow(analysis))
plot(analysis$date_time[seq_len(week_n)], Y[seq_len(week_n), 1], type = "l",
     xlab = "Time", ylab = "log(1 + Appliances energy, Wh)",
     main = "Appliances energy use: first observed week")
dev.off()

png("estimated_smooth_temperature_effects.png", width = 1200, height = 700, res = 120)
plot(main_fit$smooth_grid$T_out, main_fit$smooth_grid$m_log1p_Appliances,
     type = "l", lwd = 2, xlab = "Outdoor temperature T_out (degrees C)",
     ylab = "Estimated smooth component", main = "Estimated smooth temperature components")
lines(main_fit$smooth_grid$T_out, main_fit$smooth_grid$m_log1p_lights, lwd = 2, lty = 2)
legend("topright", legend = c("log1p(Appliances)", "log1p(lights)"),
       lty = c(1, 2), lwd = 2, bty = "n")
dev.off()

png("holdout_actual_vs_predicted.png", width = 1200, height = 700, res = 120)
# Reconstruct the main model's holdout predictions for plotting.
# The prediction is repeated here to keep the fit object compact.
df <- 8L; n_internal <- df - 4L
probs <- seq(0, 1, length.out = n_internal + 2L)[-c(1L, n_internal + 2L)]
knots <- as.numeric(quantile(X_train, probs = probs, names = FALSE, type = 7))
bounds <- range(X)
Bt <- splines::bs(X_train, knots = knots, degree = 3, intercept = TRUE, Boundary.knots = bounds)
Bv <- splines::bs(X_test, knots = knots, degree = 3, intercept = TRUE, Boundary.knots = bounds)
cy <- qr.solve(Bt, Y_train); cz <- qr.solve(Bt, Z_train)
myv <- Bv %*% cy; mzv <- Bv %*% cz
A_main <- main_fit$A
pred_test <- Z_test %*% t(A_main) + myv - mzv %*% t(A_main)
plot(Y_test[, 1], pred_test[, 1], pch = 16, cex = 0.45,
     xlab = "Observed log1p(Appliances)", ylab = "Predicted log1p(Appliances)",
     main = "Chronological holdout: observed versus predicted")
abline(0, 1, lty = 2, lwd = 2)
dev.off()

cat("Analysis completed. Output CSV and PNG files were written to:", getwd(), "\n")
cat("Rows used:", n, " | Training rows:", n_train, " | Holdout rows:", n - n_train, "\n")
cat("Training period ends:", format(main_fit$train_end),
    " | Holdout period starts:", format(main_fit$test_start), "\n")
cat("Main model: cubic B-spline df =", main_fit$spline_df,
    "; HAC lag =", main_fit$hac_lag, "\n")
