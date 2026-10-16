############################################################
# Real Data Analysis: Appliances Energy Prediction Dataset
# Subsections:
# 1. Descriptive Analysis
# 2. Estimated Parametric Component
############################################################

# Install packages if required
# install.packages("splines")  # Usually included with base R

library(splines)

# ----------------------------------------------------------
# 1. Import and prepare the data
# ----------------------------------------------------------

# Set the working directory to the folder containing the CSV
setwd("C:/Users/hp/Desktop/appliances+energy+prediction")

data <- read.csv("energydata_complete.csv",
                 stringsAsFactors = FALSE)

# Convert date and time
data$date <- as.POSIXct(data$date,
                        format = "%Y-%m-%d %H:%M:%S",
                        tz = "UTC")

# Remove observations with missing values in the variables used
# in the analysis
variables <- c(
  "date", "Appliances", "lights",
  "T_out", "RH_out", "Windspeed",
  "Press_mm_hg", "Visibility", "Tdewpoint"
)

data <- data[complete.cases(data[, variables]), ]

# Extract hour of day
data$hour <- as.numeric(format(data$date, "%H"))

# Response matrix: appliance and lighting energy consumption
Y <- as.matrix(data[, c("Appliances", "lights")])

# Parametric covariates
Z <- as.matrix(data[, c(
  "T_out", "RH_out", "Windspeed",
  "Press_mm_hg", "Visibility", "Tdewpoint"
)])

# Centre the parametric covariates
Z_mean <- colMeans(Z)
Zc <- sweep(Z, 2, Z_mean, FUN = "-")

n <- nrow(data)

# ----------------------------------------------------------
# 2. Descriptive Analysis
# ----------------------------------------------------------

cat("\nNumber of observations:", n, "\n")
cat("Study period:\n")
print(range(data$date))

# Descriptive statistics for appliance energy consumption
cat("\nDescriptive statistics: Appliances\n")
print(c(
  Mean = mean(data$Appliances),
  SD = sd(data$Appliances),
  Median = median(data$Appliances),
  Q1 = quantile(data$Appliances, 0.25),
  Q3 = quantile(data$Appliances, 0.75),
  P90 = quantile(data$Appliances, 0.90),
  P95 = quantile(data$Appliances, 0.95),
  P99 = quantile(data$Appliances, 0.99)
))

# Descriptive statistics for lighting energy consumption
cat("\nDescriptive statistics: Lights\n")
print(c(
  Mean = mean(data$lights),
  SD = sd(data$lights),
  Median = median(data$lights),
  Q1 = quantile(data$lights, 0.25),
  Q3 = quantile(data$lights, 0.75),
  P90 = quantile(data$lights, 0.90),
  P95 = quantile(data$lights, 0.95),
  P99 = quantile(data$lights, 0.99)
))

# Percentage of zero lighting observations
zero_lights <- 100 * mean(data$lights == 0)

cat("\nPercentage of zero lighting observations:",
    round(zero_lights, 2), "%\n")

# Correlation between the two response variables
cat("\nCorrelation between Appliances and Lights:\n")
print(cor(data$Appliances, data$lights))

# Hourly mean energy consumption
hourly_means <- aggregate(
  cbind(Appliances, lights) ~ hour,
  data = data,
  FUN = mean
)

cat("\nHourly mean energy consumption:\n")
print(hourly_means, row.names = FALSE)

# Hours corresponding to the minimum and maximum hourly means
min_appliance_hour <- hourly_means$hour[
  which.min(hourly_means$Appliances)
]

max_appliance_hour <- hourly_means$hour[
  which.max(hourly_means$Appliances)
]

min_lights_hour <- hourly_means$hour[
  which.min(hourly_means$lights)
]

max_lights_hour <- hourly_means$hour[
  which.max(hourly_means$lights)
]

cat("\nMinimum mean appliance consumption occurs at hour:",
    min_appliance_hour, "\n")
cat("Maximum mean appliance consumption occurs at hour:",
    max_appliance_hour, "\n")
cat("Minimum mean lighting consumption occurs at hour:",
    min_lights_hour, "\n")
cat("Maximum mean lighting consumption occurs at hour:",
    max_lights_hour, "\n")


# ----------------------------------------------------------
# 3. Estimated Parametric and Nonparametric Components
# ----------------------------------------------------------

# Construct a cubic B-spline basis with fewer basis functions.
# There are only 24 distinct hourly values (0,...,23).
B <- splines::bs(
  data$hour,
  df = 12,
  degree = 3,
  intercept = TRUE,
  Boundary.knots = c(0, 23)
)

# Compute the QR decomposition of the spline basis.
qr_B <- qr(B, tol = 1e-10)

cat("Number of spline basis columns:", ncol(B), "\n")
cat("Numerical rank of spline basis:", qr_B$rank, "\n")

if (qr_B$rank < ncol(B)) {
  stop("The spline basis is still rank-deficient. Reduce df.")
}

# Residualize the response and parametric covariates
# with respect to the spline basis.
Y_res <- qr.resid(qr_B, Y)
Z_res <- qr.resid(qr_B, Zc)

# Estimate the parametric coefficient matrix.
Gamma_hat <- solve(
  crossprod(Z_res),
  crossprod(Z_res, Y_res)
)

A_hat <- t(Gamma_hat)

rownames(A_hat) <- c("Appliances", "Lights")
colnames(A_hat) <- colnames(Z)

cat("\nEstimated parametric coefficient matrix A_hat:\n")
print(round(A_hat, 3))

# ----------------------------------------------------------
# 4. Estimate the Nonparametric Component
# ----------------------------------------------------------

# Remove the estimated parametric contribution.
Y_adjusted <- Y - Zc %*% Gamma_hat

# Estimate spline coefficients by least squares.
spline_coef <- qr.coef(qr_B, Y_adjusted)

# Check for missing coefficient estimates.
if (anyNA(spline_coef)) {
  stop("NA coefficients remain. Check the spline basis rank.")
}

rownames(spline_coef) <- paste0("Basis_", seq_len(nrow(spline_coef)))
colnames(spline_coef) <- c("Appliances", "Lights")

cat("\nEstimated spline coefficients:\n")
print(round(spline_coef, 6))

# Estimated smooth component
m_hat <- B %*% spline_coef
colnames(m_hat) <- c("m_Appliances", "m_Lights")

# Verify the dimensions and numerical validity.
cat("\nDimensions of m_hat:\n")
print(dim(m_hat))

cat("\nAny missing values in m_hat?\n")
print(anyNA(m_hat))

# Fitted responses
Y_hat <- Zc %*% Gamma_hat + m_hat

# Residual sum of squares
RSS <- colSums((Y - Y_hat)^2)

cat("\nResidual sum of squares:\n")
print(RSS)


# ----------------------------------------------------------
# 5. In-sample fit summaries
# ----------------------------------------------------------
residuals_hat<- Y-Y_hat
rmse <- sqrt(colMeans(residuals_hat^2))
mae <- colMeans(abs(residuals_hat))

SSE <- colSums(residuals_hat^2)
SST <- colSums(sweep(Y, 2, colMeans(Y), FUN = "-")^2)
R_squared <- 1 - SSE / SST

fit_summary <- data.frame(
  Response = colnames(Y),
  RMSE = as.numeric(rmse),
  MAE = as.numeric(mae),
  R_squared = as.numeric(R_squared)
)

cat("\nIn-sample fit summary:\n")
print(fit_summary, row.names = FALSE)

cat("\nCorrelation between fitted-model residuals:\n")
print(cor(residuals_hat[, 1], residuals_hat[, 2]))

# Optional: save the principal numerical results
write.csv(hourly_means, "hourly_mean_energy.csv",
          row.names = FALSE)

write.csv(A_hat, "estimated_parametric_coefficients.csv")

write.csv(fit_summary, "model_fit_summary.csv",
          row.names = FALSE)

############################################################
# End of code
############################################################
