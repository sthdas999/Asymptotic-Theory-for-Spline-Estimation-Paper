############################################################
# Simulation Study
# Asymptotic Theory for Adaptive Spline Estimation in a
# Multivariate Semiparametric Regression Model
############################################################

rm(list = ls())
set.seed(12345)

############################################################
# Required packages
############################################################

library(splines)
library(MASS)

############################################################
# Simulation settings
############################################################

n.values <- c(200, 500, 1000)

n.rep <- 1000

q <- 2                  # response dimension
d <- 2                  # dimension of linear covariates
p <- 1                  # dimension of nonparametric covariate

alpha <- 0.05           # significance level
conf.level <- 1 - alpha

# Smoothness parameter
# The manuscript does not explicitly specify a numerical value
# of r. It is therefore kept as a user-controlled parameter.
r <- 2

############################################################
# True parameter matrix
############################################################

A.true <- matrix(c(1, 0.5,
                   -0.5, 1),
                 nrow = 2,
                 byrow = TRUE)

A.true
############################################################
# True nonparametric regression function
############################################################

m.true <- function(x) {
  
  cbind(
    sin(2 * pi * x),
    cos(2 * pi * x)
  )
  
}
############################################################
# Choice of spline dimension
############################################################

choose.K <- function(n, r = 2) {
  
  K <- ceiling(n^(1 / (2 * r + 1)))
  
  return(K)
  
}

choose.K(200, r)
choose.K(500, r)
choose.K(1000, r)
############################################################
# Generate error vectors under the five models
############################################################

generate.errors <- function(n, x, model) {
  
  eps <- matrix(0, nrow = n, ncol = 2)
  
  ##########################################################
  # Model 1:
  # Homoscedastic independent Gaussian errors
  #
  # epsilon_i ~ N(0, 0.5 I_2)
  ##########################################################
  
  if (model == 1) {
    
    Sigma <- matrix(c(0.5, 0,
                      0, 0.5),
                    nrow = 2,
                    byrow = TRUE)
    
    eps <- MASS::mvrnorm(
      n = n,
      mu = c(0, 0),
      Sigma = Sigma
    )
  }
  
  ##########################################################
  # Model 2:
  # Heteroscedastic Gaussian errors
  #
  # Sigma(x) =
  # [0.3 + x^2       0       ]
  # [0                0.5+0.5x]
  ##########################################################
  
  if (model == 2) {
    
    for (i in 1:n) {
      
      Sigma.i <- matrix(
        c(
          0.3 + x[i]^2, 0,
          0,            0.5 + 0.5 * x[i]
        ),
        nrow = 2,
        byrow = TRUE
      )
      
      eps[i, ] <- MASS::mvrnorm(
        n = 1,
        mu = c(0, 0),
        Sigma = Sigma.i
      )
    }
  }
  
  ##########################################################
  # Model 3:
  # Correlated multivariate Gaussian errors
  #
  # Sigma =
  # [0.5  0.4]
  # [0.4  1.0]
  ##########################################################
  
  if (model == 3) {
    
    Sigma <- matrix(
      c(
        0.5, 0.4,
        0.4, 1.0
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    eps <- MASS::mvrnorm(
      n = n,
      mu = c(0, 0),
      Sigma = Sigma
    )
  }
  
  ##########################################################
  # Model 4:
  # Non-Gaussian heavy-tailed errors
  #
  # epsilon_i ~ t_5(0, Sigma)
  #
  # For a multivariate t distribution with df = 5,
  # Covariance = (5/(5-2)) Sigma.
  #
  # We therefore choose the scale matrix 0.3 I_2 so that
  # the resulting marginal covariance is approximately 0.5 I_2.
  ##########################################################
  
  if (model == 4) {
    
    df <- 5
    
    Sigma.scale <- matrix(
      c(
        0.3, 0,
        0,   0.3
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    # Generate multivariate t observations:
    #
    # epsilon = Z / sqrt(U/df),
    #
    # where Z ~ N(0, Sigma.scale)
    # and U ~ Chi-square(df).
    
    Z <- MASS::mvrnorm(
      n = n,
      mu = c(0, 0),
      Sigma = Sigma.scale
    )
    
    U <- rchisq(n, df = df)
    
    eps <- Z / sqrt(U / df)
  }
  
  ##########################################################
  # Model 5:
  # Nonlinear covariate-dependent covariance
  #
  # Sigma(x) =
  #
  # [0.4 + 0.5 sin^2(2 pi x)   0.3 cos(2 pi x)]
  #
  # [0.3 cos(2 pi x)            0.6 + x       ]
  ##########################################################
  
  if (model == 5) {
    
    for (i in 1:n) {
      
      Sigma.i <- matrix(
        c(
          0.4 + 0.5 * sin(2 * pi * x[i])^2,
          0.3 * cos(2 * pi * x[i]),
          
          0.3 * cos(2 * pi * x[i]),
          0.6 + x[i]
        ),
        nrow = 2,
        byrow = TRUE
      )
      
      eps[i, ] <- MASS::mvrnorm(
        n = 1,
        mu = c(0, 0),
        Sigma = Sigma.i
      )
    }
  }
  
  return(eps)
}
############################################################
# Generate one complete data set
############################################################

generate.data <- function(n, model, r = 2) {
  
  ##########################################################
  # Generate X and Z
  ##########################################################
  
  X <- runif(n, 0, 1)
  
  Z <- matrix(
    rnorm(n * d),
    nrow = n,
    ncol = d
  )
  
  ##########################################################
  # True nonparametric component
  ##########################################################
  
  M <- m.true(X)
  
  ##########################################################
  # Generate model-specific errors
  ##########################################################
  
  eps <- generate.errors(
    n = n,
    x = X,
    model = model
  )
  
  ##########################################################
  # Generate response
  #
  # Y_i = A Z_i + m(X_i) + epsilon_i
  ##########################################################
  
  Y <- Z %*% t(A.true) + M + eps
  
  return(
    list(
      X = X,
      Z = Z,
      Y = Y,
      M = M,
      eps = eps
    )
  )
}
############################################################
# Construct cubic B-spline basis
############################################################

construct.spline.basis <- function(x, K) {
  
  n <- length(x)
  
  # Number of interior knots
  n.inter <- max(1, K - 4)
  
  knots <- seq(
    from = 0,
    to = 1,
    length.out = n.inter + 2
  )
  
  knots <- knots[-c(1, length(knots))]
  
  B <- splines::bs(
    x,
    degree = 3,
    knots = knots,
    Boundary.knots = c(0, 1),
    intercept = TRUE
  )
  
  return(B)
}
############################################################
# Profile spline estimator
############################################################

profile.spline.estimator <- function(Y, Z, X, r = 2) {
  
  n <- nrow(Y)
  
  ##########################################################
  # Select spline dimension
  ##########################################################
  
  K <- choose.K(n, r)
  
  ##########################################################
  # Construct cubic B-spline basis
  ##########################################################
  
  B <- construct.spline.basis(X, K)
  
  ##########################################################
  # Estimate conditional mean of Y given X
  #
  # Theta_Y_hat =
  # (B'B)^(-1) B'Y
  ##########################################################
  
  BtB <- crossprod(B)
  
  BtB.inv <- solve(
    BtB + 1e-10 * diag(ncol(B))
  )
  
  Theta.Y <- BtB.inv %*% crossprod(B, Y)
  
  ##########################################################
  # Estimate conditional mean of Z given X
  ##########################################################
  
  Theta.Z <- BtB.inv %*% crossprod(B, Z)
  
  ##########################################################
  # Fitted conditional means
  ##########################################################
  
  mu.Y.hat <- B %*% Theta.Y
  
  mu.Z.hat <- B %*% Theta.Z
  
  ##########################################################
  # Profile residuals
  ##########################################################
  
  Y.tilde <- Y - mu.Y.hat
  
  Z.tilde <- Z - mu.Z.hat
  
  ##########################################################
  # Profile least-squares estimator
  ##########################################################
  
  S.zz <- crossprod(Z.tilde)
  
  S.yz <- crossprod(Y.tilde, Z.tilde)
  
  A.hat <- S.yz %*%
    solve(
      S.zz + 1e-10 * diag(d)
    )
  
  ##########################################################
  # Estimated nonparametric function
  ##########################################################
  
  m.hat <- function(x.new) {
    
    B.new <- construct.spline.basis(
      x.new,
      K
    )
    
    B.new %*% Theta.Y -
      (B.new %*% Theta.Z) %*% t(A.hat)
    
  }
  
  ##########################################################
  # Residuals from fitted model
  ##########################################################
  
  eps.hat <- Y -
    Z %*% t(A.hat) -
    mu.Y.hat +
    mu.Z.hat %*% t(A.hat)
  
  ##########################################################
  # Return results
  ##########################################################
  
  list(
    A.hat = A.hat,
    B = B,
    K = K,
    Theta.Y = Theta.Y,
    Theta.Z = Theta.Z,
    mu.Y.hat = mu.Y.hat,
    mu.Z.hat = mu.Z.hat,
    Y.tilde = Y.tilde,
    Z.tilde = Z.tilde,
    eps.hat = eps.hat,
    m.hat = m.hat
  )
}
############################################################
# Sandwich covariance estimator for vec(A_hat)
############################################################

sandwich.se <- function(A.hat,
                        Y.tilde,
                        Z.tilde,
                        eps.hat) {
  
  n <- nrow(Y.tilde)
  
  d <- ncol(Z.tilde)
  q <- ncol(Y.tilde)
  
  ##########################################################
  # Bread
  ##########################################################
  
  S.zz <- crossprod(Z.tilde) / n
  
  bread <- solve(
    S.zz + 1e-10 * diag(d)
  )
  
  ##########################################################
  # Meat
  ##########################################################
  
  meat <- matrix(
    0,
    nrow = q * d,
    ncol = q * d
  )
  
  for (i in 1:n) {
    
    zi <- matrix(
      Z.tilde[i, ],
      nrow = 1
    )
    
    ei <- matrix(
      eps.hat[i, ],
      nrow = 1
    )
    
    # vec(e_i z_i')
    score.i <- as.vector(
      t(ei) %*% zi
    )
    
    meat <- meat +
      tcrossprod(score.i)
  }
  
  meat <- meat / n
  
  ##########################################################
  # Kronecker-product bread
  ##########################################################
  
  bread.full <- kronecker(
    bread,
    diag(q)
  )
  
  ##########################################################
  # Sandwich covariance
  ##########################################################
  
  Sigma.A <- bread.full %*%
    meat %*%
    t(bread.full) / n
  
  ##########################################################
  # Standard errors
  ##########################################################
  
  SE <- sqrt(
    pmax(diag(Sigma.A), 0)
  )
  
  return(
    list(
      covariance = Sigma.A,
      SE = SE
    )
  )
}
############################################################
# Performance measures
############################################################

calculate.performance <- function(
    A.estimates,
    SE.estimates,
    A.true,
    alpha = 0.05) {
  
  R <- nrow(A.estimates)
  
  true.vec <- as.vector(t(A.true))
  
  bias <- colMeans(
    A.estimates
  ) - true.vec
  
  rmse <- sqrt(
    colMeans(
      (sweep(
        A.estimates,
        2,
        true.vec,
        "-"
      ))^2
    )
  )
  
  ##########################################################
  # Confidence interval coverage
  ##########################################################
  
  coverage <- numeric(
    length(true.vec)
  )
  
  z.alpha <- qnorm(
    1 - alpha / 2
  )
  
  for (j in 1:length(true.vec)) {
    
    lower <- A.estimates[, j] -
      z.alpha * SE.estimates[, j]
    
    upper <- A.estimates[, j] +
      z.alpha * SE.estimates[, j]
    
    coverage[j] <- mean(
      lower <= true.vec[j] &
        upper >= true.vec[j]
    )
  }
  
  ##########################################################
  # Average standard error
  ##########################################################
  
  avg.SE <- colMeans(
    SE.estimates,
    na.rm = TRUE
  )
  
  return(
    data.frame(
      Parameter = c(
        "A11",
        "A12",
        "A21",
        "A22"
      ),
      Bias = bias,
      RMSE = rmse,
      Coverage = coverage,
      Avg.SE = avg.SE
    )
  )
}
############################################################
# Integrated squared error
############################################################

calculate.ISE <- function(
    fit,
    grid.size = 1001) {
  
  x.grid <- seq(
    0,
    1,
    length.out = grid.size
  )
  
  ##########################################################
  # Estimated function
  ##########################################################
  
  m.hat.grid <- fit$m.hat(
    x.grid
  )
  
  ##########################################################
  # True function
  ##########################################################
  
  m.true.grid <- m.true(
    x.grid
  )
  
  ##########################################################
  # Squared Euclidean error
  ##########################################################
  
  sq.error <- rowSums(
    (m.hat.grid - m.true.grid)^2
  )
  
  ##########################################################
  # Numerical integration using trapezoidal rule
  ##########################################################
  
  ISE <- sum(
    diff(x.grid) *
      (
        sq.error[-length(sq.error)] +
          sq.error[-1]
      ) / 2
  )
  
  return(ISE)
}
############################################################
# One Monte Carlo replication
############################################################

one.replication <- function(
    n,
    model,
    r = 2) {
  
  ##########################################################
  # Generate data
  ##########################################################
  
  dat <- generate.data(
    n = n,
    model = model,
    r = r
  )
  
  ##########################################################
  # Estimate model
  ##########################################################
  
  fit <- profile.spline.estimator(
    Y = dat$Y,
    Z = dat$Z,
    X = dat$X,
    r = r
  )
  
  ##########################################################
  # Robust standard errors
  ##########################################################
  
  sandwich <- sandwich.se(
    A.hat = fit$A.hat,
    Y.tilde = fit$Y.tilde,
    Z.tilde = fit$Z.tilde,
    eps.hat = fit$eps.hat
  )
  
  ##########################################################
  # Vectorize A_hat
  ##########################################################
  
  A.vec <- as.vector(
    t(fit$A.hat)
  )
  
  ##########################################################
  # ISE
  ##########################################################
  
  ISE <- calculate.ISE(
    fit
  )
  
  ##########################################################
  # Return
  ##########################################################
  
  list(
    A.hat = A.vec,
    SE = sandwich$SE,
    ISE = ISE
  )
}
############################################################
# Main Monte Carlo simulation
############################################################

simulation.results <- list()

start.time <- Sys.time()

for (model in 1:5) {
  
  cat("\n============================================\n")
  cat("Running Model:", model, "\n")
  cat("============================================\n")
  
  model.results <- list()
  
  for (n in n.values) {
    
    cat("\nSample size n =", n, "\n")
    
    ########################################################
    # Storage
    ########################################################
    
    A.storage <- matrix(
      NA,
      nrow = n.rep,
      ncol = 4
    )
    
    SE.storage <- matrix(
      NA,
      nrow = n.rep,
      ncol = 4
    )
    
    ISE.storage <- numeric(
      n.rep
    )
    
    ########################################################
    # Monte Carlo replications
    ########################################################
    
    for (rep in 1:n.rep) {
      
      if (rep %% 100 == 0) {
        
        cat(
          "Replication:",
          rep,
          "of",
          n.rep,
          "\n"
        )
      }
      
      result <- one.replication(
        n = n,
        model = model,
        r = r
      )
      
      A.storage[rep, ] <- result$A.hat
      
      SE.storage[rep, ] <- result$SE
      
      ISE.storage[rep] <- result$ISE
    }
    
    ########################################################
    # Performance table
    ########################################################
    
    performance <- calculate.performance(
      A.estimates = A.storage,
      SE.estimates = SE.storage,
      A.true = A.true,
      alpha = alpha
    )
    
    performance$n <- n
    
    performance$Model <- model
    
    ########################################################
    # ISE summary
    ########################################################
    
    ISE.summary <- data.frame(
      Model = model,
      n = n,
      Mean.ISE = mean(
        ISE.storage,
        na.rm = TRUE
      ),
      SD.ISE = sd(
        ISE.storage,
        na.rm = TRUE
      )
    )
    
    ########################################################
    # Store results
    ########################################################
    
    model.results[[paste0("n", n)]] <- list(
      performance = performance,
      ISE = ISE.summary,
      A.estimates = A.storage,
      SE.estimates = SE.storage,
      ISE.values = ISE.storage
    )
  }
  
  simulation.results[[paste0("Model", model)]] <-
    model.results
}

end.time <- Sys.time()

cat("\nSimulation completed.\n")
cat("Starting time:", start.time, "\n")
cat("Ending time:", end.time, "\n")
############################################################
# Combine parameter results
############################################################

all.performance <- do.call(
  rbind,
  lapply(
    simulation.results,
    function(model.list) {
      
      do.call(
        rbind,
        lapply(
          model.list,
          function(x) x$performance
        )
      )
      
    }
  )
)

rownames(all.performance) <- NULL

all.performance
############################################################
# Combine ISE results
############################################################

all.ISE <- do.call(
  rbind,
  lapply(
    simulation.results,
    function(model.list) {
      
      do.call(
        rbind,
        lapply(
          model.list,
          function(x) x$ISE
        )
      )
      
    }
  )
)

rownames(all.ISE) <- NULL

all.ISE
############################################################
# Table 1 -- Model 1
############################################################

table1 <- subset(
  all.performance,
  Model == 1
)

table1 <- table1[
  order(table1$n),
  c(
    "n",
    "Parameter",
    "Bias",
    "RMSE",
    "Coverage",
    "Avg.SE"
  )
]

table1
############################################################
# Table 2 -- Model 2
############################################################

table2 <- subset(
  all.performance,
  Model == 2
)

table2 <- table2[
  order(table2$n),
  c(
    "n",
    "Parameter",
    "Bias",
    "RMSE",
    "Coverage",
    "Avg.SE"
  )
]

table2
############################################################
# Table 3 -- Model 3
############################################################

table3 <- subset(
  all.performance,
  Model == 3
)

table3 <- table3[
  order(table3$n),
  c(
    "n",
    "Parameter",
    "Bias",
    "RMSE",
    "Coverage",
    "Avg.SE"
  )
]

table3
############################################################
# Table 4 -- Model 4
############################################################

table4 <- subset(
  all.performance,
  Model == 4
)

table4 <- table4[
  order(table4$n),
  c(
    "n",
    "Parameter",
    "Bias",
    "RMSE",
    "Coverage",
    "Avg.SE"
  )
]

table4
############################################################
# Table 5 -- Model 5
############################################################

table5 <- subset(
  all.performance,
  Model == 5
)

table5 <- table5[
  order(table5$n),
  c(
    "n",
    "Parameter",
    "Bias",
    "RMSE",
    "Coverage",
    "Avg.SE"
  )
]

table5
############################################################
# ISE by model and sample size
############################################################

ISE.table <- all.ISE[
  order(
    all.ISE$Model,
    all.ISE$n
  ),
]

ISE.table
############################################################
# Round numerical results
############################################################

all.performance.rounded <- all.performance

all.performance.rounded$Bias <-
  round(
    all.performance.rounded$Bias,
    3
  )

all.performance.rounded$RMSE <-
  round(
    all.performance.rounded$RMSE,
    3
  )

all.performance.rounded$Coverage <-
  round(
    all.performance.rounded$Coverage,
    3
  )

all.performance.rounded$Avg.SE <-
  round(
    all.performance.rounded$Avg.SE,
    3
  )

all.performance.rounded
############################################################
# Save results
############################################################

write.csv(
  all.performance.rounded,
  file = "simulation_parameter_results.csv",
  row.names = FALSE
)

write.csv(
  ISE.table,
  file = "simulation_ISE_results.csv",
  row.names = FALSE
)

save(
  simulation.results,
  all.performance,
  all.ISE,
  file = "simulation_results.RData"
)
############################################################
# RMSE versus sample size
############################################################

library(graphics)

rmse.A11 <- subset(
  all.performance,
  Parameter == "A11"
)

plot(
  rmse.A11$n,
  rmse.A11$RMSE,
  type = "b",
  log = "xy",
  xlab = "Sample size n",
  ylab = "RMSE",
  main = "RMSE of A11"
)
############################################################
# ISE versus sample size
############################################################

plot(
  ISE.table$n[ISE.table$Model == 1],
  ISE.table$Mean.ISE[ISE.table$Model == 1],
  type = "b",
  log = "xy",
  xlab = "Sample size n",
  ylab = "Mean ISE",
  main = "ISE under Model 1"
)