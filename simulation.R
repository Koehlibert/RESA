library(discSurv)
library(Matrix)
set.seed(1)

n_ints <- c(4, 6, 8, 10)
cens_rates <- c(0.25, 0.5, 0.75)
n_n <- c(1000, 2500)
noise_sds <- c(1)
p <- 5
des_mat <- expand.grid(seq_along(n_ints),
                       seq_along(cens_rates),
                       seq_along(n_n),
                       seq_along(noise_sds))
beta <- c(-1, -0.5, 0, 0.5, 1)
time_betas <- readRDS("time_coefs.rds")
n_sim <- 1000
for (i_setting in 1:nrow(des_mat)) {
  n_int <- n_ints[des_mat[i_setting, 1]]
  cens_rate <- cens_rates[des_mat[i_setting, 2]]
  n <- n_n[des_mat[i_setting, 3]]
  noise_sd <- noise_sds[des_mat[i_setting, 4]]
  time_beta_t <- time_betas[[des_mat[i_setting, 1]]][[1]]
  time_probs_c <-
    time_betas[[des_mat[i_setting, 1]]][[2]][[des_mat[i_setting, 2]]]
  res_tmp <- matrix(nrow = 0,
                    ncol = 2 * (n_int + length(beta)))
  actual_cens_rate <- rep(0, n_sim)
  for (i in 1:n_sim) {
    X <- matrix(rnorm(n * p), nrow = n)
    eta <-
      t(apply(matrix(rep(X %*% beta, n_int), nrow = n), 1, function(x)
        x + time_beta_t + rnorm(n_int, 0, noise_sd)))
    probs <- apply(eta, 1:2, function(x)
      1 / (1 + exp(-x)))
    event <- apply(eta, 1:2, function(x)
      as.numeric(runif(1) < x))

    right_censored <- which(rowSums(event) == 0)
    event[right_censored, ncol(event)] <- 1
    ev_time <-
      apply(event, 1, function(x)
        min(which(x != 0)))

    cens_time <- sample(1:(n_int),
                        size = n,
                        prob = time_probs_c[1:(n_int)],
                        replace = TRUE)

    cens_indic <- as.numeric(ev_time <= cens_time)
    actual_cens_rate[i] <- sum(cens_indic == 0) / n

    time <- cbind(c(cens_time), ev_time)[cens_indic * n + 1:n]

    data <- data.frame(time_int = ev_time,
                       status = cens_indic,
                       X)

    data_long <-
      dataLong(
        data,
        timeColumn = "time_int",
        eventColumn = "status",
        timeAsFactor = TRUE
      )

    surv_mod <-
      glm(y ~ timeInt + X1 + X2 + X3 + X4 + X5 + 0,
          data = data_long,
          family = binomial(link = "logit"))

    beta_hat <- coef(surv_mod)[c(paste("timeInt",
                                       1:n_int, sep = ""),
                                 paste("X", 1:length(beta), sep = ""))]

    w_i <- rep(1, nrow(data_long))
    w_1 <- 1
    # W <- diag(n)
    W <-
      Diagonal(nrow(data_long),
               surv_mod$fitted.values * (1 - surv_mod$fitted.values))
    dummy_time <- t(sapply(data_long$timeInt, function(x) {
      tmp <- rep(0, n_int)
      tmp[x] <- 1
      tmp
    }))
    X_bg <-
      Matrix(cbind(dummy_time, as.matrix(data_long[, paste("X",
                                                           1:length(beta),
                                                           sep = "")])))
    tmp <- solve(t(X_bg) %*% W %*% X_bg)
    Q <- X_bg %*% tmp %*% t(X_bg)
    pi <- surv_mod$fitted.values
    xi <- 0.5 * diag(Q) * ((1 + w_1) * pi - w_1)
    bias <- tmp %*% t(X_bg) %*% W %*% xi
    beta_tilde <- beta_hat - bias

    res_tmp <-
      rbind(res_tmp, c(beta_hat, beta_tilde@x))
  }
  saveRDS(actual_cens_rate, file = paste("c_rate", i_setting, ".RDS", sep = ""))
  saveRDS(res_tmp, file = paste("res", i_setting, ".RDS", sep = ""))
}
