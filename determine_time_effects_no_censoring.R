library(readxl)
library(tidyverse)
library(janitor)
library(MASS)
library(survival)
library(VGAM)
library(pec)
library(discSurv)
library(rngtools)
###############################################################################
getSurvData <- function(seed,
                        n,
                        r_scores = 0.5,
                        mu = 0,
                        gumbel_loc = 0,
                        gumbel_scale = 1,
                        lambda = 1,
                        beta = c(-1, -0.5, 0, 0.5, 1)) {
  set.seed(seed)
  d <- length(beta)
  sigma_mat <- matrix(rep(r_scores, d ^ 2), nrow = d)
  diag(sigma_mat) <- 1
  X <-
    mu + mvrnorm(n = n,
                 mu = rep(0, d),
                 Sigma = sigma_mat)
  err <- rgumbel(n, gumbel_loc, gumbel_scale)
  temp <- data.frame(model.matrix(~ X + 0)) # data frame
  temp$time <-
    as.vector(exp(-(X %*% beta - err) / lambda))

  temp$id     <- 1:nrow(temp)

  t_admin_cens <- #Inf
    quantile(temp$time, 0.975)

  # random censoring times, exponentially distributed,
  temp$status <- 1
  admin_cens <- temp$time > t_admin_cens
  temp$status[admin_cens] <- 0
  temp$time[admin_cens] <- t_admin_cens
  return(temp)
}

temp <-
  getSurvData(
    seed = 1,
    n = 100000,
    r_scores = 0,
    mu = 0,
    gumbel_loc = 0,
    gumbel_scale = 1,
    lambda = 5
  )

table(temp$status)
plot(density(temp$time))

get_coefs_from_data <- function(data)
{
  n_ints <- c(4, 6, 8, 10)
  ints_ls <- lapply(n_ints, function(n_int) {
    seq(0, max(data$time), length.out = n_int + 1)[-1]
  })
  disc_data_ls <- lapply(ints_ls, function(ints) {
    contToDisc(data,
               "time",
               ints,
               timeAsFactor = TRUE) %>%
      mutate(timeDisc = ifelse(is.na(timeDisc), tail(levels(timeDisc), 1), timeDisc))
  })
  time_coefs <- lapply(disc_data_ls, function(data_disc) {
    n_int <- max(as.numeric(data_disc$timeDisc))
    data_long <-
      dataLong(data_disc, "timeDisc", "status", timeAsFactor = TRUE)
    T_glm_surv = glm(
      formula = y ~ timeInt + X1 + X2 + X3 + X4 + X5 - 1,
      data = data_long,
      family = "binomial"(link = "logit")
    )
    coef_time <-
      coef(T_glm_surv)[paste("timeInt", unique(data_long$timeInt), sep = "")]
    cens_rates <- c(0.25, 0.5, 0.75)
    rel_ev_table <-
      sapply(2:n_int, function(t) {
        sum(data_disc$timeDisc == t) / nrow(data_disc)
      })
    c_eff_ls <- lapply(cens_rates, function(x) {
      rep((x) / (n_int - 1), (n_int - 1))
    })
    probs_cens <- lapply(c_eff_ls, function(c_eff) {
      a <- rep(1, n_int)
      for (t in (n_int - 1):1) {
        a[t] <- c_eff[t] / rel_ev_table[t]
        if (a[t] > a[t + 1]) {
          tmp <- a[t] - a[t + 1] + 0.001
          a[t] <- a[t] - tmp
          c_eff[t] <- c_eff[t] - tmp * rel_ev_table[t]
          c_eff[(t - 1):1] <- c_eff[(t - 1):1] +
            tmp * (rel_ev_table[t] / (t - 1))
        }
      }
      # a[1:(n_int - 1)] <- c_eff / rel_ev_table
      # for (t in (n_int - 1) : 1) {
      #   if (a[t] > a[t + 1]) {
      #     tmp <- a[t] - a[t + 1] + 0.025
      #     a[t] <- a[t] - tmp
      #     a[1:(t - 1)] <- a[1:(t - 1)] + tmp / (t - 1)
      #   }
      # }
      for (t in 2:(n_int - 1)) {
        a[t] <- a[t] - sum(a[1:(t - 1)])
      }
      a[n_int] <- 1 - sum(a[-n_int])
      a
    })
    time_coefs =
      list(t = coef_time,
           c = probs_cens)
  })
  time_coefs
}
time_coefs <- get_coefs_from_data(temp)
saveRDS(time_coefs, "time_coefs.rds")
