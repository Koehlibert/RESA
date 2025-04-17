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
                        beta = c(-1, -0.5, 0, 0.5, 1),
                        censrate = 0.1) {
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
  temp$Ttrue <-
    as.vector(exp(-(X %*% beta - err) / lambda))

  temp$id     <- 1:nrow(temp)

  t_admin_cens <- #Inf
    quantile(temp$Ttrue, 0.975)

  # random censoring times, exponentially distributed,
  temp$Ctrue <- rexp(n, rate = censrate)
  temp$status <- (temp$Ttrue <= temp$Ctrue) + 0
  temp$time <- pmin(temp$Ttrue, temp$Ctrue)
  admin_cens <- temp$time > t_admin_cens
  temp$status[admin_cens] <- 0
  temp$time[admin_cens] <- t_admin_cens
  return(temp)
}

temp_c50 <-
  getSurvData(
    seed = 1,
    n = 100000,
    r_scores = 0,
    mu = 0,
    gumbel_loc = 0,
    gumbel_scale = 1,
    lambda = 5,
    censrate = 0.625
  )
temp_c25 <-
  getSurvData(
    seed = 1,
    n = 100000,
    r_scores = 0,
    mu = 0,
    gumbel_loc = 0,
    gumbel_scale = 1,
    lambda = 5,
    censrate = 0.425
  )
temp_c75 <-
  getSurvData(
    seed = 1,
    n = 100000,
    r_scores = 0,
    mu = 0,
    gumbel_loc = 0,
    gumbel_scale = 1,
    lambda = 5,
    censrate = 1.5
  )

# table(temp$status)
# par(mfrow = c(1, 3))
# plot(density(temp$Ttrue))
# plot(density(temp$Ctrue))
# plot(density(temp$time))

get_coefs_from_data <- function(data)
{
  n_ints <- c(4, 6, 8, 10)
  ints_ls <- lapply(n_ints, function(n_int) {
    seq(0, max(data$time), length.out = n_int)
  })
  disc_data_ls <- lapply(ints_ls, function(ints) {
    lapply(c("Ttrue", "Ctrue"), function(ev_col) {
      contToDisc(data,
                 ev_col,
                 ints,
                 timeAsFactor = TRUE) %>%
        mutate(timeDisc = ifelse(is.na(timeDisc), tail(levels(timeDisc), 1), timeDisc))
    })
  })
  time_coefs <- lapply(disc_data_ls, function(data_ls) {
    lapply(data_ls, function(data) {
      data_long <-
        dataLong(data, "timeDisc", "status", timeAsFactor = TRUE)
      T_glm_surv = glm(
        formula = y ~ timeInt + X1 + X2 + X3 + X4 + X5 - 1,
        data = data_long,
        family = "binomial"(link = "logit")
      )
      coef(T_glm_surv)[paste("timeInt", unique(data_long$timeInt), sep = "")]
    })
  })
  time_coefs
}
time_coefs_c25 <- get_coefs_from_data(temp_c25)
saveRDS(time_coefs_c25, "time_coefs_c25.rds")
time_coefs_c50 <- get_coefs_from_data(temp_c50)
saveRDS(time_coefs_c50, "time_coefs_c50.rds")
time_coefs_c75 <- get_coefs_from_data(temp_c75)
saveRDS(time_coefs_c75, "time_coefs_c75.rds")
