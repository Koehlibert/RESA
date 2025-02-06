library(ggplot2)
library(Matrix)
library(dplyr)

n1 <- 9750
n2 <- 250
n <- n1 + n2

beta_0 <- -10
beta_1 <- -3
beta_2 <- -1
beta <- c(beta_1, beta_2)

b <- 1000
res_df <- matrix(nrow = b, ncol = 9)
for(i in 1:b)
{
  print(i)
  x1 <- c(rnorm(n1), rnorm(n2, 0, 5))
  x2 <- c(rnorm(n1), rnorm(n2, 0, 7))
  p <- 1/(1+exp(beta_0+x1*beta_1+x2*beta_2+rnorm(1)))
  y <- sapply(p, function(x) sample(c(0,1), 1, prob = c(x,1-x)))
  # plot(density(p))
  # table(y)

  mod <- glm(y~x1+x2, family = binomial(link="logit"))
  # summary(mod)
  beta_hat <- coef(mod)[1:3]
  res_df[i, 1:3] <- beta_hat

  w_i <- rep(1, n)
  w_1 <- 1
  # W <- diag(n)
  W <- Diagonal(n, mod$fitted.values*(1-mod$fitted.values))
  X <- Matrix(cbind(rep(1, n), x1, x2))
  tmp <- solve(t(X)%*%W%*%X)
  Q <- X %*% tmp %*% t(X)
  pi <- mod$fitted.values
  xi <- 0.5*diag(Q)*((1+w_1)*pi-w_1)
  bias <- tmp %*% t(X) %*% W %*% xi
  beta_tilde <- beta_hat-bias
  res_df[i, 4:6] <- as.vector(beta_tilde)
  res_df[i, 7:9] <- as.vector(bias)
}
plot_data <-
  data.frame(sqe = c((res_df[,2]-beta[1])^2,
                     (res_df[,3]-beta[2])^2,
                     (res_df[,5]-beta[1])^2,
                     (res_df[,6]-beta[2]))^2,
             coef = rep(rep(1:2, each = b), 2),
             type = rep(c("Pre", "Post"), each = 2*b))
ggplot(plot_data, aes(x = sqe, color = type)) +
  geom_density() + facet_grid(~coef)
plot(density(res_df[,9]))
plot_data_2 <-
  data.frame(e = c((res_df[,2]-beta[1]),
                   (res_df[,3]-beta[2]),
                   (res_df[,5]-beta[1]),
                   (res_df[,6]-beta[2])),
             coef = rep(rep(1:2, each = b), 2),
             type = rep(c("Pre", "Post"), each = 2*b))
ggplot(plot_data_2, aes(x = e, color = type)) +
  geom_density() + facet_grid(~coef)
mean(plot_data_2 %>% filter(type == "Pre") %>% select(e) %>% unlist())
mean(plot_data_2 %>% filter(type == "Post") %>% select(e) %>% unlist())
