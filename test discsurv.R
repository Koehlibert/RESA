library(discSurv)
library(dplyr)
library(eha)
library(Matrix)
data("scania")
scania <- scania %>%
  mutate(exit = ceiling(exit),
         birthdate = floor(birthdate),
         spell = exit - enter) %>% #spell refers to the observed duration of a person
  mutate(enter = enter - 50,
         exit = exit - 50)

scania_long <- scania %>%
  dataLong(timeColumn = "spell", "event", timeAsFactor = FALSE) %>%
  mutate(sex = as.numeric(sex == "female"),
         ses = as.numeric(ses == "upper"),
         immigrant = as.numeric(immigrant == "yes"),)

mod <- glm(y ~ timeInt + sex + ses + immigrant, data = scania_long,
           family = binomial("logit"))

n <- nrow(scania_long)
w_i <- rep(1, n)
w_1 <- 1
W <- Diagonal(n, mod$fitted.values*(1-mod$fitted.values))
X <- Matrix(cbind(rep(1, n), as.matrix(scania_long[,c("timeInt", "sex",
                                            "ses", "immigrant")])))
tmp <- solve(t(X)%*%W%*%X)
Q <- X %*% tmp %*% t(X)
pi <- mod$fitted.values
xi <- 0.5*diag(Q)*((1+w_1)*pi-w_1)
bias <- tmp %*% t(X) %*% W %*% xi
beta_tilde <- beta_hat-bias
