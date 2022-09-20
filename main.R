library(MASS)

source("data.R")
source("newton.R")

res <- newton(X, y, matrix(c(0,0,0,0), 4), 3)
res2 <- newton2(X, y, 3)
res3 <- newton3(X, y, 10, ctrl=list(subs=0.1))

glmres <- glm.fit(X, y, family = binomial())

irls.sgd::get_deviance(res2$theta, X, y, binomial())
irls.sgd::get_deviance(res3$theta, X, y, binomial())
irls.sgd::get_deviance(glmres$coefficients, X, y, binomial())

glmres$coefficients
