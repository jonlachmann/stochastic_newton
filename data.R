nobs <- 300

set.seed(2016)
#simulate data
#independent variables
x1 <- rnorm(nobs, 3, 2) + 0.1 * (1:nobs)
x2 <- rbinom(nobs, 1, 0.3)
x3 <- rpois(n = nobs, lambda = 4)
x3[16:nobs] <- x3[16:nobs] - rpois(n = 15, lambda = 2)

#dependent variable
y <- c(rbinom(nobs/6, 1, 0.1), rbinom(nobs/3, 1, 0.25), rbinom(nobs/3, 1, 0.75), rbinom(nobs/6, 1, 0.9))

x0 <- rep(1, nobs) #bias
X <- cbind(x0, x1, x2, x3)