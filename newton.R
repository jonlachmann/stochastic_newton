sigmoid <- function (x) {
  1 / (1 + exp(-x))
}

log_hessian <- function (X, theta) {
  eta <- t(theta) %*% t(X)
  mu <- sigmoid(eta)
  mu_var <- mu * (1 - mu)
  Z_sig <- diag(as.numeric(mu_var)) # This is W in WLS in IRLS
  hess <- t(X) %*% Z_sig %*% X
}

log_grad <- function (X, y, theta) {
  t(X) %*% (sigmoid(X %*% theta) - y)
}

log_weights <- function (X, theta) {
  eta <- X %*% theta
  opexp <- 1 + exp(eta)
  mu_eta <- exp(eta) / opexp^2
  mu <- sigmoid(eta)
  var_mu <- mu * (1 - mu)
  w <- (mu_eta ^ 2) / var_mu
  return(w)
}

newton <- function (X, y, theta, num_iter) {
  # the loop for the iterative process
  for (i in seq_len(num_iter)) {
    grad <- log_grad(X, y, theta)
    hess <- log_hessian(X, theta)

    theta <- theta - ginv(hess) %*% grad
  }
  return(theta)
}

newton2 <- function (X, y, num_iter) {
  theta <- matrix(rep(0, ncol(X)), nrow=ncol(X))

  for (i in seq_len(num_iter)) {
    eta <- X %*% theta # Linear predictor
    mu <- sigmoid(eta) # Fitted values
    mu_var <- mu * (1 - mu) # This is W in WLS in IRLS

    grad <- t(X) %*% (mu - y) # Gradient

    mu_var_mat <- diag(as.numeric(mu_var)) # Variance of fitted values
    hess <- t(X) %*% mu_var_mat %*% X # Hessian

    theta <- theta - ginv(hess) %*% grad # Newton step
  }
  return(list(theta=theta, weights=mu_var))
}

newton3 <- function (X, y, num_iter, ctrl=list(subs=0.1)) {
  theta <- matrix(rep(0, ncol(X)), nrow=ncol(X))
  nobs <- nrow(X)
  w <- matrix(0.43, nobs, 1)
  sub_size <- nobs * ctrl$subs
  if (ctrl$subs != 1) {
    subsi <- sample.int(nobs, sub_size, replace = T)
  } else {
    subsi <- 1 : nobs
  }
  for (i in seq_len(num_iter)) {
    eta <- X[subsi, , drop=F] %*% theta # Linear predictor
    mu <- sigmoid(eta) # Fitted values
    mu_var <- mu * (1 - mu) # This is W in WLS in IRLS

    w[subsi] <- mu_var # Add new weights to weights vector

    grad <- t(X[subsi, , drop=F]) %*% (mu - y[subsi]) # Gradient

    mu_var_mat <- diag(as.numeric(mu_var)) # Variance of fitted values
    hess <- t(X[subsi, , drop=F]) %*% mu_var_mat %*% X[subsi, , drop=F] # Hessian

    theta <- theta - 0.5 * ginv(hess) %*% grad # Newton step

    if (ctrl$subs != 1) subsi <- wrswoR::sample_int_expj(nobs, sub_size, prob = (w + 0.1))
  }
  return(list(theta=theta, weights=mu_var))
}
