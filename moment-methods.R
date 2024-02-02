# Code for implementing the moment-based methods

# plm: pseudolikelihood model
plm <- function(y, X, Z, betastart, alphastart, tol=1e-6, max_iter=100, verbose=FALSE) {
  result_list <- list()
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  get_mu <- function(beta) c(exp(X %*% beta))
  get_phi <- function(alpha) c(exp(Z %*% alpha))
  
  get_grad <- function(alpha, beta) {
    mu <- get_mu(beta)
    phi <- get_phi(alpha)
    error <- y - mu
    vals1 <- (error / phi) * X
    grad1 <- colMeans(vals1)
    
    r2 <- error^2 / mu
    vals2 <- (r2 / phi - 1) * Z
    grad2 <- colMeans(vals2)
    
    grad <- c(grad1, grad2)
    grad
  }
  
  get_hess <- function(alpha, beta) {
    mu <- get_mu(beta)
    phi <- get_phi(alpha)
    error <- y - mu
    r2 <- error^2 / mu
    
    # 11
    weight11 <- sqrt(mu / phi)
    hess11 <- -crossprod(weight11 * X)
      
    # 12
    hess12 <- -(t((error / phi) * X) %*% Z)
      
    # 21
    hess21 <- -t(((y^2 - mu^2) / (phi * mu)) * Z) %*% X
      
    # 22
    weight22 <- sqrt(r2 / phi)
    hess22 <- -crossprod(weight22 * Z)
      
    # Arrange in block matrix
    hess <- rbind(
      cbind(hess11, hess12),
      cbind(hess21, hess22)) / n
    hess
  }
  
  result_list$grad <- get_grad(alphastart, betastart)
  result_list$hess <- get_hess(alphastart, betastart)
  result_list
}

Z <- model.matrix(~x3)
mod_quasipois <- glm(y ~ x1 + x2, data=dat, family=quasipoisson())
phi0 <- summary(mod_quasipois)$dispersion
alpha_start <- rep(0, ncol(Z))
alpha_start[1] <- log(phi0)

mod_plm <- plm(y, X, Z, betastart = beta_start, alphastart = alpha_start)

# delta <- 1e-8
# beta_start2 <- beta_start
# beta_start2[3] <- beta_start2[3] + delta
# alpha_start2 <- alpha_start
# alpha_start2[2] <- alpha_start2[2] + delta
# mod_plm2 <- plm(y, X, Z, betastart = beta_start2, alphastart = alpha_start2)
# (mod_plm2$grad - mod_plm$grad) / delta
# mod_plm$hess[,5]