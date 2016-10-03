sample.correlated.data <- function(p, n, rho, df, sigma) {
  s <- matrix(rho, p, p)
  s <- s + diag(rep(1 - rho, p))
  cov.mat <- riwish(df, s)
  cor.mat <- cov2cor(cov.mat)
  eig <- eigen(cor.mat)
  L <- eig[[2]] %*% diag(sqrt(eig[[1]])) %*% t(eig[[2]])
  x <- matrix(rnorm(p * n, 0, sigma), p, n)
  xx <- L %*% x
  return(xx)
}
