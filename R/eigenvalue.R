eigenvalue <- function(x) {
  d <- svd(x, nu = 0, nv = 0)$d[1]
  d <- d / sqrt(nrow(x))
  return(d)
}
