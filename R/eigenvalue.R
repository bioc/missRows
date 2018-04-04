eigenvalue <- function(X) {
    d <- svd(X, nu=0, nv=0)$d[1]
    d <- d / sqrt(nrow(X))
    return(d)
}
