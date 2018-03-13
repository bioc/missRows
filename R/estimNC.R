estimNC <- function (X, minNC=0, maxNC) {

    crit <- NULL
    p <- ncol(X)
    n <- nrow(X)
    maxNC <- min(c(n - 2, p - 1, maxNC))
    X <- scale(X)

    if (minNC == 0) {
        crit <- mean(X^2, na.rm=TRUE) * (n * p)/(p * (n - 1))
    }

    rr <- svd(X, nu=maxNC, nv=maxNC)

    for (q in max(minNC, 1):maxNC) {

        if (q > 1) {
            rec <- tcrossprod(sweep(rr$u[, 1:q, drop=FALSE], 2, rr$d[1:q], 
                            FUN = "*"), rr$v[, 1:q, drop=FALSE])
        }

        if (q == 1) {
            rec <- tcrossprod(rr$u[, 1, drop=FALSE] * rr$d[1], 
                            rr$v[, 1, drop=FALSE])
        }

        crit <- c(crit, mean((n * p * (X - rec)/((n - 1) * 
                                                p - q * (n + p - q - 1)))^2,
                    na.rm=TRUE))
    }

    if (any(diff(crit) > 0)) {
    ncp <- which(diff(crit) > 0)[1]
    } else { ncp <- which.min(crit) }

    return(ncp + minNC - 1)
}
