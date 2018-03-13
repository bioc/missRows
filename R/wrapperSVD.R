wrapperSVD <- function (X, rWeights=NULL, cWeights=NULL, ncp=Inf) {
    
    tryCatchWE <- function(expr) {
        W <- NULL
        
        wHandler <- function(w) {
            W <- w
            invokeRestart("muffleWarning")
        }
        
        list(value = withCallingHandlers(tryCatch(expr, error=function(e) e),
                                        warning=wHandler), warning=W)
    }
    
    if (is.null(rWeights))
        rWeights <- rep(1/nrow(X), nrow(X))
    if (is.null(cWeights))
        cWeights <- rep(1, ncol(X))
    ncp <- min(ncp, nrow(X) - 1, ncol(X))
    rWeights <- rWeights/sum(rWeights)
    X <- t(t(X) * sqrt(cWeights)) * sqrt(rWeights)
    if (ncol(X) < nrow(X)) {
        usuelleSVD <- tryCatchWE(svd(X, nu=ncp, nv=ncp))$val
        if (names(usuelleSVD)[[1]] == "message") {
            usuelleSVD <- tryCatchWE(svd(t(X), nu=ncp, nv=ncp))$val
            if (names(usuelleSVD)[[1]] == "d") {
                aux <- usuelleSVD$u
                usuelleSVD$u <- usuelleSVD$v
                usuelleSVD$v <- aux
            }
            else {
                bb <- eigen(crossprod(X, X), symmetric=TRUE)
                usuelleSVD <- vector(mode="list", length=3)
                usuelleSVD$d[usuelleSVD$d < 0] = 0
                usuelleSVD$d <- sqrt(usuelleSVD$d)
                usuelleSVD$v <- bb$vec[, 1:ncp]
                usuelleSVD$u <- t(t(crossprod(t(X), usuelleSVD$v))/
                                    usuelleSVD$d[1:ncp])
            }
        }
        U <- usuelleSVD$u
        V <- usuelleSVD$v
        if (ncp > 1) {
            mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
            mult[mult == 0] <- 1
            U <- t(t(U) * mult)
            V <- t(t(V) * mult)
        }
        U <- U/sqrt(rWeights)
        V <- V/sqrt(cWeights)
    }
    else {
        usuelleSVD <- tryCatchWE(svd(t(X), nu=ncp, nv=ncp))$val
        if (names(usuelleSVD)[[1]] == "message") {
            usuelleSVD <- tryCatchWE(svd(X, nu=ncp, nv=ncp))$val
            if (names(usuelleSVD)[[1]] == "d") {
                aux <- usuelleSVD$u
                usuelleSVD$u <- usuelleSVD$v
                usuelleSVD$v <- aux
            }
            else {
                bb <- eigen(crossprod(t(X), t(X)), symmetric=TRUE)
                usuelleSVD <- vector(mode="list", length=3)
                usuelleSVD$d[usuelleSVD$d < 0] <- 0
                usuelleSVD$d <- sqrt(usuelleSVD$d)
                usuelleSVD$v <- bb$vec[, 1:ncp]
                usuelleSVD$u <- t(t(crossprod(X, usuelleSVD$v)) / 
                                    usuelleSVD$d[1:ncp])
            }
        }
        U <- usuelleSVD$v
        V <- usuelleSVD$u
        mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
        mult[mult == 0] <- 1
        V <- t(t(V) * mult)/sqrt(cWeights)
        U <- t(t(U) * mult)/sqrt(rWeights)
    }
    vs <- usuelleSVD$d[1:min(ncol(X), nrow(X) - 1)]
    num <- which(vs[1:ncp] < 1e-15)
    if (length(num) == 1) {
        U[, num] <- U[, num, drop = FALSE] * vs[num]
        V[, num] <- V[, num, drop = FALSE] * vs[num]
    }
    if (length(num) > 1) {
        U[, num] <- t(t(U[, num]) * vs[num])
        V[, num] <- t(t(V[, num]) * vs[num])
    }
    res <- list(vs = vs, U = U, V = V)
    return(res)
}
