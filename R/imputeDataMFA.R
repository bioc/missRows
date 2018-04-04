imputeDataMFA <- function(datasets, U, missRows, comp,
                        maxIter=500, tol=1e-10) {
    
    X <- do.call(cbind, datasets)
    id.na <- is.na(X)
    Xold <- X
    Xold[id.na] <- 0
    U <- as.matrix(U[, comp])
    Vold <- matrix(999, nrow=ncol(X), ncol=length(comp))
    
    iter <- 0
    
    repeat {
        Xold <- scale(Xold)
        center <- attr(Xold, "scaled:center")
        sigma <- attr(Xold, "scaled:scale")
        
        V <- t(solve(crossprod(U)) %*% t(U) %*% Xold)
        Xnew <- U %*% t(V)
        Xnew <- sweep(Xnew, 2, sigma, "*")
        Xnew <- sweep(Xnew, 2, center, "+")
        Xnew[!id.na] <- X[!id.na]
        
        if (sqrt(sum((V - Vold)^2)) < tol) break
        
        Xold <- Xnew
        
        if (iter == maxIter) {
            warning("maximum number of iterations reached in data imputation",
                    call.=FALSE)
            break
        }
        
        Vold <- V
        iter <- iter + 1
    }
    
    rownames(Xnew) <- rownames(X)
    
    imputData <- vector("list", length(missRows))
    names(imputData) <- names(missRows)
    nl <- c(0, unlist(lapply(datasets, ncol)))
    
    for (nm in names(imputData)) {
        id <- which(names(nl[-1]) == nm)
        from <- sum(nl[seq_len(id)]) + 1
        to <- sum(nl[seq_len(id + 1)])
        imputData[[nm]] <- t(Xnew[missRows[[nm]], from:to])
    }
    
    return(imputData)
}
