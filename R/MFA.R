MFA <- function(dataTables, ncomp, nbRows, nbTables, ncTables) {
    
    weights <- rep(0, nbTables)
    Kmat <- matrix(0, nbRows, 0)
    center <- sigma <- NULL
    
    for (k in seq_len(nbTables)) {
        scaledTable <- scale(dataTables[[k]]) * sqrt(nbRows/(nbRows - 1))
        ## calculer sd, if(sd < 1e-08) sd <- 1 puis reduir
        center <- c(center, attr(scaledTable, "scaled:center"))
        sigma <- c(sigma, attr(scaledTable, "scaled:scale") * 
                    sqrt((nbRows - 1)/nbRows))
        sigma[sigma < 1e-08] <- 1
        weights[k] <- 1/(eigenvalue(scaledTable))^2
        Kmat <- cbind(Kmat, scaledTable)
    }
    
    ##- global PCA
    tmp <- wrapperSVD(Kmat, cWeights=rep(weights, ncTables), ncp=ncomp)
    
    eig <- tmp$vs[seq_len(ncomp)]
    U <- sweep(tmp$U, 2, eig, "*")
    
    res <- list(U = U)
    return(res)
}
