MFA <- function(dataTables, ncomp, nb.rows, nb.tables, nc.table) {
    
    weights <- rep(0, nb.tables)
    K.mat <- matrix(0, nb.rows, 0)
    center <- sigma <- NULL
    
    for (k in 1:nb.tables) {
        scaled.table <- scale(dataTables[[k]]) * sqrt(nb.rows/(nb.rows - 1))
        ## calculer sd, if(sd < 1e-08) sd <- 1 puis reduir
        center <- c(center, attr(scaled.table, "scaled:center"))
        sigma <- c(sigma, attr(scaled.table, "scaled:scale") * 
                    sqrt((nb.rows - 1)/nb.rows))
        sigma[sigma < 1e-08] <- 1
        weights[k] <- 1/(eigenvalue(scaled.table))^2
        K.mat <- cbind(K.mat, scaled.table)
    }
    
    ## global PCA
    tmp <- wrapperSVD(K.mat, col.w = rep(weights, nc.table), ncp = ncomp)
    
    eig <- tmp$vs[1:ncomp]
    U <- sweep(tmp$U, 2, eig, "*")
    
    res <- list(U = U)
    return(res)
}
