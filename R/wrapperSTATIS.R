wrapperSTATIS <- function (Ktab, nf = 3, tol = 1e-07) {
  
  nlig <- nrow(Ktab[[1]]) 
  ntab <- length(Ktab)
  result <- sep <- list()
  
  sep <- lapply(Ktab, function(wk) { wk <- as.matrix(wk)
                                     wk <- t(wk) * sqrt(1/ncol(wk))
                                     wk <- t(wk) %*% wk
                                     wk })
  
  sep <- matrix(unlist(sep), nlig * nlig, ntab)
  
  RV <- t(sep) %*% sep
  ak <- sqrt(diag(RV))
  RV <- sweep(RV, 1, ak, "/")
  RV <- sweep(RV, 2, ak, "/")
  
  eig1 <- eigen(RV, symmetric = TRUE)
  
  if (any(eig1$vectors[, 1] < 0)) {
    eig1$vectors[, 1] <- -eig1$vectors[, 1]
  }
  
  tabw <- eig1$vectors[, 1]
  
  sep <- t(t(sep)/ak)
  C.ro <- rowSums(sweep(sep, 2, tabw, "*"))
  C.ro <- matrix(unlist(C.ro), nlig, nlig)
  eig1 <- eigen(C.ro, symmetric = TRUE)
  result$C.ro <- C.ro
  rm(C.ro)
  eig <- eig1$values
  rank <- sum((eig/eig[1]) > tol)
  
  if (nf <= 0) { nf <- 2 }
  
  if (nf > rank) { nf <- rank }
  
  wref <- eig1$vectors[, 1:nf]
  rm(eig1)
  w <- data.frame(t(t(wref) * sqrt(eig[1:nf])))
  names(w) <- paste("C", 1:nf, sep = "")
  result$C.li <- w
  
  return(result)
}
