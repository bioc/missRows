imputeDataMFA <- function(datasets, U, missRows, comp,
                          max.iter = 500, tol = 1e-10) {

  X <- do.call(cbind, datasets)
  id.na <- is.na(X)
  X.old <- X
  X.old[id.na] <- 0
  U <- as.matrix(U[, comp])
  V.old <- matrix(999, nrow = ncol(X), ncol = length(comp))

  iter <- 0

  repeat {
    X.old <- scale(X.old)
    center <- attr(X.old, "scaled:center")
    sigma <- attr(X.old, "scaled:scale")

    V <- t(solve(crossprod(U)) %*% t(U) %*% X.old)
    X.new <- U %*% t(V)
    X.new <- sweep(X.new, 2, sigma, "*")
    X.new <- sweep(X.new, 2, center, "+")
    X.new[!id.na] <- X[!id.na]

    if (sqrt(sum((V - V.old)^2)) < tol) break

    X.old <- X.new

    if (iter == max.iter) {
      warning("maximum number of iterations reached in data imputation",
              call. = FALSE)
      break
    }

    V.old <- V
    iter <- iter + 1
  }

  rownames(X.new) <- rownames(X)

  imputData <- vector("list", length(missRows))
  names(imputData) <- names(missRows)
  nl <- c(0, unlist(lapply(datasets, ncol)))

  for (nm in names(imputData)) {
    id <- which(names(nl[-1]) == nm)
    from <- sum(nl[1:id]) + 1
    to <- sum(nl[1:(id + 1)])
    imputData[[nm]] <- X.new[missRows[[nm]], from:to]
  }

  return(imputData)
}
