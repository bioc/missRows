searchsComb <- function(args, idx) {

  nargs <- length(args)
  nx <- vapply(args, length, 1)
  rep.fac <- c(1, cumprod(nx)[-nargs])

  iArgs <- seq_len(nargs)
  id <- ceiling(idx/rep.fac[iArgs]) %% nx[iArgs]
  id0 <- (id == 0)
  id[id0] <- nx[id0]

  comb <- NULL
  for (i in iArgs) {
    comb <- c(comb, args[[i]][id[i]])
  }

  comb
}
