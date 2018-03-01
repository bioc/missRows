tuneM <- function(object, ncomp = 2, M.max = 30, inc = 5, N = 10, tol = 1e-06,
                  show.plot = TRUE) {

  #-- initialization of variables ------------------------------------------------#
  #-------------------------------------------------------------------------------#

  #-- object
  if (class(object) != "MIDTList") {
    stop("'object' must be of class 'MIDTList'.", call. = FALSE)
  }

  #-- strata
  strt <- strata(object)
  strtLevels <- levels(strt)

  #-- datasets
  datasets <- incompleteData(object)
  datasets <- lapply(datasets, function(x, ind) {rownames(x) <- ind; return(x)},
                     names(strt))
  names(datasets) <- paste0("data", seq(length(datasets)))

  nb.cols <- sapply(datasets, ncol)
  nb.tables <- length(datasets)
  nb.rows <- length(strt)

  #-- checking general input parameters ------------------------------------------#
  #-------------------------------------------------------------------------------#

  #-- ncomp
  if (is.null(ncomp) || !is.numeric(ncomp) || (ncomp < 2) || !is.finite(ncomp)) {
    stop("invalid number of components, 'ncomp'.", call. = FALSE)
  }

  ncomp <- round(ncomp)

  #-- M.max (number max of imputations)
  if (is.null(M.max) || !is.numeric(M.max) || M.max < 1 || !is.finite(M.max)) {
    stop("invalid maximum number of imputations, 'M.max'.", call. = FALSE)
  }

  #-- increment of M
  if (is.null(inc) || !is.numeric(inc) || inc < 1 || !is.finite(inc)) {
    stop("invalid increment of M, 'inc'.", call. = FALSE)
  }

  if (inc > round(M.max/2 + 0.5)) {
    stop("'inc' must be less than or equal to ", round(M.max/2 + 0.5), ".",
         call. = FALSE)
  }

  M <- N * M.max

  #-- end checking ---------------------------------------------------------------#


  #-- creation of posible imputations in each stratum of each data ---------------#
  #-------------------------------------------------------------------------------#
  perm <- miss.row <- list()
  k <- 1
  id.data <- NULL

  for (j in seq_along(datasets)) {

    if (any(apply(is.na(datasets[[j]]), 1, function(x) { all(x) }))) {
      perm[[k]] <- miss.row[[k]] <- list()
      id.data <- c(id.data, names(datasets)[j])
      i <- 1

      for (s in seq_along(strtLevels)) {
        id.stratum <- (strt == strtLevels[s])
        tmp <- apply(is.na(datasets[[j]][strt == strtLevels[s], ]), 1, all)

        if (any(tmp)) {
          donors <- setdiff(names(tmp), names(tmp)[tmp])
          tmp2 <- t(permutations(length(donors), sum(tmp), donors))
          perm[[k]][[i]] <- tmp2 ## permutations per stratum
          miss.row[[k]][[i]]<- names(tmp)[tmp] ## missing rows
          i <- i + 1
        }
      }
      k <- k + 1
    }
  }

  #-- number of posible imputations
  tmp <- lapply(perm, function(x) vapply(x, ncol, 1))
  nb.miss.str <- unlist(lapply(tmp, length))
  id.data.miss <- list()
  from <- 1

  for (i in seq_along(nb.miss.str)) {
    to <- sum(nb.miss.str[1:i])
    id.data.miss[[i]] <- seq(from, to)
    from <- to + 1
  }

  M.total <- prod(unlist(tmp))

  #-- cheking N * M.max < M.total
  if (M > M.total) {
    stop("'N * M.max' must be less than ", M.total, ".", call. = FALSE)
  }

  seq.perm.data <- alply(matrix(unlist(tmp)), 1, seq)

  #-- selection of the donor indexes ---------------------------------------------#
  #-------------------------------------------------------------------------------#
  M.idx <- sample.int(min(M.total, 1e15), M)
  id.donor <- NULL

  for (i in seq_along(M.idx)) {
    id.donor <- rbind(id.donor, searchsComb(seq.perm.data, M.idx[i]))
  }

  #-- iterative approach ---------------------------------------------------------#
  #-------------------------------------------------------------------------------#
  variates.MFA <- list()
  Ml <- seq(inc, M.max, by = inc)
  nbMl <- length(Ml)
  ave.RV.coef <- sd.RV.coef <- NULL

  #-- initial configuration M_0
  conf0 <- list()
  m <- 1
  In <- split(1:(m * inc * N), rep(1:N, length = m * inc * N))

  #-- realisation of the MFA on the imputated data
  for (i in unlist(In)) {
    imput.data <- datasets

    for (j in seq_along(id.data)) {
      k <- 1
      for (s in id.data.miss[[j]]) {
        imput.ind <- perm[[j]][[k]][, id.donor[i, s]]

        #-- create imputate data
        imput.data[[id.data[j]]][miss.row[[j]][[k]], ] <-
          datasets[[id.data[j]]][imput.ind, ]
        k <- k + 1
      }
    }

    #-- realisation of the MFA
    result <- MFA(imput.data, ncomp, nb.rows, nb.tables, nb.cols)
    variates.MFA[[i]] <- data.frame(result$U)
  }

  for (n in 1:N) {
    #-- calculation of the compromise space (STATIS method)
    conf0[[n]] <- STATIS(variates.MFA[In[[n]]], nf = ncomp)$C.li
  }

  #-- configurations for M_l, l > 1
  old.ave <- -1
  m <- 2

  repeat {
    In <- split(1:(m * inc * N), rep(1:N, length = m * inc * N))
    subIn <- split((((m - 1) * inc * N) + 1):(m * inc * N),
                   rep(1:N, length = inc * N))
    RV.coef <- NULL

    #-- realisation of the MFA on the imputated data
    for (i in unlist(subIn)) {
      imput.data <- datasets

      for (j in seq_along(id.data)) {
        k <- 1
        for (s in id.data.miss[[j]]) {
          imput.ind <- perm[[j]][[k]][, id.donor[i, s]]

          #-- create imputate data
          imput.data[[id.data[j]]][miss.row[[j]][[k]], ] <-
            datasets[[id.data[j]]][imput.ind, ]
          k <- k + 1
        }
      }

      #-- realisation of the MFA
      result <- MFA(imput.data, ncomp, nb.rows, nb.tables, nb.cols)
      variates.MFA[[i]] <- data.frame(result$U)
    }

    for (n in 1:N) {
      #-- calculation of the compromise space (STATIS method)
      conf <- STATIS(variates.MFA[In[[n]]], nf = ncomp)$C.li
      RV.coef <- c(RV.coef, RVcoeff(conf0[[n]], conf))
      conf0[[n]] <- conf
    }

    ave.RV.coef <- c(ave.RV.coef, mean(RV.coef))
    sd.RV.coef <- c(sd.RV.coef, sqrt(var(RV.coef)))

    if (m >= nbMl | abs(old.ave - ave.RV.coef[m - 1]) < tol) break

    old.ave <- ave.RV.coef[m - 1]
    m <- m + 1
  }

  #-- graphic representation ------------------------------------------------------#
  #--------------------------------------------------------------------------------#
  df <- data.frame(x = 1:length(ave.RV.coef), avg = ave.RV.coef, sd = sd.RV.coef)
  lab <- paste0("(", Ml[1:(nbMl - 1)], ",", Ml[2:nbMl], ")")
  res <- list(stats = df[, -1])
  res$stats <- data.frame(imputations = lab, res$stats)

  g <- ggplot(df, aes(x = df$x, y = df$avg)) +
    geom_point(size = 2.5) + theme_bw() +
    geom_errorbar(aes(ymax = df$avg + df$sd, ymin = df$avg - df$sd), width = 0.15) +
    theme(panel.spacing = unit(2, "lines")) +
    labs(x = expression(paste("number of imputations (",
                              italic(M[l]), ", ", italic(M[l + 1]), ")")),
         y = 'RV coefficient\n') +
    scale_x_continuous(breaks = seq(nbMl - 1), labels = lab) +
    theme(axis.title = element_text(size = 16)) +
    theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1)) +
    theme(axis.text.y = element_text(size = 14))

  if (show.plot) { print(g) }

  #-- results --------------------------------------------------------------------#
  #-------------------------------------------------------------------------------#
  res <- c(res, list(ggp = g))
  class(res) <- "tuneM"
  return(invisible(res))
}
