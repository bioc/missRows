MIMFA <- function(object, ncomp = 2, M = NULL, estim.ncp = FALSE,
                  max.iter = 500, tol = 1e-10) {


  #-- checking general input arguments -------------------------------------------#
  #-------------------------------------------------------------------------------#

  #-- object
  if (class(object) != "MIDTList") {
    stop("'object' must be an object of class 'MIDTList'.", call. = FALSE)
  }

  #-- ncomp
  if (is.null(ncomp) || !is.numeric(ncomp) || (ncomp < 2) || !is.finite(ncomp)) {
    stop("invalid number of components, 'ncomp'. It must be a positive integer",
         call. = FALSE)
  }

  if (ncomp > min(length(strata(object)) - 1,
                  sum(sapply(incompleteData(object), ncol)))) {
    stop("'ncomp' must be a numeric value lower or equal to ",
         min(length(strata(object)) - 1, sum(sapply(incompleteData(object), ncol))),
         call. = FALSE)
  }

  ncomp <- round(ncomp)

  #-- M (number of imputations)
  if (is.null(M) || !is.numeric(M) || M <= 1 || !is.finite(M)) {
    stop("invalid number of imputations, 'M'. It must be a positive integer",
         call. = FALSE)
  }

  M <- round(M)

  #-- estim.ncp
  if (length(estim.ncp) != 1 || !is.logical(estim.ncp)) {
    stop("'estim.ncp' must be either TRUE or FALSE", call. = FALSE)
  }

  #-- max.iter
  if (length(max.iter) != 1) {
    stop("'max.iter' must be a single value", call. = FALSE)
  }

  if (is.null(max.iter) || max.iter < 1 || !is.finite(max.iter) ||
      !is.numeric(max.iter)) {
    stop("invalid value for 'max.iter'. It must be a positive integer",
         call. = FALSE)
  }

  max.iter <- round(max.iter)

  #-- tol
  if (length(tol) != 1) {
    stop("'tol' must be a single value", call. = FALSE)
  }

  if (is.null(tol) || tol < 0 || !is.finite(tol)) {
    stop("invalid value for 'tol'. It must be a positive number",
         call. = FALSE)
  }

  #-- end checking ---------------------------------------------------------------#


  #-- initialization of variables ------------------------------------------------#
  #-------------------------------------------------------------------------------#
  incompData <- incompleteData(object)
  strt <- strata(object)
  missRows <- missingRows(object)

  nb.cols <- sapply(incompData, ncol)
  nb.tables <- length(incompData)
  nb.rows <- length(strt)
  str.levels <- levels(strt)

  #-- creation of posible imputations in each stratum of each data ---------------#
  #-------------------------------------------------------------------------------#
  perm <- strMissRows <- strMiss <- vector("list", length(missRows))
  nmMissRows <- names(missRows)
  names(perm) <- names(strMissRows) <- names(strMiss) <- nmMissRows

  for (nm in names(strMissRows)) {
    perm[[nm]] <- strMissRows[[nm]] <- list()

    for (s in str.levels) {
      tmp <- apply(is.na(incompData[[nm]][strt == s, ]), 1, all)

      if (any(tmp)) {
        donors <- setdiff(names(tmp), names(tmp)[tmp])
        ## permutations per stratum
        perm[[nm]][[s]] <- t(permutations(length(donors), sum(tmp), donors))
        ## missing rows per stratum
        strMissRows[[nm]][[s]] <- names(tmp)[tmp]
      }

    }
  }

  rm(tmp)

  #-- number of posible imputations
  tmp <- lapply(perm, function(x) vapply(x, ncol, 1))
  nbStrMiss <- unlist(lapply(tmp, length))
  from <- 1

  for (i in seq_along(nbStrMiss)) {
    to <- sum(nbStrMiss[1:i])
    strMiss[[i]] <- seq(from, to)
    from <- to + 1
  }

  M.total <- prod(unlist(tmp))
  seq.perm.data <- alply(matrix(unlist(tmp)), 1, seq)

  if (is.null(M)) { M <- 30 }

  M <- min(M, M.total)

  #-- selection of the donor indexes ---------------------------------------------#
  #-------------------------------------------------------------------------------#
  M.idx <- sample.int(min(M.total, 1e15), M)
  id.donor <- NULL

  for (i in seq_along(M.idx)) {
    id.donor <- rbind(id.donor, searchsComb(seq.perm.data, M.idx[i]))
  }

  #-- realisation of the MFA on the imputed data ---------------------------------#
  #-------------------------------------------------------------------------------#
  U <- list()
  center <- sigma <- matrix(0, nrow = M, ncol = sum(nb.cols))

  for (m in 1:M) { ## nb. of imputations M
    completeData <- incompData

    for (nm in nmMissRows) {
      k <- 1

      for (s in strMiss[[nm]]) { ## nb. of strata with missing rows
        imputInd <- perm[[nm]][[k]][, id.donor[m, s]]

        #-- create imputate data
        completeData[[nm]][strMissRows[[nm]][[k]], ] <-
          incompData[[nm]][imputInd, ]
        k <- k + 1
      }
    }

    #-- realisation of the MFA
    result <- MFA(completeData, ncomp, nb.rows, nb.tables, nb.cols)
    U[[m]] <- data.frame(result$U)
  }

  rm(perm)

  #-- calculation of the compromise space (STATIS method) ------------------------#
  #-------------------------------------------------------------------------------#
  tmp <- wrapperSTATIS(U, nf = ncomp)
  colnames(tmp$C.li) <- paste0("comp ", 1:ncomp)

  #-- estimation of the number of components for data imputation -----------------#
  #-------------------------------------------------------------------------------#
  if (estim.ncp) {
    ncp <- estimNCP(tmp$C.ro, min.ncp = 2, ncomp)
    attr(ncp, "estimated") <- TRUE
  } else {
    ncp <- ncomp
    attr(ncp, "estimated") <- FALSE
  }

  #-- data imputation ------------------------------------------------------------#
  #-------------------------------------------------------------------------------#
  impD <- imputeDataMFA(incompData, tmp$C.li, missRows, comp = 1:ncp,
                        max.iter = max.iter, tol = tol)

  #-- results: MIDTList S4 class -------------------------------------------------#
  #-------------------------------------------------------------------------------#
  object <- new("MIDTList",
                object,
                compromise = tmp$C.li[, 1:ncp],
                configurations = U,
                imputedRows = impD,
                MIparam = list(method = "MFA",
                               ncomp = ncp,
                               M = M,
                               M.total = M.total))

  return(object)
}
