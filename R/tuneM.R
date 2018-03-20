tuneM <- function(object, ncomp=2, Mmax=30, inc=5, N=10, tol=1e-06,
                    showPlot=TRUE) {
    
    ##- initialization of variables ------------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- object
    if (class(object) != "MIDTList") {
        stop("'object' must be of class 'MIDTList'.", call.=FALSE)
    }
    
    ##- strata
    strt <- strata(object)
    strtLevels <- levels(strt)
    
    ##- datasets
    datasets <- incompleteData(object)
    datasets <- lapply(datasets, function(x, ind) {rownames(x) <- ind; x},
                        names(strt))
    names(datasets) <- paste0("data", seq(length(datasets)))
    
    nbCols <- vapply(datasets, ncol, 1L)
    nbTables <- length(datasets)
    nbRows <- length(strt)
    
    ##- checking general input parameters ------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- ncomp
    if (is.null(ncomp) || !is.numeric(ncomp) || (ncomp < 2) || 
        !is.finite(ncomp)) {
        stop("invalid number of components, 'ncomp'.", call.=FALSE)
    }
    
    ncomp <- round(ncomp)
    
    ##- Mmax (number max of imputations)
    if (is.null(Mmax) || !is.numeric(Mmax) || Mmax < 1 || 
        !is.finite(Mmax)) {
        stop("invalid maximum number of imputations, 'Mmax'.", call.=FALSE)
    }
    
    ##- increment of M
    if (is.null(inc) || !is.numeric(inc) || inc < 1 || !is.finite(inc)) {
        stop("invalid increment of M, 'inc'.", call.=FALSE)
    }
    
    if (inc > round(Mmax/2 + 0.5)) {
        stop("'inc' must be less than or equal to ", round(Mmax/2 + 0.5), ".",
            call.=FALSE)
    }
    
    M <- N * Mmax
    
    ##- end checking ---------------------------------------------------------#
    
    
    ##- creation of posible imputations in each stratum of each data ---------#
    ##------------------------------------------------------------------------#
    perm <- missRow <- list()
    k <- 1
    idData <- NULL
    
    for (j in seq_along(datasets)) {
        
        if (any(apply(is.na(datasets[[j]]), 1, function(x) { all(x) }))) {
            perm[[k]] <- missRow[[k]] <- list()
            idData <- c(idData, names(datasets)[j])
            i <- 1
            
            for (s in seq_along(strtLevels)) {
                id.stratum <- (strt == strtLevels[s])
                tmp <- apply(is.na(datasets[[j]][strt == strtLevels[s], ]), 
                            1, all)
                
                if (any(tmp)) {
                    donors <- setdiff(names(tmp), names(tmp)[tmp])
                    tmp2 <- t(permutations(length(donors), sum(tmp), donors))
                    perm[[k]][[i]] <- tmp2 ## permutations per stratum
                    missRow[[k]][[i]]<- names(tmp)[tmp] ## missing rows
                    i <- i + 1
                }
            }
            k <- k + 1
        }
    }
    
    ##- number of posible imputations
    tmp <- lapply(perm, function(x) vapply(x, ncol, 1))
    nbMissStr <- unlist(lapply(tmp, length))
    idDataMiss <- list()
    from <- 1
    
    for (i in seq_along(nbMissStr)) {
        to <- sum(nbMissStr[seq_len(i)])
        idDataMiss[[i]] <- seq(from, to)
        from <- to + 1
    }
    
    Mtotal <- prod(unlist(tmp))
    
    ##- cheking N * Mmax < Mtotal
    if (M > Mtotal) {
        stop("'N * Mmax' must be less than ", Mtotal, ".", call.=FALSE)
    }
    
    seqPermData <- alply(matrix(unlist(tmp)), 1, seq)
    
    ##- selection of the donor indexes ---------------------------------------#
    ##------------------------------------------------------------------------#
    Midx <- sample.int(min(Mtotal, 1e15), M)
    idDonor <- NULL
    
    for (i in seq_along(Midx)) {
        idDonor <- rbind(idDonor, searchsComb(seqPermData, Midx[i]))
    }
    
    ##- iterative approach ---------------------------------------------------#
    ##------------------------------------------------------------------------#
    variatesMFA <- list()
    Ml <- seq(inc, Mmax, by=inc)
    nbMl <- length(Ml)
    aveRVcoef <- sdRVcoef <- NULL
    
    ##- initial configuration M_0
    conf0 <- list()
    m <- 1
    In <- split(seq_len(m * inc * N), rep(seq_len(N), length = m * inc * N))
    
    ##- realisation of the MFA on the imputated data
    for (i in unlist(In)) {
        imputData <- datasets
        
        for (j in seq_along(idData)) {
            k <- 1
            for (s in idDataMiss[[j]]) {
                imputInd <- perm[[j]][[k]][, idDonor[i, s]]
                
                ##- create imputate data
                imputData[[idData[j]]][missRow[[j]][[k]], ] <-
                    datasets[[idData[j]]][imputInd, ]
                k <- k + 1
            }
        }
        
        ##- realisation of the MFA
        result <- MFA(imputData, ncomp, nbRows, nbTables, nbCols)
        variatesMFA[[i]] <- data.frame(result$U)
    }
    
    for (n in seq_len(N)) {
        ##- calculation of the compromise space (STATIS method)
        conf0[[n]] <- STATIS(variatesMFA[In[[n]]], nf=ncomp)$Cli
    }
    
    ##- configurations for M_l, l > 1
    oldAve <- -1
    m <- 2
    
    repeat {
        In <- split(seq_len(m * inc * N), rep(seq_len(N), length=m * inc * N))
        subIn <- split((((m - 1) * inc * N) + 1):(m * inc * N),
                        rep(seq_len(N), length=inc * N))
        RV <- NULL
        
        ##- realisation of the MFA on the imputated data
        for (i in unlist(subIn)) {
            imputData <- datasets
            
            for (j in seq_along(idData)) {
                k <- 1
                for (s in idDataMiss[[j]]) {
                    imputInd <- perm[[j]][[k]][, idDonor[i, s]]
                    
                    ##- create imputate data
                    imputData[[idData[j]]][missRow[[j]][[k]], ] <-
                        datasets[[idData[j]]][imputInd, ]
                    k <- k + 1
                }
            }
            
            ##- realisation of the MFA
            result <- MFA(imputData, ncomp, nbRows, nbTables, nbCols)
            variatesMFA[[i]] <- data.frame(result$U)
        }
        
        for (n in seq_len(N)) {
            ##- calculation of the compromise space (STATIS method)
            conf <- STATIS(variatesMFA[In[[n]]], nf=ncomp)$Cli
            RV <- c(RV, RVcoeff(conf0[[n]], conf))
            conf0[[n]] <- conf
        }
        
        aveRVcoef <- c(aveRVcoef, mean(RV))
        sdRVcoef <- c(sdRVcoef, sqrt(var(RV)))
        
        if (m >= nbMl | abs(oldAve - aveRVcoef[m - 1]) < tol) break
        
        oldAve <- aveRVcoef[m - 1]
        m <- m + 1
    }
    
    ##- graphic representation -----------------------------------------------#
    ##------------------------------------------------------------------------#
    df <- data.frame(x = seq_along(aveRVcoef), avg = aveRVcoef, 
                    sd = sdRVcoef)
    lab <- paste0("(", Ml[seq_len(nbMl - 1)], ",", Ml[2:nbMl], ")")
    res <- list(stats = df[, -1])
    res$stats <- data.frame(imputations = lab, res$stats)
    
    g <- ggplot(df, aes(x=df$x, y=df$avg)) +
        geom_point(size=2.5) + theme_bw() +
        geom_errorbar(aes(ymax=df$avg + df$sd, ymin=df$avg - df$sd), 
                    width=0.15) +
        theme(panel.spacing=unit(2, "lines")) +
        labs(x=expression(paste("number of imputations (",
                            italic(M[l]), ", ", italic(M[l + 1]), ")")),
            y = 'RV coefficient\n') +
        scale_x_continuous(breaks=seq(nbMl - 1), labels=lab) +
        theme(axis.title=element_text(size=16)) +
        theme(axis.text.x=element_text(size=14, angle=45, hjust=1, 
                                        vjust=1)) +
        theme(axis.text.y=element_text(size=14))
    
    if (showPlot) { print(g) }
    
    ##- results --------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    res <- c(res, list(ggp = g))
    class(res) <- "tuneM"
    return(invisible(res))
}
