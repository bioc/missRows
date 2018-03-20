plotInd <- function(object,
                    comp=1:2,
                    colStrata=NULL,
                    colMissing="white",
                    confAreas=c('none', 'ellipse', 'convex.hull'),
                    confLevel=0.95,
                    ellipseType=c("norm", "t"),
                    alpha=0.1,
                    lwd=0.3,
                    cex=3,
                    legTitle="Strata") {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- object
    if (class(object) != "MIDTList") {
        stop("'object' must be an object of class 'MIDTList'", call.=FALSE)
    }
    
    if (is.null(object@compromise)) {
        stop("No 'compromise' slot found in the MIDTList object.",
            " Run MIMFA first", call.=FALSE)
    }
    
    ##- comp
    if (is.null(object@MIparam)) {
        stop("No 'MIparam' slot found in the MIDTList object.",
            " Run MIMFA first", call.=FALSE)
    }
    
    ncomp <- MIparam(object)$ncomp
    
    if (length(comp) != 2) {
        stop("'comp' must be a vector of length equal to 2", call.=FALSE)
    } else {
        if (any(!is.finite(comp)))
            stop("the elements of 'comp' must be positive integers", 
                call.=FALSE)
        
        if (!is.numeric(comp) || any(comp < 1))
            stop("the elements of 'comp' must be positive integers", 
                call.=FALSE)
        
        if (any(comp > ncomp))
            stop("the elements of 'comp' must be smaller or equal than ", 
                ncomp, ".", call.=FALSE)
    }
    
    comp <- round(comp)
    
    ##- internal function for character color checking -----------#
    ##------------------------------------------------------------#
    isColor <- function(x) { vapply(x, function(x) {
        tryCatch(is.matrix(col2rgb(x)), error=function(e) FALSE) },
        TRUE)
    }
    ##------------------------------------------------------------#
    
    ##- colStrata
    if (is.null(colStrata)) {
        colStrata <- rainbow(length(levels(strata(object))))
        names(colStrata) <- levels(strata(object))
    } else {
        if (length(colStrata) != length(levels(strata(object)))) {
            stop("'colStrata' must be a color names vector of length ",
                length(levels(strata(object))), ".", call.=FALSE)
        } else {
            if (any(!isColor(colStrata))) {
                stop("'colStrata' must be a character vector of recognized",
                    " colors.", call.=FALSE)
            }
        }
    }
    
    if (is.null(names(colStrata))) {
        names(colStrata) <- levels(strata(object))
    } else {
        if (any(!(names(colStrata) %in% levels(strata(object)))))
            stop("names of 'colStrata' must be a character from: ",
                toString(levels(strata(object))), call.=FALSE)
    }
    
    ##- colMissing
    if (length(as.vector(colMissing)) != 1) {
        stop("'colMissing' must be a character of recognized colors.",
            call.=FALSE)
    }
    
    if (!isColor(colMissing) | is.na(colMissing)) {
        stop("'colMissing' must be a character of recognized colors.",
            call.=FALSE)
    }
    
    ##- confAreas
    choices <- c('none', 'ellipse', 'convex.hull')
    confAreas <- confAreas[1]
    confAreas <- choices[pmatch(confAreas, choices)]
    
    if (is.na(confAreas)) {
        stop("'confAreas' should be one of 'none', 'ellipse' or",
            " 'convex.hull'.", call.=FALSE)
    }
    
    ##- confLevel
    
    
    ##- ellipse type
    choices <- c("norm", "t")
    ellipseType <- ellipseType[1]
    ellipseType <- choices[pmatch(ellipseType, choices)]
    
    if (is.na(ellipseType)) {
        stop("'ellipseType' should be one of 'norm' or 't'.", call.=FALSE)
    }
    
    ##- alpha
    if (length(alpha) != 1) {
        stop("'alpha' must be of length 1", call.=FALSE)
    }
    
    if ((alpha < 0) | (alpha > 1)) {
        stop("alpha transparency value ", alpha, 
            ", outside the interval [0,1]", call.=FALSE)
    }
    
    ##- lwd
    if (length(lwd) != 1) {
        stop("'lwd' must be of length 1", call.=FALSE)
    }
    
    if (lwd < 0) {
        stop("'lwd' must be zero or a positive number", call.=FALSE)
    }
    
    ##- cex
    if (!is.numeric(cex)) {
        stop("'cex' must be a positive number", call.=FALSE)
    }
    
    if (length(cex) != 1) {
        stop("'cex' must be of length 1", call.=FALSE)
    }
    
    if (cex <= 0) {
        stop("'cex' must be a positive number", call.=FALSE)
    }
    
    ##- legTitle
    legTitle <- as.graphicsAnnot(legTitle)
    
    ##- end checking ---------------------------------------------------------#
    
    
    ##- individuals scatter plot ---------------------------------------------#
    ##------------------------------------------------------------------------#
    comprConf <- compromise(object)
    n <- nrow(comprConf)
    miss <- rep("not", n)
    miss[names(strata(object)) %in% unlist(missingRows(object))] <- "yes"
    pch <- 21
    
    ##- none confidence areas ------------------------------------------------#
    ##------------------------------------------------------------------------#
    if (confAreas == 'none') {
        df <- data.frame(x = comprConf[, comp[1]], y = comprConf[, comp[2]],
                        ind = names(strata(object)), stratum = strata(object),
                        missing = miss)
        
        df$ind.miss <- paste(df$ind, df$missing, sep = ".")
        df$ind.miss <- as.factor(df$ind.miss)
        
        indCols <- colStrata[df$stratum]
        
        g <- ggplot() + theme_bw() +
            geom_hline(yintercept=0, color='grey30', size=0.5, 
                        linetype=2) +
            geom_vline(xintercept=0, color='grey30', size=0.5, 
                        linetype=2) +
            geom_point(data=NULL, aes(x=df$x[df$missing == "not"],
                                    y=df$y[df$missing == "not"],
                                    fill=df$stratum[df$missing == "not"],
                                    color=df$stratum[df$missing == "not"]),
                        size=cex, shape=pch) +
            geom_point(data=NULL, aes(x=df$x[df$missing == "yes"],
                                    y=df$y[df$missing == "yes"],
                                    color=df$stratum[df$missing == "yes"]),
                        size=cex, shape=pch, fill=colMissing) +
            scale_fill_manual(name=legTitle, values=indCols) +
            guides(colour = "none") +
            labs(x=paste0('Comp ', comp[1]), y=paste0('Comp ', comp[2]))
    }
    
    ##- confidence ellipses --------------------------------------------------#
    ##------------------------------------------------------------------------#
    if (confAreas == 'ellipse') {
        df <- data.frame(x = comprConf[, comp[1]], y = comprConf[, comp[2]],
                        ind = names(strata(object)), stratum = strata(object),
                        missing = miss, conf = rep("compromise", n))
        
        m <- length(configurations(object))
        
        for (j in seq_len(m)) {
            X <- as.matrix(configurations(object)[[j]][, comp])
            P <- X %*% solve(t(X) %*% X) %*% t(X)
            traj <- P %*% as.matrix(comprConf[, comp])
            
            temp <- data.frame(x = traj[, comp[1]], y = traj[, comp[2]],
                                ind = names(strata(object)), 
                                stratum = strata(object),
                                missing = miss, conf = rep("imputed", n))
            df <- rbind(df, temp)
        }
        
        df$ind.conf <- paste(df$ind, df$conf, sep = ".")
        df$ind.conf[df$conf == "compromise"] <-
            as.character(df$stratum[df$conf == "compromise"])
        df$ind.conf <- as.factor(df$ind.conf)
        
        indCols <- colStrata[df$stratum]
        names(indCols) <- df$ind.conf
        
        g <- ggplot() + theme_bw() +
            geom_hline(yintercept=0, color='grey30', size=0.5, 
                        linetype=2) +
            geom_vline(xintercept=0, color='grey30', size=0.5, 
                        linetype=2) +
            stat_ellipse(data=NULL, 
                        mapping=aes(x=df$x[df$conf == "imputed"],
                                    y=df$y[df$conf == "imputed"],
                                    fill=df$ind.conf[df$conf == "imputed"],
                                    color=df$ind.conf[df$conf == "imputed"]),
                        geom='polygon', type=ellipseType, level=confLevel,
                        alpha=alpha, size=lwd) +
            geom_point(data=NULL,
                        aes(x=df$x[df$missing == "not" & 
                                        df$conf == "compromise"],
                            y=df$y[df$missing == "not" & 
                                        df$conf == "compromise"],
                            fill=df$stratum[df$missing == "not" &
                                                df$conf == "compromise"],
                            color=df$stratum[df$missing == "not" & 
                                                df$conf == "compromise"]),
                        size=cex, shape=pch) +
            geom_point(data=NULL,
                        aes(x=df$x[df$missing == "yes" & 
                                        df$conf == "compromise"],
                            y=df$y[df$missing == "yes" & 
                                        df$conf == "compromise"],
                            color=df$stratum[df$missing == "yes" & 
                                                df$conf == "compromise"]),
                        size=cex, shape=pch, fill=colMissing) +
            scale_colour_manual(breaks=df$ind.cong[df$conf == "compromise"], 
                                values=indCols) +
            scale_fill_manual(name=legTitle, 
                                breaks=df$ind.conf[df$conf == "compromise"],
                                values=indCols) +
            labs(x=paste0('Comp ', comp[1]), y=paste0('Comp ', comp[2]))
    }
    
    ##- convex hulls ---------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if (confAreas == 'convex.hull') {
        df <- data.frame(x = comprConf[, comp[1]], y = comprConf[, comp[2]],
                        ind = names(strata(object)), stratum = strata(object),
                        missing = miss, conf = rep("compromise", n))
        
        m <- length(configurations(object))
        
        for (j in seq_len(m)) {
            X <- as.matrix(configurations(object)[[j]][, comp])
            P <- X %*% solve(t(X) %*% X) %*% t(X)
            traj <- P %*% as.matrix(comprConf[, comp])
            
            temp <- data.frame(x = traj[, comp[1]], y = traj[, comp[2]],
                                ind = names(strata(object)), 
                                stratum = strata(object),
                                missing = miss, conf = rep("imputed", n))
            df <- rbind(df, temp)
        }
        
        df$ind.conf <- paste(df$ind, df$conf, sep = ".")
        df$ind.conf[df$conf == "compromise"] <-
                as.character(df$stratum[df$conf == "compromise"])
        df$ind.conf <- as.factor(df$ind.conf)
        
        find_hull <- function(df) df[chull(df$x, df$y), ]
        hulls <- ddply(df, "ind", find_hull)
        
        indCols <- colStrata[df$stratum]
        names(indCols) <- df$ind.conf
        
        g <- ggplot() + theme_bw() +
            geom_hline(yintercept=0, color='grey30', size=0.5, 
                        linetype=2) +
            geom_vline(xintercept=0, color='grey30', size=0.5, 
                        linetype=2) +
            geom_polygon(data=NULL,
                        aes(x=hulls$x[hulls$conf == "imputed"],
                            y=hulls$y[hulls$conf == "imputed"],
                            fill=hulls$ind.conf[hulls$conf == "imputed"],
                            color=hulls$ind.conf[hulls$conf == "imputed"]),
                        alpha=alpha, size=lwd) +
            geom_point(data=NULL,
                        aes(x=df$x[df$missing == "not" & 
                                        df$conf == "compromise"],
                            y=df$y[df$missing == "not" & 
                                        df$conf == "compromise"],
                            fill=df$stratum[df$missing == "not" & 
                                                df$conf == "compromise"],
                            color=df$stratum[df$missing == "not" & 
                                                df$conf == "compromise"]),
                        size=cex, shape=pch) +
            geom_point(data=NULL,
                        aes(x=df$x[df$missing == "yes" & 
                                        df$conf == "compromise"],
                            y=df$y[df$missing == "yes" & 
                                        df$conf == "compromise"],
                            color=df$stratum[df$missing == "yes" & 
                                                df$conf == "compromise"]),
                        size=cex, shape=pch, fill=colMissing) +
            scale_colour_manual(breaks=df$ind.cong[df$conf == "compromise"], 
                                values=indCols) +
            scale_fill_manual(name=legTitle, 
                                breaks=df$ind.conf[df$conf == "compromise"],
                                values=indCols) +
            labs(x=paste0('Comp ', comp[1]), y=paste0('Comp ', comp[2]))
    }
    
    print(g)
    return(invisible(g))
}
