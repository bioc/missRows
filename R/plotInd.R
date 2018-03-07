plotInd <- function(object,
                    comp = 1:2,
                    col.strata = NULL,
                    col.missing = "white",
                    conf.areas = c('none', 'ellipse', 'convex.hull'),
                    conf.level = 0.95,
                    ellipse.type = c("norm", "t"),
                    alpha = 0.1,
                    lwd = 0.3,
                    cex = 3,
                    legend.title = "Strata") {

  #-- checking general input arguments ---------------------------------------#
  #---------------------------------------------------------------------------#

  #-- object
  if (class(object) != "MIDTList") {
    stop("'object' must be an object of class 'MIDTList'", call. = FALSE)
  }

  if (is.null(object@compromise)) {
    stop("No 'compromise' slot found in the MIDTList object. Run MIMFA first",
         call. = FALSE)
  }

  #-- comp
  if (is.null(object@MIparam)) {
    stop("No 'MIparam' slot found in the MIDTList object. Run MIMFA first",
         call. = FALSE)
  }

  ncomp <- MIparam(object)$ncomp

  if (length(comp) != 2) {
    stop("'comp' must be a vector of length equal to 2", call. = FALSE)
  } else {
    if (any(!is.finite(comp)))
      stop("the elements of 'comp' must be positive integers", call. = FALSE)

    if (!is.numeric(comp) || any(comp < 1))
      stop("the elements of 'comp' must be positive integers", call. = FALSE)

    if (any(comp > ncomp))
      stop("the elements of 'comp' must be smaller or equal than ", ncomp, ".",
           call. = FALSE)
  }

  comp <- round(comp)

  #-- internal function for character color checking -----------#
  #-------------------------------------------------------------#
  isColor <- function(x) { sapply(x, function(x) {
    tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE) })
  }
  #-------------------------------------------------------------#

  #-- col.strata
  if (is.null(col.strata)) {
    col.strata <- rainbow(length(levels(strata(object))))
    names(col.strata) <- levels(strata(object))
  } else {
    if (length(col.strata) != length(levels(strata(object)))) {
      stop("'col.strata' must be a color names vector of length ",
           length(levels(strata(object))), ".", call. = FALSE)
    } else {
      if (any(!isColor(col.strata))) {
        stop("'col.strata' must be a character vector of recognized colors.",
             call. = FALSE)
      }
    }
  }

  if (is.null(names(col.strata))) {
    names(col.strata) <- levels(strata(object))
  } else {
    if (any(!(names(col.strata) %in% levels(strata(object)))))
      stop("names of 'col.strata' must be a character from: ",
           toString(levels(strata(object))), call. = FALSE)
  }

  #-- col.missing
  if (length(as.vector(col.missing)) != 1) {
    stop("'col.missing' must be a character of recognized colors.",
         call. = FALSE)
  }

  if (!isColor(col.missing) | is.na(col.missing)) {
    stop("'col.missing' must be a character of recognized colors.",
         call. = FALSE)
  }

  #-- conf.areas
  choices <- c('none', 'ellipse', 'convex.hull')
  conf.areas <- conf.areas[1]
  conf.areas <- choices[pmatch(conf.areas, choices)]

  if (is.na(conf.areas)) {
    stop("'conf.areas' should be one of 'none', 'ellipse' or 'convex.hull'.",
         call. = FALSE)
  }

  #-- ellipse type
  choices <- c("norm", "t")
  ellipse.type <- ellipse.type[1]
  ellipse.type <- choices[pmatch(ellipse.type, choices)]

  if (is.na(ellipse.type)) {
    stop("'ellipse.type' should be one of 'norm' or 't'.",
         call. = FALSE)
  }

  #-- alpha
  if (length(alpha) != 1) {
    stop("'alpha' must be of length 1", call. = FALSE)
  }

  if ((alpha < 0) | (alpha > 1)) {
    stop("alpha transparency value ", alpha, ", outside the interval [0,1]",
         call. = FALSE)
  }

  #-- lwd
  if (length(lwd) != 1) {
    stop("'lwd' must be of length 1", call. = FALSE)
  }

  if (lwd < 0) {
    stop("'lwd' must be zero or a positive number",
         call. = FALSE)
  }

  #-- cex
  if (!is.numeric(cex)) {
    stop("'cex' must be a positive number", call. = FALSE)
  }

  if (length(cex) != 1) {
    stop("'cex' must be of length 1", call. = FALSE)
  }

  if (cex <= 0) {
    stop("'cex' must be a positive number", call. = FALSE)
  }

  #-- legend.title
  legend.title <- as.graphicsAnnot(legend.title)

  #-- end checking -----------------------------------------------------------#


  #-- individuals scatter plot -----------------------------------------------#
  #---------------------------------------------------------------------------#
  comprConf <- compromise(object)
  n <- nrow(comprConf)
  miss <- rep("not", n)
  miss[names(strata(object)) %in% unlist(missingRows(object))] <- "yes"
  pch <- 21

  #-- none confidence areas --------------------------------------------------#
  #---------------------------------------------------------------------------#
  if (conf.areas == 'none') {
    df <- data.frame(x = comprConf[, comp[1]], y = comprConf[, comp[2]],
                     ind = names(strata(object)), stratum = strata(object),
                     missing = miss)

    df$ind.miss <- paste(df$ind, df$missing, sep = ".")
    df$ind.miss <- as.factor(df$ind.miss)

    ind.cols <- col.strata[df$stratum]

    g <- ggplot() + theme_bw() +
      geom_hline(yintercept = 0, color = 'grey30', size = 0.5, linetype = 2) +
      geom_vline(xintercept = 0, color = 'grey30', size = 0.5, linetype = 2) +
      geom_point(data = NULL, aes(x = df$x[df$missing == "not"],
                                  y = df$y[df$missing == "not"],
                                  fill = df$stratum[df$missing == "not"],
                                  color = df$stratum[df$missing == "not"]),
                 size = cex, shape = pch) +
      geom_point(data = NULL, aes(x = df$x[df$missing == "yes"],
                                  y = df$y[df$missing == "yes"],
                                  color = df$stratum[df$missing == "yes"]),
                 size = cex, shape = pch, fill = col.missing) +
      scale_fill_manual(name = legend.title, values = ind.cols) +
      guides(colour = "none") +
      labs(x = paste0('Comp ', comp[1]), y = paste0('Comp ', comp[2]))
  }

  #-- confidence ellipses ----------------------------------------------------#
  #---------------------------------------------------------------------------#
  if (conf.areas == 'ellipse') {
    df <- data.frame(x = comprConf[, comp[1]], y = comprConf[, comp[2]],
                     ind = names(strata(object)), stratum = strata(object),
                     missing = miss, conf = rep("compromise", n))

    m <- length(configurations(object))

    for (j in 1:m) {
      X <- as.matrix(configurations(object)[[j]][, comp])
      P <- X %*% solve(t(X) %*% X) %*% t(X)
      traj <- P %*% as.matrix(comprConf[, comp])

      temp <- data.frame(x = traj[, comp[1]], y = traj[, comp[2]],
                         ind = names(strata(object)), stratum = strata(object),
                         missing = miss, conf = rep("imputed", n))
      df <- rbind(df, temp)
    }

    df$ind.conf <- paste(df$ind, df$conf, sep = ".")
    df$ind.conf[df$conf == "compromise"] <-
      as.character(df$stratum[df$conf == "compromise"])
    df$ind.conf <- as.factor(df$ind.conf)

    ind.cols <- col.strata[df$stratum]
    names(ind.cols) <- df$ind.conf

    g <- ggplot() + theme_bw() +
      geom_hline(yintercept = 0, color = 'grey30', size = 0.5, linetype = 2) +
      geom_vline(xintercept = 0, color = 'grey30', size = 0.5, linetype = 2) +
      stat_ellipse(data = NULL, 
                   mapping = aes(x = df$x[df$conf == "imputed"],
                                 y = df$y[df$conf == "imputed"],
                                 fill = df$ind.conf[df$conf == "imputed"],
                                 color = df$ind.conf[df$conf == "imputed"]),
                   geom = 'polygon', type = ellipse.type, alpha = alpha, 
                   size = lwd) +
      geom_point(data = NULL,
                 aes(x = df$x[df$missing == "not" & df$conf == "compromise"],
                     y = df$y[df$missing == "not" & df$conf == "compromise"],
                     fill = df$stratum[df$missing == "not" & 
                                         df$conf == "compromise"],
                     color = df$stratum[df$missing == "not" & 
                                          df$conf == "compromise"]),
                 size = cex, shape = pch) +
      geom_point(data = NULL,
                 aes(x = df$x[df$missing == "yes" & df$conf == "compromise"],
                     y = df$y[df$missing == "yes" & df$conf == "compromise"],
                     color = df$stratum[df$missing == "yes" & 
                                          df$conf == "compromise"]),
                 size = cex, shape = pch, fill = col.missing) +
      scale_colour_manual(breaks = df$ind.cong[df$conf == "compromise"], 
                          values = ind.cols) +
      scale_fill_manual(name = legend.title, 
                        breaks = df$ind.conf[df$conf == "compromise"],
                        values = ind.cols) +
      labs(x = paste0('Comp ', comp[1]), y = paste0('Comp ', comp[2]))
  }

  #-- convex hulls -----------------------------------------------------------#
  #---------------------------------------------------------------------------#
  if (conf.areas == 'convex.hull') {
    df <- data.frame(x = comprConf[, comp[1]], y = comprConf[, comp[2]],
                     ind = names(strata(object)), stratum = strata(object),
                     missing = miss, conf = rep("compromise", n))

    m <- length(configurations(object))

    for (j in 1:m) {
      X <- as.matrix(configurations(object)[[j]][, comp])
      P <- X %*% solve(t(X) %*% X) %*% t(X)
      traj <- P %*% as.matrix(comprConf[, comp])

      temp <- data.frame(x = traj[, comp[1]], y = traj[, comp[2]],
                         ind = names(strata(object)), stratum = strata(object),
                         missing = miss, conf = rep("imputed", n))
      df <- rbind(df, temp)
    }

    df$ind.conf <- paste(df$ind, df$conf, sep = ".")
    df$ind.conf[df$conf == "compromise"] <-
      as.character(df$stratum[df$conf == "compromise"])
    df$ind.conf <- as.factor(df$ind.conf)

    find_hull <- function(df) df[chull(df$x, df$y), ]
    hulls <- ddply(df, "ind", find_hull)

    ind.cols <- col.strata[df$stratum]
    names(ind.cols) <- df$ind.conf

    g <- ggplot() + theme_bw() +
      geom_hline(yintercept = 0, color = 'grey30', size = 0.5, linetype = 2) +
      geom_vline(xintercept = 0, color = 'grey30', size = 0.5, linetype = 2) +
      geom_polygon(data = NULL,
                   aes(x = hulls$x[hulls$conf == "imputed"],
                       y = hulls$y[hulls$conf == "imputed"],
                       fill = hulls$ind.conf[hulls$conf == "imputed"],
                       color = hulls$ind.conf[hulls$conf == "imputed"]),
                   alpha = alpha, size = lwd) +
      geom_point(data = NULL,
                 aes(x = df$x[df$missing == "not" & df$conf == "compromise"],
                     y = df$y[df$missing == "not" & df$conf == "compromise"],
                     fill = df$stratum[df$missing == "not" & 
                                         df$conf == "compromise"],
                     color = df$stratum[df$missing == "not" & 
                                          df$conf == "compromise"]),
                 size = cex, shape = pch) +
      geom_point(data = NULL,
                 aes(x = df$x[df$missing == "yes" & df$conf == "compromise"],
                     y = df$y[df$missing == "yes" & df$conf == "compromise"],
                     color = df$stratum[df$missing == "yes" & 
                                          df$conf == "compromise"]),
                 size = cex, shape = pch, fill = col.missing) +
      scale_colour_manual(breaks = df$ind.cong[df$conf == "compromise"], 
                          values = ind.cols) +
      scale_fill_manual(name = legend.title, 
                        breaks = df$ind.conf[df$conf == "compromise"],
                        values = ind.cols) +
      labs(x = paste0('Comp ', comp[1]), y = paste0('Comp ', comp[2]))
  }

  print(g)
  return(invisible(g))
}
