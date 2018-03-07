plotVar <- function(object, comp = 1:2, col = NULL, var.names = FALSE,
                    cex = 3, pch = 19, alpha = 0.7, spty = TRUE, cutoff = 0,
                    rad.in = 0.5, overlap = TRUE, ncols = 2,
                    legend.title = "Tables") {

  #-- checking general input arguments & initialization of variables --------#
  #--------------------------------------------------------------------------#

  #-- object
  if (class(object) != "MIDTList") {
    stop("'object' must be an object of class 'MIDTList'.", call. = FALSE)
  }

  if (is.null(object@compromise)) {
    stop("No 'compromise' slot found in the MIDTList object. Run MIMFA first",
         call. = FALSE)
  }

  completeData <- imputedData(object)
  comprConf <- compromise(object)
  dtab <- 1:length(completeData)
  nmTables <- tableNames(object)
  nbTables <- length(completeData)

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

  #-- col

  #-- internal function for character color checking -----------#
  #-------------------------------------------------------------#
  isColor <- function(x) { sapply(x, function(x) {
    tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE) })
  }
  #-------------------------------------------------------------#

  #-- brewer palette 'Set1' ------------------------------------#
  brewerPal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                 "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  #-------------------------------------------------------------#

  if (is.null(col)) {
    col <- brewerPal[dtab]
    names(col) <- nmTables
  } else {
    if (length(col) != nbTables) {
      stop("'col' must be a color names vector of length ",
           nbTables, ".", call. = FALSE)
    } else {
      if (any(!isColor(col))) {
        stop("'col' must be a character vector of recognized colors.",
             call. = FALSE)
      }
    }
  }

  if (is.null(names(col))) {
    names(col) <- nmTables
  } else {
    if (any(!(names(col) %in% nmTables))) {
      stop("names of 'col' must be a character from: ",
           toString(nmTables), call. = FALSE)
    }
  }

  #-- var.names
  flag.var.names <- rep(FALSE, nbTables)

  if (is.logical(var.names)) {
    if (var.names) {
      var.names <- nmTables
      flag.var.names <- rep(TRUE, nbTables)
    }
  } else {
    if (is.character(var.names)) {
      if (!any(var.names %in% nmTables)) {
        stop("One element of 'var.names' does not match with the names ",
             "of the data tables. \n'var.names' should be one of: ",
             toString(nmTables), call. = FALSE)
      } else {
        id <- (nmTables %in% var.names)
        flag.var.names[id] <- TRUE
      }
    } else {
      stop("Incorrect value for 'var.names'", call. = FALSE)
    }
  }

  #-- cex
  if (length(cex) == 1) {
    if (is.null(cex) || !is.numeric(cex) || !is.finite(cex) || (cex < 0)) {
      stop("'cex' must be a positive number", call. = FALSE)
    }
    cex <- rep(cex, nbTables)
  } else {
    if (length(cex) > nbTables) {
      stop("'cex' must be a numeric vector of either length 1 or length ",
           nbTables, call. = FALSE)
    }
    if (any(is.null(cex)) | any(!is.numeric(cex)) | any(!is.finite(cex)) |
        any(cex < 0)) {
      stop("'cex' must be a numeric vector of either length 1 or length ",
           nbTables, call. = FALSE)
    }
  }

  #-- pch
  if (length(pch) == 1) {
    pch <- rep(pch, nbTables)
  } else {
    if (length(pch) > nbTables) {
      stop("'pch' must be a vector of either length 1 or length ",
           nbTables, call. = FALSE)
    }
  }

  #-- alpha
  if (length(alpha) != 1) {
    stop("'alpha' must be a single value", call. = FALSE)
  }

  if (is.null(alpha) || !is.numeric(alpha) || !is.finite(alpha)) {
    stop("'alpha' transparency value must be numeric", call. = FALSE)
  }

  if ((alpha < 0) || (alpha > 1)) {
    stop("'alpha' transparency value ", alpha, ", outside the interval [0,1]",
         call. = FALSE)
  }

  #-- spty
  if (length(spty) != 1 || !is.logical(spty)) {
    stop("'spty' must be either TRUE or FALSE", call. = FALSE)
  }

  #-- cutoff
  if (length(cutoff) != 1) {
    stop("'cutoff' must be a single value", call. = FALSE)
  }

  if (is.null(cutoff) || !is.numeric(cutoff) || !is.finite(cutoff)) {
    stop("'cutoff' value must be numeric", call. = FALSE)
  }

  if ((cutoff > 1) || (cutoff < 0)) {
    stop("'cutoff' value ", cutoff, ", outside the interval [0,1]",
         call. = FALSE)
  }

  #-- rad.in
  if (length(rad.in) != 1) {
    stop("'rad.in' must be a single value", call. = FALSE)
  }

  if (is.null(rad.in) || !is.numeric(rad.in) || !is.finite(rad.in)) {
    stop("'rad.in' value must be numeric", call. = FALSE)
  }

  if ((rad.in > 1) || (rad.in < 0)) {
    stop("'rad.in' value ", rad.in, ", outside the interval [0,1]",
         call. = FALSE)
  }

  #-- overlap
  if (length(overlap) != 1 || !is.logical(overlap)) {
    stop("'overlap' must be a logical value", call. = FALSE)
  }

  #-- ncols
  if (length(ncols) != 1) {
    stop("'ncols' must be of length 1", call. = FALSE)
  }

  if (is.null(ncols) || !is.numeric(ncols) || !is.finite(ncols) ||
      ncols < 1) {
    stop("invalid number of columns 'ncols' for facetting data tables",
         call. = FALSE)
  }

  #-- legend.title
  legend.title <- as.graphicsAnnot(legend.title)

  #-- end checking -----------------------------------------------------------#


  #-- variables scatter plot -------------------------------------------------#
  #---------------------------------------------------------------------------#
  #-- individuals coordinates --#
  coord <- lapply(dtab, 
                  function(x) { cor(completeData[[x]], comprConf[, comp]) })
  idx <- vector("list", nbTables)

  if (cutoff > 0) {
    for (i in 1:(nbTables - 1)) {
      for (j in 2:nbTables) {
        mat <- coord[[i]] %*% t(coord[[j]])
        tmp <- (abs(mat) >= cutoff)

        if (all(!tmp)) {
          stop("'cutoff' must be < ", round(max(abs(mat)), 2), call. = FALSE)
        }

        tmp <- which(tmp, arr.ind = TRUE)
        idx[[i]] <- rbind(idx[[i]], tmp[, 1, drop = FALSE])
        idx[[j]] <- rbind(idx[[j]], tmp[, 2, drop = FALSE])
      }
    }

    idx <- lapply(idx, function(x) unique(x))
    namevar <- NULL

    for (i in 1:nbTables) {
      coord[[i]] <- matrix(coord[[i]][idx[[i]], ], ncol = length(comp))
      namevar <- c(namevar, colnames(completeData[[i]])[idx[[i]]])
    }
  }

  #-- data frame for ggplot
  nbVar <- sapply(coord, nrow)
  df <- do.call(rbind, coord)
  if (cutoff == 0) { namevar <- rownames(df) }
  rownames(df) <- NULL

  df <- data.frame(df, namevar, rep(nmTables, nbVar))
  names(df) <- c("x", "y", "Var", "Tables")
  df$Tables <- factor(df$Tables, levels = nmTables)

  x <- y <- Tables <- Var <- NULL
  circle <- data.frame(x = cos(seq(0, 2 * pi, l = 1000)),
                       y = sin(seq(0, 2 * pi, l = 1000)))

  circle.int <- data.frame(x = rad.in * cos(seq(0, 2 * pi, l = 1000)),
                           y = rad.in * sin(seq(0, 2 * pi, l = 1000)))

  #-- carried out ggplot
  g <- ggplot(df, aes(x = x, y = y)) + theme_bw() +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +
    geom_path(data = circle, aes(x = x, y = y), color = "black") +
    geom_path(data = circle.int, aes(x = x, y = y), color = "black") +
    geom_vline(aes(xintercept = 0), linetype = 2, colour = "gray40") +
    geom_hline(aes(yintercept = 0), linetype = 2, colour = "gray40")

  if (spty) { g <- g + coord_equal(ratio = 1) }

  for (i in 1:nbTables) {
    id <- (df$Tables == nmTables[i])

    if (flag.var.names[i]) {
      g <- g + geom_text(data = df[id, ],
                         aes(x = x, y = y, colour = Tables, label = Var),
                         size = cex[i], alpha = alpha, show.legend = FALSE) +
        geom_point(data = df[id, ][1, ],
                   aes(x = 0, y = 0, colour = Tables),
                   shape = pch[i], size = 0, alpha = alpha)
    } else {
      g <- g + geom_point(data = df[id, ],
                          aes(x = x, y = y, colour = Tables),
                          shape = pch[i], size = cex[i], alpha = alpha)
    }
  }

  g <- g + scale_colour_manual(name = legend.title, values = col,
                               breaks = nmTables) +
    guides(colour = guide_legend(override.aes = list(shape = pch, 
                                                     size = cex))) +
    labs(x = paste0('Comp ', comp[1]), y = paste0('Comp ', comp[2])) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))

  if (!overlap) {
    g <- g + facet_wrap(~ Tables, ncol = ncols, as.table = TRUE) +
      theme(strip.text = element_text(size = 12))
  }

  print(g)
  return(invisible(list(df = df, ggp = g)))
}
