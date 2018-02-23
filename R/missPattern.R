missPattern <- function(object, col.strata = NULL, col.missing = "grey70",
                        cex.ttitle = 12, legend.title = "Strata",
                        miss.lab = "miss", show.plot = TRUE) {

  #-- checking general input arguments -------------------------------------------#
  #-------------------------------------------------------------------------------#

  #-- object
  if (class(object) != "MIDTList") {
    stop("'object' must be an object of class 'MIDTList'.", call. = FALSE)
  }

  #-- col.strata

  #-- internal function for character color checking -----------#
  #-------------------------------------------------------------#
  isColor <- function(x) { sapply(x, function(x) {
    tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE) })
  }
  #-------------------------------------------------------------#

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
    stop("'col.missing' must be a single character of recognized colors.",
         call. = FALSE)
  }

  if (!isColor(col.missing) | is.na(col.missing)) {
    stop("'col.missing' must be a single character of recognized colors.",
         call. = FALSE)
  }

  #-- end checking ---------------------------------------------------------------#


  #-- initialization of variables ------------------------------------------------#
  #-------------------------------------------------------------------------------#
  missData <- incompleteData(object)
  strt <- strata(object)
  nmTables <- tableNames(object)

  mrp <- matrix(nrow = length(strt), ncol = length(missData))
  rownames(mrp) <- names(strt)
  colnames(mrp) <- nmTables

  nbmr <- matrix(nrow = length(levels(strt)), ncol = length(missData))
  rownames(nbmr) <- levels(strt)
  colnames(nbmr) <- nmTables

  #-- finds the missing data pattern ---------------------------------------------#
  #-------------------------------------------------------------------------------#
  for (i in 1:length(missData)) {
    id.miss <- apply(is.na(missData[[i]]), 1, all)
    tmp <- split(id.miss, strt)
    nbmr[, i] <- sapply(tmp, sum)
    tmp <- unlist(split(id.miss, strt))
    mrp[, i] <- tmp
  }

  nbmr <- rbind(nbmr, c(colSums(nbmr)))
  nbmr <- cbind(nbmr, c(rowSums(nbmr)))
  colnames(nbmr)[ncol(nbmr)] <- "    "
  mrp <- data.frame(mrp, strata = strt)

  #-- plot missing rows pattern --------------------------------------------------#
  #-------------------------------------------------------------------------------#
  ymax <- cumsum(table(mrp$strata))
  ymin <- c(0, ymax[-length(ymax)])

  ymin.m <- 0:(length(strt) - 1)
  ymax.m <- 1:length(strt)

  df <- data.frame(ymin = NA, ymax = NA, strata = NA, Table = NA)
  df.m <- data.frame(ymin = 0, ymax = 0, strata = "", Table = nmTables[1])

  for (i in 1:length(missData)) {
    tmp <- data.frame(ymin = ymin, ymax = ymax, strata = levels(strt),
                      Table = rep(nmTables[i], length(levels(strt))))
    df <- rbind(df, tmp)

    id <- mrp[, i]
    tmp <- data.frame(ymin = ymin.m, ymax = ymax.m,
                      strata = rep(miss.lab, nrow(mrp)),
                      Table = rep(nmTables[i], nrow(mrp)))
    tmp <- tmp[id, ]
    df.m <- rbind(df.m, tmp)
  }

  df <- df[-1, ]
  df$strata <- factor(df$strata,
                      levels = c(sort(unique(df$strata)), "", miss.lab))
  df.m$strata <- factor(df.m$strata,
                        levels = c("", miss.lab))
  df$Table <- factor(df$Table, levels = nmTables)
  df.m$Table <- factor(df.m$Table, levels = nmTables)

  ind.cols <- c(col.strata, "transparent", col.missing)
  names(ind.cols)[length(col.strata) + 1:2] <- c("", miss.lab)

  g <- ggplot() + theme_bw() + facet_wrap(~ Table) +
    scale_y_reverse(expand = c(0, 1.5)) +
    geom_rect(data = df, aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax,
                             fill = strata), colour = "white") +
    geom_rect(data = df.m, aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax,
                               fill = strata),
              colour = "white") +
    scale_fill_manual(values = ind.cols, breaks = names(ind.cols),
                      name = legend.title) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_text(size = cex.ttitle)) +
    ggtitle("Missingness pattern") +
    theme(plot.title = element_text(hjust = 0.5))

  if (show.plot) { print(g) }

  #-- results --------------------------------------------------------------------#
  #-------------------------------------------------------------------------------#
  res <- list(nbMissing = nbmr, isMissing = mrp, ggp = g)
  class(res) = "missPattern"
  return(invisible(res))
}
