imputedData <- function(object) {

  #-- checking general input arguments -------------------------------------------#
  #-------------------------------------------------------------------------------#

  #-- object is of 'MIDTList' S4 class
  if (class(object) != "MIDTList") {
    stop("'object' must be an object of class 'MIDTList'.",
         call. = FALSE)
  }

  #-- end checking ---------------------------------------------------------------#

  miss <- missingRows(object)
  X <- incompleteData(object)
  imput <- imputedRows(object)

  for (nm in names(missingRows(object))) {
    X[[nm]][miss[[nm]], ] <- imput[[nm]][miss[[nm]], ]
  }

  return(X)
}
