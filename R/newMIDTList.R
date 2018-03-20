newMIDTList <- function(..., strata=NULL, tableNames=NULL) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- input tables
    K <- list(...)
    if (is.list(K[[1]])) { K <- K[[1]] }
    
    if (length(K) <= 1) {
        stop("at least two data tables must be passed as arguments in '...'",
            call.=FALSE)
    }
    
    ##- strata
    if (is.null(strata)) {
        stop("'strata' must be a named vector or data frame", call.=FALSE)
    }
    
    if (is.data.frame(strata)) {
        tmp <- rownames(strata)
        strata <- strata[, 1]
        names(strata) <- tmp
    }
    
    ##- table names
    if (is.null(tableNames)) { tableNames <- names(K) }
    
    if (is.null(tableNames)) { tableNames <- paste0("Table ", 1:length(K)) }
    
    if (length(tableNames) != length(K)) {
        stop("the length of 'tableNames' must be equal to ", length(K),
            call.=FALSE)
    }
    
    if (any(duplicated(tableNames))) {
        stop("non-unique values in 'tableNames'", call.=FALSE)
    }
    
    names(K) <- tableNames
    
    ##- the individual tables
    for (i in 1:length(K)) {
        if (length(dim(K[[i]])) != 2) {
            stop("the '", names(K)[i], "' data table must be a matrix or",
                " data frame.", call.=FALSE)
        }
        
        if (!is.numeric(as.matrix(K[[i]]))) {
            stop("the '", names(K)[i], "' data table must be a matrix or",
                " data frame.", call.=FALSE)
        }
    }
    
    ### Inf values in tables
    for (i in 1:length(K)) {
        if (any(apply(K[[i]], 1, is.infinite))) {
            stop("infinite values in '", names(K)[i], "' data table.", 
                call.=FALSE)
        }
    }
    
    ##- equal row numbers among tables
    nr <- c(vapply(K, nrow, 1L), length(strata))
    if (length(unique(nr)) != 1)
        stop("non equal row numbers among tables and/or strata", 
            call.=FALSE)
    
    ##- the tables are rows named
    for (i in 1:length(K)) {
        if (is.null(rownames(K[[i]]))) {
            stop("the '", names(K)[i], "' data table must be rows named.", 
                call.=FALSE)
        }
    }
    
    ##- strata is named vector
    if (is.null(names(strata))) {
        stop("'strata' must be a named vector.", call.=FALSE)
    }
    
    ##- samples are ordered correctly
    rnames <- names(strata)
    for (i in 1:length(K)) {
        if (!identical(rnames, rownames(K[[i]]))) {
            stop("non equal row names among tables and/or strata", 
                call.=FALSE)
        }
    }
    
    strata <- as.factor(strata)
    
    ##- there are no missing individuals
    miss <- lapply(K, function(X) {
            rownames(X)[which(apply(is.na(X), 1, function(x) { all(x) }))] })
    
    if (length(unlist(miss)) == 0)
        stop("no missing rows in the data tables, MI is not useful.", 
            " Perform MFA.", call.=FALSE)
    
    miss <- miss[vapply(miss, length, 1L) > 0]
    
    ##- end checking ---------------------------------------------------------#
    
    
    ##- the S4 class ---------------------------------------------------------#
    ##------------------------------------------------------------------------#
    object <- MIDTList(incompleteData=K, strata=strata, 
                        tableNames=tableNames, missingRows=miss)
    return(object)
}
