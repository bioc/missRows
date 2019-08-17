##- Helper function for check general input arguments ------------------------#
##----------------------------------------------------------------------------#
checkInputArg <- function(dTables, cData, strata, assayNames) {
    
    ##- colData
    if (is.null(cData)) {
        stop("'colData' must be a DataFrame", call.=FALSE)
    }

    if (!is(cData, "DataFrame")) {
        stop("'colData' must be a DataFrame", call.=FALSE)
    } else {
        tmp <- rownames(cData)
        
        if (is.null(tmp)) {
            stop("'colData' must be a DataFrame rows named.", call.=FALSE)
        }
    }
    
    ##- strata
    if (is.null(strata) & (ncol(cData) == 1)) {
        strata <- colnames(cData)
    }
    
    if (is.null(strata)) {
        strata <- colnames(cData) <- "strata"
    }
    
    if (length(strata) > 1) {
        stop("'strata' must be a vector of length 1.", call.=FALSE)
    }
    
    if (is.character(strata)) {
        cond <- (strata %in% colnames(cData))

        if (!cond || length(cond) == 0) {
            stop("'strata' do not match colnames in colData.", call.=FALSE)
        }
    } else {
        if (strata %in% seq_len(ncol(cData))) {
            stop("'strata' must be a integer value less than or equal to ", 
                ncol(cData), call.=FALSE)
        }
    }
    
    ##- table names
    if (is.null(assayNames)) { assayNames <- names(dTables) }

    if (is.null(assayNames)) { 
        assayNames <- paste0("Table ", seq_along(dTables)) 
    }

    if (length(assayNames) != length(dTables)) {
        stop("the length of 'assayNames' must be equal to ", 
            length(dTables), call.=FALSE)
    }
    
    if (any(duplicated(assayNames))) {
        stop("non-unique names in 'assayNames'", call.=FALSE)
    }

    names(dTables) <- assayNames
    
    ##- the individual tables
    for (i in seq_along(dTables)) {
        if (length(dim(dTables[[i]])) != 2) {
            stop("the '", names(dTables)[i], "' data table must be a matrix",
                " or data frame.", call.=FALSE)
        }
        
        dTables[[i]] <- as.matrix(dTables[[i]])
        
        if (!is.numeric(dTables[[i]])) {
            stop("the '", names(dTables)[i], "' data table must be a matrix",
                " or data frame.", call.=FALSE)
        }
    }

    ##- the tables are columns named
    for (i in seq_along(dTables)) {
        if (is.null(colnames(dTables[[i]]))) {
            stop("the '", names(dTables)[i], 
                "' data table must be columns named.", call.=FALSE)
        }
    }
    
    ##- samples are named correctly
    cnames <- rownames(cData)
    
    for (i in seq_along(dTables)) {
        if (!(all(colnames(dTables[[i]]) %in% cnames))) {
            stop("columns names in '", names(dTables)[i], 
                "' data table do not match row names in colData.", call.=FALSE)
        }
    }

    ### Inf values in tables
    for (i in seq_along(dTables)) {
        if (any(apply(dTables[[i]], 1, is.infinite))) {
            stop("infinite values in '", names(dTables)[i], "' data table.", 
                call.=FALSE)
        }
    }
    
    ##- missing individuals
    cnames <- rownames(cData)
    miss <- vector("list", length = length(dTables))
    names(miss) <- names(dTables)
    
    for (i in names(dTables)) {
        idx <- !(cnames %in% colnames(dTables[[i]]))
        miss[[i]] <- cnames[idx]
    }
    
    miss <- miss[vapply(miss, length, 1L) > 0]
    
    return(list(tablesList = dTables, strata = strata, missingIndv = miss))
}
