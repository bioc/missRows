newMIDTList <- function(..., object) UseMethod("newMIDTList")


##- newMIDTList default method -----------------------------------------------#
##----------------------------------------------------------------------------#
newMIDTList.default <- function(..., colData=NULL, strata=NULL, 
                                assayNames=NULL) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- input tables
    K <- list(...)
    
    if (length(K) == 1) {
        K <- K[[1]]

        if (!is.list(K) || is.data.frame(K)) {
            stop("at least two data tables must be passed as arguments",
                " in '...'", call.=FALSE)
        }
    }

    if (length(K) <= 1) {
        stop("at least two data tables must be passed as arguments in '...'",
            call.=FALSE)
    }
    
    ##- colData
    if (is.null(colData)) {
        stop("'colData' must be a data frame", call.=FALSE)
    }
    
    if (!is.data.frame(colData)) {
        stop("'colData' must be a data frame", call.=FALSE)
    } else {
        tmp <- rownames(colData)
        
        if (is.null(tmp)) {
            stop("'colData' must be a data frame with row names.", call.=FALSE)
        }
        
        colData <- DataFrame(strata=colData[, 1], row.names=tmp)
    }
    
    ##- strata
    if (is.null(strata) & (ncol(colData) == 1)) {
        strata <- colnames(colData)
    }
    
    if (is.null(strata)) {
        strata <- colnames(colData) <- "strata"
    }
    
    if (length(strata) > 1) {
        stop("'strata' must be a vector of length 1.", call.=FALSE)
    }
    
    if (is.character(strata)) {
        cond <- (strata %in% colnames(colData))
        
        if (!cond || length(cond) == 0) {
            stop("'strata' do not match colnames in colData.", call.=FALSE)
        }
    } else {
        if (strata %in% seq_len(ncol(colData))) {
            stop("'strata' must be a integer value less than or equal to ", 
                ncol(colData), call.=FALSE)
        }
    }
    
    ##- table names
    if (is.null(assayNames)) { assayNames <- names(K) }
    
    if (is.null(assayNames)) { assayNames <- paste0("Table ", seq_along(K)) }
    
    if (length(assayNames) != length(K)) {
        stop("the length of 'assayNames' must be equal to ", length(K),
            call.=FALSE)
    }
    
    if (any(duplicated(assayNames))) {
        stop("non-unique names in 'assayNames'", call.=FALSE)
    }
    
    names(K) <- assayNames
    
    ##- the individual tables
    for (i in seq_along(K)) {
        if (length(dim(K[[i]])) != 2) {
            stop("the '", names(K)[i], "' data table must be a matrix or",
                " data frame.", call.=FALSE)
        }
        
        K[[i]] <- as.matrix(K[[i]])
        
        if (!is.numeric(K[[i]])) {
            stop("the '", names(K)[i], "' data table must be a matrix or",
                " data frame.", call.=FALSE)
        }
    }
    
    ### Inf values in tables
    for (i in seq_along(K)) {
        if (any(apply(K[[i]], 1, is.infinite))) {
            stop("infinite values in '", names(K)[i], "' data table.", 
                call.=FALSE)
        }
    }
    
    ##- the tables are columns named
    for (i in seq_along(K)) {
        if (is.null(colnames(K[[i]]))) {
            stop("the '", names(K)[i], "' data table must be columns named.", 
                call.=FALSE)
        }
    }
    
    ##- samples are named correctly
    cnames <- rownames(colData)
    
    for (i in seq_along(K)) {
        if (!(all(colnames(K[[i]]) %in% cnames))) {
            stop("columns names in '", names(K)[i], "' data table do not ",
                "match row names in colData.", call.=FALSE)
        }
    }
    
    ##- missing individuals
    miss <- list()
    
    for (i in seq_along(K)) {
        idx <- !(cnames %in% colnames(K[[i]]))
        miss[[i]] <- cnames[idx]
    }
    
    names(miss) <- names(K)
    miss <- miss[vapply(miss, length, 1L) > 0]
    
    ##- end checking ---------------------------------------------------------#
    
    
    ##- the MIDTList S4 class ------------------------------------------------#
    ##------------------------------------------------------------------------#
    midt <- MultiAssayExperiment(experiments=K, colData=colData)
    midt <- new("MIDTList", midt, strata=strata, missingIndv=miss)
    
    return(midt)
}


##- newMIDTList method for MultiAssayExperiment class ------------------------#
##----------------------------------------------------------------------------#
newMIDTList.MultiAssayExperiment <- function(object, strata=NULL, ...) {
    
    ##- initialization of variables ------------------------------------------#
    ##------------------------------------------------------------------------#
    K <- assays(object)
    
    ##- strata
    if (is.null(strata) & (ncol(colData(object)) == 1)) {
        strata <- colnames(colData(object))
        
        if (is.null(strata)) {
            strata <- colnames(colData(object)) <- "strata"
        }
    }
    
    if (length(strata) > 1) {
        stop("'strata' must be a vector of length 1.", call.=FALSE)
    }
    
    if (is.character(strata)) {
        cond <- (strata %in% colnames(colData(object)))
        
        if (!cond || length(cond) == 0) {
            stop("'strata' do not match colnames in colData.", call.=FALSE)
        }
    } else {
        if (strata %in% seq_len(ncol(colData(object)))) {
            stop("'strata' must be a integer value less than or equal to ", 
                ncol(colData(object)), call.=FALSE)
        }
    }
    
    cnames <- rownames(colData(object))
    
    ##- the individual tables
    for (i in seq_along(K)) {
        if (length(dim(K[[i]])) != 2) {
            stop("the '", names(K)[i], "' data table must be a matrix or",
                " data frame.", call.=FALSE)
        }
        
        K[[i]] <- as.matrix(K[[i]])
        
        if (!is.numeric(K[[i]])) {
            stop("the '", names(K)[i], "' data table must be a matrix or",
                " data frame.", call.=FALSE)
        }
    }
    
    ### Inf values in tables
    for (i in seq_along(K)) {
        if (any(apply(K[[i]], 1, is.infinite))) {
            stop("infinite values in '", names(K)[i], "' data table.", 
                call.=FALSE)
        }
    }
    
    ##- missing individuals
    miss <- vector("list", length = length(K))
    names(miss) <- names(K)
    mapList <- mapToList(sampleMap(object), "assay")
    
    for (i in names(K)) {
        idx <- !(cnames %in% mapList[[i]]$primary)
        miss[[i]] <- cnames[idx]
    }
    
    miss <- miss[vapply(miss, length, 1L) > 0]

    
    ##- the MIDTList S4 class ------------------------------------------------#
    ##------------------------------------------------------------------------#
    midt <- new("MIDTList", object, strata=strata, missingIndv=miss)
    
    return(midt)
}
