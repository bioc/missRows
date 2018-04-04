imputedData <- function(object) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- object is of 'MIDTList' S4 class
    if (class(object) != "MIDTList") {
        stop("'object' must be an object of class 'MIDTList'.",
            call.=FALSE)
    }
    
    ##- end checking ---------------------------------------------------------#
    
    ##- internal function for inserting columns in a matrix ------#
    ##------------------------------------------------------------#
    insertCols <- function(mat, id, cols) {
        
        mat <- as.data.frame(mat)
        
        for(i in seq_along(id)) {
            colSeq <- seq(from=id[i], to=ncol(mat))
            mat[, colSeq + 1] <- mat[, colSeq]
            mat[, id[i]] <- cols[, i]
        }
        
        return(as.matrix(mat))
    }
    ##------------------------------------------------------------#
    
    X <- assays(object)
    dfmap <- sampleMap(object)
    dfmap <- mapToList(dfmap, "assay")
    cnames <- rownames(colData(object))
    
    for (i in names(X)) {
        colnames(X[[i]]) <- dfmap[[i]]$primary
        idx <- (cnames %in% colnames(X[[i]]))
        X[[i]] <- X[[i]][, cnames[idx]]
        missColId <- which(!idx)
        X[[i]] <- insertCols(X[[i]], missColId, imputedIndv(object)[[i]])
        colnames(X[[i]]) <- cnames
    }
    
    return(X)
}
