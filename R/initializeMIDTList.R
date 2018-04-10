##- MIDTList constructor from a list -----------------------------------------#
##----------------------------------------------------------------------------#
MIDTListFromTablesList <- function(tablesList, colData=NULL, strata=NULL, 
                                    assayNames=NULL) {
    

    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    resCheck <- checkInputArg(tablesList, colData, strata, assayNames)
        
        
    ##- output for MIDTList S4 class -----------------------------------------#
    ##------------------------------------------------------------------------#
    mae <- MultiAssayExperiment(ExperimentList(resCheck$tablesList), 
                                colData=colData)
    
    return(list(mae = mae, strata = resCheck$strata,
                missingIndv = resCheck$miss))
}


##- MIDTList constructor from a MultiAssayExperiment -------------------------#
##----------------------------------------------------------------------------#
MIDTListFromMultiAssayExperiment <- function(object, strata) {
    
    ##- initialization of variables ------------------------------------------#
    ##------------------------------------------------------------------------#
    tablesList <- assays(object)
    cData <- colData(object)
    cnames <- rownames(cData)
    mapList <- mapToList(sampleMap(object), "assay")
    
    for (i in names(tablesList)) {
        idx <- (cnames %in% mapList[[i]]$primary)
        colnames(tablesList[[i]]) <- cnames[idx]
    }
    

    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    resCheck <- checkInputArg(tablesList, cData, strata, assayNames=NULL)
        
        
    ##- output for MIDTList S4 class -----------------------------------------#
    ##------------------------------------------------------------------------#
    mae <- MultiAssayExperiment(ExperimentList(resCheck$tablesList), 
                                colData=cData)

    return(list(mae = mae, strata = resCheck$strata, 
                missingIndv = resCheck$missingIndv))
}
