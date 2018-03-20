###############################################################################
### MIDTList S4 class definition
###############################################################################

##- setClass -----------------------------------------------------------------#
##----------------------------------------------------------------------------#
MIDTList <- setClass("MIDTList",
                    slots=c(incompleteData="ANY",
                    strata="ANY",
                    tableNames="ANY",
                    missingRows="ANY",
                    compromise="ANY",
                    configurations="ANY",
                    imputedRows="ANY",
                    MIparam="ANY"))


##- setGeneric ---------------------------------------------------------------#
##----------------------------------------------------------------------------#

##- incompleteData
setGeneric(name="incompleteData",
            def=function(object) standardGeneric("incompleteData"))

##- strata
setGeneric(name="strata",
            def=function(object) standardGeneric("strata"))

##- tableNames
setGeneric(name="tableNames",
            def=function(object, ...) standardGeneric("tableNames"))

setGeneric("tableNames<-",
            def =function(object, value) standardGeneric("tableNames<-"))

##- missingRows
setGeneric(name="missingRows",
            def=function(object) standardGeneric("missingRows"))

##- compromise
setGeneric(name="compromise",
            def=function(object) standardGeneric("compromise"))

##- configurations
setGeneric(name="configurations",
            def=function(object, ...) standardGeneric("configurations"))

##- imputedRows
setGeneric(name="imputedRows",
            def=function(object) standardGeneric("imputedRows"))

##- MIparam
setGeneric(name="MIparam",
            def=function(object) standardGeneric("MIparam"))



##- setMethod ----------------------------------------------------------------#
##----------------------------------------------------------------------------#
##- show
setMethod("show",
        signature="MIDTList",
        definition=function(object) {

            nbMiss <- NULL

            for (j in seq_along(object@incompleteData)) {
                idMiss <- apply(is.na(object@incompleteData[[j]]), 1, all)
                if (any(idMiss)) {
                    nbMiss <- c(nbMiss, sum(idMiss))
                } else {
                    nbMiss <- c(nbMiss, 0)
                }
            }

            nt <- length(object@incompleteData)

            cat("An object of class ", class(object), ".",
                "\n\nTables:\n", sep = "")
            info <- data.frame(names(object@incompleteData),
                            vapply(object@incompleteData, nrow, 1L),
                            vapply(object@incompleteData, ncol, 1L),
                            nbMiss,
                            row.names = paste0("Table ", 1:nt, " "))
            colnames(info) <- c("name", "rows", "columns", "missing")
            print(info)

            cat("\nStrata:")
            print(table(object@strata))

            if (!is.null(object@MIparam)) {
                cat("\nMultiple imputation in", object@MIparam$method)
                cat("\n---------------------------")
                cat("\nTotal number of possible imputations:",
                    object@MIparam$Mtotal)
                cat("\nNumber of multiple imputations:", object@MIparam$M)

                if (attr(object@MIparam$ncomp, "estimated")) {
                cat("\nEstimated number of components for data imputation:",
                    object@MIparam$ncomp)
                } else {
                cat("\nNo estimated number of components for data imputation.")
                cat(" Defaults", object@MIparam$ncomp)
                }
            }
        })

##- incompleteData
setMethod(f="incompleteData", signature="MIDTList",
        definition=function(object) object@incompleteData)

##- strata
setMethod(f="strata", signature="MIDTList",
        definition=function(object) object@strata)

##- tableNames
setMethod(f="tableNames", signature="MIDTList",
        definition=function(object, ...) object@tableNames)

setReplaceMethod(f="tableNames", signature="MIDTList",
                definition=function(object, value) {
                nt <- length(object@incompleteData)

                if (length(value) != nt | is.matrix(value) | 
                    is.list(value)) {
                    stop("'tableNames<-' accessor is only valid for vectors",
                        " of length ", nt, call.=FALSE)
                }

                value <- as.character(value)

                if (any(duplicated(value))) {
                    stop("'tableNames<-' accessor is only valid for vectors",
                        " with unique values", call.=FALSE)
                }

                object@tableNames  <- value
                names(object@incompleteData) <- value
                return(object)
            })

##- missingRows
setMethod(f="missingRows", signature="MIDTList",
        definition=function(object) object@missingRows)

##- compromise
setMethod(f="compromise", signature="MIDTList",
        definition=function(object) {
            if (is.null(object@compromise)) {
                cat("No 'compromise' slot found in the MIDTList object.",
                    "Run MI first.")
            } else {
                object@compromise
            }
        })

##- configurations
setMethod(f="configurations", signature="MIDTList",
        definition=function(object, M="all") {
            if (is.null(object@configurations)) {
                cat("No 'configurations' slot found in the MIDTList object.",
                    "Run MI first.")
            } else {
                if (M == "all") {
                    object@configurations
                } else {
                    object@configurations[[M]]
                }
            }
        })

##- imputedRows
setMethod(f="imputedRows", signature="MIDTList",
        definition=function(object) {
            if (is.null(object@imputedRows)) {
                cat("No 'imputedRows' slot found in the MIDTList object.",
                    "Run MI first.")
            } else {
                object@imputedRows
            }
        })

##- MIparam
setMethod(f="MIparam", signature="MIDTList",
        definition=function(object) {
            if (is.null(object@MIparam)) {
                cat("No 'MIparam' slot found in the MIDTList object.",
                    "Run MI first.")
            } else {
                object@MIparam
            }
        })

