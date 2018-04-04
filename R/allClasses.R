###############################################################################
### MIDTList S4 class definition
###############################################################################

##- setClass -----------------------------------------------------------------#
##----------------------------------------------------------------------------#
MIDTList <- setClass("MIDTList",
                    contains="MultiAssayExperiment",
                    slots=c(
                    strata="ANY",
                    missingIndv="ANY",
                    compromise="ANY",
                    configurations="ANY",
                    imputedIndv="ANY",
                    MIparam="ANY"))


##- setGeneric ---------------------------------------------------------------#
##----------------------------------------------------------------------------#

##- strata
setGeneric(name="strata",
            def=function(object) standardGeneric("strata"))

##- missingIndv
setGeneric(name="missingIndv",
            def=function(object) standardGeneric("missingIndv"))

##- compromise
setGeneric(name="compromise",
            def=function(object) standardGeneric("compromise"))

##- configurations
setGeneric(name="configurations",
            def=function(object, ...) standardGeneric("configurations"))

##- imputedIndv
setGeneric(name="imputedIndv",
            def=function(object) standardGeneric("imputedIndv"))

##- MIparam
setGeneric(name="MIparam",
            def=function(object) standardGeneric("MIparam"))


##- setMethod ----------------------------------------------------------------#
##----------------------------------------------------------------------------#
##- show
setMethod("show",
        signature="MIDTList",
        definition=function(object) {

            nbMiss <- unlist(lapply(object@missingIndv, length))
            nt <- length(assays(object))

            cat("An object of class ", class(object), ".",
                "\n\nTables:\n", sep = "")
            info <- data.frame(
                            vapply(assays(object), nrow, 1L),
                            vapply(assays(object), ncol, 1L),
                            nbMiss,
                            row.names = names(assays(object))
                                )
            colnames(info) <- c("features", "individuals", "miss.indv")
            print(info)

            cat("\nStrata:")
            print(table(colData(object)[, object@strata]))

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

##- missingIndv
setMethod(f="missingIndv", signature="MIDTList",
        definition=function(object) object@missingIndv)

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

##- imputedIndv
setMethod(f="imputedIndv", signature="MIDTList",
        definition=function(object) {
            if (is.null(object@imputedIndv)) {
                cat("No 'imputedIndv' slot found in the MIDTList object.",
                    "Run MI first.")
            } else {
                object@imputedIndv
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
