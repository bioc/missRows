#-- transcriptomic data --------------------------------------------------------#
#-------------------------------------------------------------------------------#

#-- load data
library(omicade4)
data(NCI60_4arrays)

#-- a subset of microarray gene expression of the NCI 60 cell lines
trans <- t(NCI60_4arrays$agilent)


#-- proteomic data -------------------------------------------------------------#
#-------------------------------------------------------------------------------#

#-- load data
library(rcellminer)
library(rcellminerData)
data(molData)

#-- proteomic data of the NCI 60 cell lines
prote <- t(exprs(molData@eSetList$pro))


#-- cell line information ------------------------------------------------------#
#-------------------------------------------------------------------------------#
samples <- getSampleData(molData)$Name
strata <- data.frame(samples = samples, stringsAsFactors = FALSE)
rownames(strata) <- samples

for (i in 1:length(samples)) {
  strata[i, ] <- strsplit(samples[i], "[:]")[[1]][1]
  rownames(trans)[i] <- paste(strsplit(rownames(trans)[i], "[.]")[[1]],
                              collapse = ":")
}


#-- complete data --------------------------------------------------------------#
#-------------------------------------------------------------------------------#

#-- liste containing the data tables
completeData <- list(trans = trans, prote = prote)

#-- check whether samples are ordered correctly
all(apply(x <- sapply(completeData, rownames), 2,
          function(y) identical(y, x[, 1])))


#-- incomplete data ------------------------------------------------------------#
#-------------------------------------------------------------------------------#

#-- incomplete trascriptomic data
trans.incompl <- trans
trans.incompl[c(1, 6, 12, 13, 19, 25, 26, 35, 36, 44, 53, 55), ] <- NA

#-- incomplete proteomic data
prote.incompl <- prote
prote.incompl[c(8, 20, 27, 28, 39, 40, 45, 54), ] <- NA

#-- liste containing the data tables
incompleteData <- list(trans = trans.incompl, prote = prote.incompl)


#-- liste containing both data sets and cell line info -------------------------#
#-------------------------------------------------------------------------------#
NCI60 <- list(completeData = completeData, incompleteData = incompleteData,
              cell.line = strata)
