##- required libraries -------------------------------------------------------#
##----------------------------------------------------------------------------#
library(S4Vectors)
library(MultiAssayExperiment)
library(omicade4)
library(rcellminer)
library(rcellminerData)


##- transcriptomic data ------------------------------------------------------#
##----------------------------------------------------------------------------#
##- load data
data(NCI60_4arrays)

##- a subset of microarray gene expression from the NCI-60 cell lines
trans <- NCI60_4arrays$agilent

##- incomplete trascriptomic data
trans <- trans[, -c(1, 6, 12, 13, 19, 25, 26, 35, 36, 44, 53, 55)]

##-- map
colname <- colnames(trans)
primary <- apply(as.matrix(colname), 1, 
                 function(x) { strsplit(x, "[.]")[[1]][2] })
transmap <- DataFrame(primary = primary, colname = colname)


##- proteomic data -----------------------------------------------------------#
##----------------------------------------------------------------------------#
##- load data
data(molData)

##-- proteomic data from the NCI 60 cell lines
prote <- getESetList(molData)$pro
colData <- getSampleData(molData)

##-- incomplete proteomic data
prote <- prote[, -c(8, 20, 27, 28, 39, 40, 45, 54)]

##-- map
colname <- sampleNames(prote)
primary <- apply(as.matrix(colname), 1, 
                 function(x) { strsplit(x, "[:]")[[1]][2] })
protemap <- DataFrame(primary = primary, colname = colname)


##- cell line information ----------------------------------------------------#
##----------------------------------------------------------------------------#
tmp <- getSampleData(molData)$Name
tmp <- data.frame(tmp, stringsAsFactors = FALSE)
cell.type <- apply(tmp, 1, function(x) { strsplit(x, "[:]")[[1]][1] })
samples <- apply(tmp, 1, function(x) { strsplit(x, "[:]")[[1]][2] })
cell.line <- DataFrame(type = cell.type, row.names = samples)


##- MultiAssayExperiment for the incomplete data -----------------------------#
##----------------------------------------------------------------------------#
listmap <- list("trans" = transmap, "prote" = protemap)
dfmap <- listToMap(listmap)

explist <- list("trans" = trans, "prote" = prote)

mae <- MultiAssayExperiment(experiments=explist,
                            colData=cell.line,
                            sampleMap=dfmap)


##- liste containing both data tables and cell line info ---------------------#
##----------------------------------------------------------------------------#
transTable <- assays(mae)$trans
colnames(transTable) <- transmap$primary

proteTable <- assays(mae)$prote
colnames(proteTable) <- protemap$primary

dataTables <- list(trans = transTable, prote = proteTable, 
                    cell.line = cell.line)


##- the NCI60 data for missRows ----------------------------------------------#
##----------------------------------------------------------------------------#
NCI60 <- list("dataTables" = dataTables, "mae" = mae)
