context("test class constructors")

data(NCI60)

test_that("MIDTList object from a list of data tables", {
    tableList <- NCI60$dataTables[1:2]
    colData <- NCI60$dataTables$cell.line
    
    midt <- MIDTList(tableList, colData=colData, 
                        assayNames=c("trans", "prote"))
    expect_is(midt, "MIDTList")
})

test_that("MIDTList object from separate data tables", {
    table1 <- NCI60$dataTables$trans
    table2 <- NCI60$dataTables$prote
    colData <- NCI60$dataTables$cell.line
    
    midt <- MIDTList(table1, table2, colData=colData,
                        assayNames=c("transcrip", "proteome"))
    expect_is(midt, "MIDTList")
})

test_that("MIDTList object from 'MultiAssayExperiment'", {
    midt <- MIDTList(NCI60$mae)
    expect_is(midt, "MIDTList")
})


context("tests on inputs")

test_that("tests for data tables", {
    table1 <- NCI60$dataTables$trans
    table2 <- NCI60$dataTables$prote
    colData <- NCI60$dataTables$cell.line
    
    expect_error(MIDTList(table1, colData=colData), 
            "at least two data tables must be passed as arguments in '...'")
    
    tableNrn <- table1
    colnames(tableNrn) <- NULL
    expect_error(MIDTList(tableNrn, table2, colData=colData), 
                "the 'Table 1' data table must be columns named.")
    
    expect_error(MIDTList(1:nrow(table1), table2, colData=colData), 
                "the 'Table 1' data table must be a matrix or data frame.")
    
    tableInf <- table1
    tableInf[1, 1] <- Inf
    expect_error(MIDTList(tableInf, table2, colData=colData), 
                 "infinite values in 'Table 1' data table.")
})
