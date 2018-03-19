context("test class constructors")

data(NCI60)

test_that("MIDTList object from a list of data tables", {
    object <- newMIDTList(NCI60$incompleteData, strata=NCI60$cell.line)
    expect_is(object, "MIDTList")
})

test_that("MIDTList object from separate data tables", {
    table1 <- NCI60$incompleteData$trans
    table2 <- NCI60$incompleteData$prote
    
    object <- newMIDTList(table1, table2, strata=NCI60$cell.line)
    expect_is(object, "MIDTList")
})


context("tests on inputs")

test_that("tests for data tables", {
    table1 <- NCI60$incompleteData$trans
    table2 <- NCI60$incompleteData$prote
    
    expect_error(newMIDTList(table1, strata=NCI60$cell.line), 
                "at least two data tables must be passed as arguments in '...'")
    
    expect_error(newMIDTList(table1, table2[-1, ], strata=NCI60$cell.line), 
                "non equal row numbers among tables and/or strata")
    
    tableNrn <- table1
    rownames(tableNrn) <- NULL
    expect_error(newMIDTList(tableNrn, table2, strata=NCI60$cell.line), 
                "the 'Table 1' data table must be rows named.")
    
    expect_error(newMIDTList(1:nrow(table1), table2, strata=NCI60$cell.line), 
                "the 'Table 1' data table must be a matrix or data frame.")
    
    tableInf <- table1
    tableInf[1, 1] <- Inf
    expect_error(newMIDTList(tableInf, table2, strata=NCI60$cell.line), 
                 "infinite values in 'Table 1' data table.")
})
