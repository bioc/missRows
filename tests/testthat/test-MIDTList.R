context("Checking class constructors")

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