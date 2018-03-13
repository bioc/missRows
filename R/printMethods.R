
##- print method for missPattern ----------------------------------------------#
##-----------------------------------------------------------------------------#
print.missPattern <- function(x, ...) {
    
    cat("Number of missing/available rows in ")
    cat("each stratum per data table.\n")
    cat("Table:\n")
    print(x$nbMissing)
}


##- print method for tuneM ----------------------------------------------------#
##-----------------------------------------------------------------------------#
print.tuneM <- function(x, ...) {
    
    print(x$stats)
}
