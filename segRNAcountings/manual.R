library(roxygen2)
roxygenise("/home/danilo/Documents/Stage_sources/segRNAcountings")

pack <- "segRNAcountings"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))
