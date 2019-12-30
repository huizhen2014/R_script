##ChIPQC
##1. Experiment sample sheet
library(ChIPQC)
samples <- read.csv(file.path(system.file("extdata",package = "ChIPQC"),
                              "example_QCexperiment.csv"))

##2. Constructing a ChIPQCexperiment object
exampleExp <- ChIPQC(samples,annotation="hg19")

##for example
data("example_QCexperiment")
exampleExp

##4. Generationg a summary QC report for experimental sample groups
ChIPQCreport(exampleExp)
















##ChIPpeak Anno

library(ChIPpeakAnno)
path <- system.file("extdata","Tead4.broadPeak",package="ChIPpeakAnno")
peaks <- toGRanges(path, format="broadPeak")
