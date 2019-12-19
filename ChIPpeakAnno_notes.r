##

##ChIPpeak Anno

library(ChIPpeakAnno)
path <- system.file("extdata","Tead4.broadPeak",package="ChIPpeakAnno")
peaks <- toGRanges(path, format="broadPeak")
