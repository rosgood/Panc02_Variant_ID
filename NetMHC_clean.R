#need to import files as list
files <- list.files(path = "~/data/Panc02_exome/tmp/", pattern = "*.reallyclean$")
for(i in files){
  path <- paste0("~/data/Panc02_exome/tmp/", i)
  predicts <- read.table(path, comment.char="", stringsAsFactors=FALSE, quote="\"") 
  #predicts <- read.table(path, quote="\"", comment.char="", stringsAsFactors=FALSE)
#need to get length of epitope for trimming
  predicts[,7] <- nchar(as.character(predicts$V2))
#sort by length then position
#need to correct for increasing peptide length from 19 to 20
  predicts.ordered <- predicts[order(predicts$V7, predicts$V1),]
#trim epitopes which do not contain mutation
  trims <- predicts.ordered[!(predicts.ordered$V1<3 & predicts.ordered$V7 == 8),]
  trims <- trims[!(trims$V1 > 10 & trims$V7 == 8),]
  trims <- trims[!(trims$V1 < 2 & trims$V7 == 9),]
  trims <- trims[!(trims$V1 > 10 & trims$V7 == 9),]
  trims <- trims[!(trims$V1 < 1 & trims$V7 ==10),]
  predicts.best <- trims[trims$V4 < 1000,]
  write.table(predicts.best, file = paste0("~/data/Panc02_exome/tmp/", sub(".reallyclean", "_best.predicted", i)), sep = "\t", col.names = FALSE, row.names = FALSE)
}

