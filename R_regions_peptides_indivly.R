library("Biostrings", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library("seqinr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library("IRanges", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
mut_locs <- list.files(path = "~/data/Panc02_exome/tmp/", pattern = "*_mut_locs_w_ID")
for(i in mut_locs){
  path <- paste0("~/data/Panc02_exome/tmp/", i)
  pt_mut <- read.table(path, sep='\t', quote="\"", comment.char="", stringsAsFactors=FALSE)
  pt_mut$pep_start <- pt_mut[,2]-10
  pt_mut$pep_stop <- pt_mut[,2]+9
  write.table(file=paste0("~/data/Panc02_exome/tmp/", sub("_mut_locs_w_ID", "_regs", i)), pt_mut, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}
regs <- read.csv(file = "/home/heatherkinkead/data/Panc02_exome/tmp/bt2_Panc02_regs", header = TRUE, sep = "\t")
f <- readAAStringSet("/home/heatherkinkead/data/Panc02_exome/tmp/bt2_Panc02_best_ID_R.fasta", format = "fasta", use.names = TRUE, skip = 0L)

#get peptides
getAAseqs <- function(x) {
  idx <- which(names(f) == x[1])
  o <- character(length(idx))
  n <- character(length(idx))
  for(a in seq(length(idx))){ # this is awful... who wrote this shit?
    i <- which(names(f) == x[a,1])
    l <- length(f[[i]])
    s <- x[a,3]
    e <- x[a,4]
    if(s<0){s <- 0}
    if(e>l){e <- l}
    if(s>e){
	o[a] <- NA
	warning(paste('Something is messed up... Check', as.character(x[a,1])))
    }else{o[a] <- as.character(f[[i]][s:e])}
    n[a] <- as.character(x[a,1])
  }

  ds <- data.frame(name=n,peptide=o)

}

ds <- getAAseqs(regs)

write.table(ds, file = "~/data/Panc02_exome/tmp/bt2_Panc02_peptides.csv", 
	quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

write.fasta(sequences = as.list(ds$peptide), 
names = ds$name, file = "~/data/Panc02_exome/tmp/bt2_Panc02_peptides.fasta")
