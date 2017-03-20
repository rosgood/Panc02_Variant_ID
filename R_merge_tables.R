files <- list.files(path = "~/data/Panc02_exome/tmp/", pattern = "*_muts_prots_only_fixed")
#this code seems to duplicate the bash code where I create a unique ID
for(i in files){
  path <- paste0("~/data/Panc02_exome/tmp/", i)
  tables <- read.table(path, quote="\"", comment.char="", stringsAsFactors=FALSE)
  tables$unique_ID <- paste(tables$V1, tables$V2)
  tables$mutation <- paste(tables$V4, tables$V5, tables$V6, tables$V7, tables$V8, tables$V9, tables$V10)
  tables[,"unique_ID"] <- tables$unique_ID
  tables[,"mut_type"] <- tables$V3
  tables[,"mutation"] <- tables$mutation
  tables[,"sequence"] <- tables$V11
  tables <- subset(tables, select = -c(V1:V2, V3, V4:V10, V11))
  write.table(file=paste0("~/data/Panc02_exome/tmp/", sub("_muts_prots_only_fixed", "_muts_prots_R", i)), tables, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}
#get UniqueID for matched names and prot mut
refs <- list.files(path = "~/data/Panc02_exome/tmp/", pattern = "*_matched_names_ref_gene")
for(i in refs){
  path2 <- paste0("~/data/Panc02_exome/tmp/", i)
  names <- read.delim(path2, header=FALSE, stringsAsFactors=FALSE)
  names$unique_ID <- paste(names$V2, names$V4)
  names[, "unique_ID"] <- names$unique_ID
  names <- subset(names, select = -c(V2,V4))
  names <- names[c(4,1,3)]
  
  write.table(file = paste0("~/data/Panc02_exome/tmp/", sub("_matched_names_ref_gene", "_ref_genes_R", i)), names, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}

#merge tables to get names and seqs together
##start here
files <- list.files(path = "~/data/Panc02_exome/tmp/", pattern = "*_muts_prots_R")
refs <- list.files(path = "~/data/Panc02_exome/tmp/", pattern = "*_ref_genes_R.txt")
for(i in seq(length(files))){
  seqs <- read.delim(paste0("~/data/Panc02_exome/tmp/", files[i]))
  IDs <- read.delim(paste0("~/data/Panc02_exome/tmp/", refs[i]))
  merged <- merge(IDs, seqs, by = "unique_ID", all.y = TRUE, sort = TRUE)
  
  write.table(merged, file = paste0("~/data/Panc02_exome/tmp/", sub("_muts_prots_R", "_merged_muts_names_R", files[i])), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
#get better IDs
files3 <- list.files(path = "~/data/Panc02_exome/tmp/", pattern = "*_merged_muts_names_R$")
for(i in files3){
  path3 <- paste0("~/data/Panc02_exome/tmp/", i)
  names2 <- read.delim(path3, header=TRUE, stringsAsFactors=FALSE)
  names2$prot_ID <- paste0(names2[,2],":",names2[,3])
  names2[, "prot_ID"] <- names2$prot_ID
  names2 <- subset(names2, select = -c(V1))
  names2 <- names2[c(6,1,2,5)]
  #write.table(names2, file = paste0("~/data/Panc02_exome/tmp/", sub("_merged_muts_names_R", "_best_ID_R", files[i])), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(file=paste0("~/data/Panc02_exome/tmp/", sub("_merged_muts_names_R", "_best_ID_R", i)), names2, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  }