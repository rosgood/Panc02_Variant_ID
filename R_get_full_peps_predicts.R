files <- list.files(path = "~/data/Panc02_exome/tmp/", pattern = "*_uniq_predicts_w_header")
full_peps <- list.files(path = "~/data/Panc02_exome/tmp/", pattern = "*_Panc02_peptides.csv")
for(i in seq(length(files))){
  predicts <- read.delim(paste0("~/data/Panc02_exome/tmp/", files[i]))
  twentymers <- read.delim(paste0("~/data/Panc02_exome/tmp/", full_peps[i]))
  merged <- merge(predicts, twentymers, by = "name", all.x = TRUE, all.y=FALSE, sort = TRUE)

  write.table(merged, file = paste0("~/data/Panc02_exome/tmp/", sub("_uniq_predicts_w_header", "_uniq_predicts_w_full_peps", files[i])), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
