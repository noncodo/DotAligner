# Rscript auc.R bench/

suppressMessages(library(pROC))
args <- commandArgs(trailingOnly = TRUE)

cat("Generating classification matrix...")
ids <- read.delim("structure_ids.txt", header = F)
mids <- NULL
for (i in 1:nrow(ids)) {
    for (j in i:nrow(ids)) {
    	ifelse( substr(ids[i,],1,7) == substr(ids[j,],1,7) , mids <- rbind(mids, 1) , mids <- rbind(mids, 0)) 
  	}
}
cat("Done\n")
cat("Importing benchmarking data and performing ROC...")
path <- args[1]
file.names <- dir(path, pattern =".short")
out <- NULL
for(i in 1:length(file.names)){
	m <- read.delim( paste0(path,file.names[i]) , header =F, sep =" ")
	#print( summary(m) )
	a <- auc( as.numeric(mids), m$V3)
	rbind( out, paste(file.names[i],a) )
}
cat("Done\nWriting to file...")
write.table(out, file="auc.tsv", quote=F, row.names=F, col.names=F)
cat("Done\nExiting R\n")