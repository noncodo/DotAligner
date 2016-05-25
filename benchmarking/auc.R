# Rscript auc.R bench/

suppressMessages(library(pROC))
args <- commandArgs(trailingOnly = TRUE)

cat("[R] Generating classification matrix...")
ids <- read.delim("structure_ids.txt", header = F)
mids <- NULL
for (i in 1:nrow(ids)) {
    for (j in i:nrow(ids)) {
    	ifelse( substr(ids[i,],1,7) == substr(ids[j,],1,7) , mids <- rbind(mids, 1) , mids <- rbind(mids, 0)) 
  	}
}
cat("Done\n")
cat("[R] Importing benchmarking data and performing ROC...\n")
cat(" Percentage complete: 0%")
path <- args[1]
file.names <- dir(path, pattern =".short")
out <- NULL
 for(i in 1:length(file.names)) {
 	m <- read.delim( paste0(path,file.names[i]) , header =F, sep =" ")
 	# ran this in bash a priori 
 	# sed -i.old '/egm/d' ./bench/*.time
 	t <- read.delim( paste0(path, gsub(".score.short",".time",file.names[i]) ), header =F, sep =" ")
 	a <- auc( as.numeric(mids), m$V3)
 	row <- cbind( file.names[i], a , mean(t$V1), sd(t$V1))
	out <- rbind( out, row )
	# print progress counter
 	if ( i %% as.integer( length(file.names)/20 ) == 0 ) {
 		pc <-  as.integer(i/(length(file.names)) * 100)+1
 		cat( paste0("..",pc,"%"))
 	}
}
cat("\n[R] Writing to file...")
write.table(out, file="auc_time.tsv", quote=F, row.names=F, col.names=F)
cat("[R] Done\nExiting R\n")
#plot(roc1) 
#plot(roc2, add=TRUE, col='red')


## in bash 
# cut -d " " -f 1 auc_time.tsv | sed -e 's/T-//g' -e 's/_s-/       /g' -e 's/_k-/  /g' -e 's/_t-/  /g' -e 's/_o-/       /g' -e 's/_e-/  /g' -e 's/\.dot.*//g' > auc_time_params.tsv
#( echo "file     auc     time    sd      T       s       k       t       o       e"; paste auc_time.tsv auc_time_params.tsv | sed 's/ /  /g' ) > benchmark_results.tsv