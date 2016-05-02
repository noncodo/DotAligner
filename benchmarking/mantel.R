#install.packages("ape", "ggplot2")
library(ape)
library(ggplot2) 

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
file.names <- dir(path, pattern =".matrix")

#for(i in 1:length(file.names)){ }
# b <- as.matrix( read.table("binary.matrix", header=T))
ids <- read.delim("structure_ids.txt", header=F)
DA <- read.delim("dotaligner_scores.matrix", header=F)
L <- read.delim("locarna_scores.matrix", header=F)
C <- read.delim("carna_scores.matrix", header=F)

# make a binary similarity matrix
b <- matrix(ncol=nrow(ids), nrow=nrow(ids))
for (i in 1:nrow(ids)) {
    for (j in 1:nrow(ids)) {
      ifelse( substr(ids[i,],1,7) == substr(ids[j,],1,7) , b[i,j] <- 1 , b[i,j] <- 0 ) 
  	}
}
rownames(b) <- ids[,]
colnames(b) <- ids[,]

colnames(DA) <- c('x','y','file1','file2','score')
mDA <- matrix(nrow= max(DA$y), ncol=max(DA$y))
mDA[cbind(DA$x,DA$y)] <- as.numeric(DA$score)
mDA[cbind(DA$y,DA$x)] <- as.numeric(DA$score)
rownames(mDA) <- ids[,]
colnames(mDA) <- ids[,]

colnames(C) <- c('x','y','file1','file2','score')
mC <- matrix(nrow= max(C$y), ncol=max(C$y))
mC[cbind(C$x,C$y)] <- as.numeric(C$score)
mC[cbind(C$y,C$x)] <- as.numeric(C$score)
rownames(mC) <- ids[,]
colnames(mC) <- ids[,]

colnames(L) <- c('x','y','file1','file2','score')
mL <- matrix(nrow= max(L$y), ncol=max(L$y))
mL[cbind(L$x,L$y)] <- as.numeric(L$score)
mL[cbind(L$y,L$x)] <- as.numeric(L$score)
rownames(mL) <- ids[,]
colnames(mL) <- ids[,]

#Matel test
DA.mt <- mantel.test(mDA,b,nperm=ncol(mDA)*100) 
C.mt <- mantel.test(mC,b,nperm=ncol(mC)*100)
L.mt <- mantel.test(mL,b,nperm=ncol(mL)*100) 
# mt$z.stat, mt$p

#Pearson correlation
DA.pc <- cor.test(mDA,b)
C.pc <- cor.test(mC,b)
L.pc <- cor.test(mL,b) 
# pc$p.value pc$estimate

#heatmap.2(b, Rowv="Colv" , symm=T, tracecol=NA)
#heatmap.2(b, Rowv=NA, Colv=NA , symm=T, trace='none', density.info='none', col=redblue)

# DA.pc 
# t = 169.43, df = 39998, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.6406532 0.6520641
# sample estimates:
#       cor
# 0.6463948

# C.pc
# t = 244.46, df = 39998, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.7700199 0.7778790
# sample estimates:
#       cor
# 0.7739793

# L.pc
# t = 181.4, df = 39998, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.6664137 0.6771676
# sample estimates:
#       cor
# 0.6718261

# > DA.mt
# $z.stat
# [1] 1311.425

# > L.mt
# $z.stat
# [1] 1264.826

# > C.mt
# $z.stat
# [1] 780.7338


