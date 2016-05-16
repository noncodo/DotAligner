#install.packages("ape", "ggplot2", "RColorBrewer")
library(ape)
library(ggplot2) 
library(RColorBrewer)
library(gplots)
library(grid)
library(pROC)

# args <- commandArgs(trailingOnly = TRUE)
# path <- args[1]
# file.names <- dir(path, pattern =".matrix.low")

#for(i in 1:length(file.names)){ }
# b <- as.matrix( read.table("binary.matrix", header=T))
ids <- read.delim("structure_ids.txt_low", header=F)
lDA <- read.delim("dotaligner_scores.matrix_low", header=F)
lL <- read.delim("locarna_scores.matrix_low", header=F)
lC <- read.delim("carna_scores.matrix_low", header=F)

hids <- read.delim("structure_ids.txt_high", header=F)
hDA <- read.delim("dotaligner_scores.matrix_high", header=F)
hL <- read.delim("locarna_scores.matrix_high", header=F)
hC <- read.delim("carna_scores.matrix_high", header=F)

# make a binary similarity matrix
b <- matrix(ncol=nrow(ids), nrow=nrow(ids))
for (i in 1:nrow(ids)) {
    for (j in 1:nrow(ids)) {
      ifelse( substr(ids[i,],1,7) == substr(ids[j,],1,7) , b[i,j] <- 1 , b[i,j] <- 0 ) 
  	}
}
rownames(b) <- ids[,]
colnames(b) <- ids[,]

hb <- matrix(ncol=nrow(hids), nrow=nrow(hids))
for (i in 1:nrow(hids)) {
    for (j in 1:nrow(hids)) {
      ifelse( substr(hids[i,],1,7) == substr(hids[j,],1,7) , hb[i,j] <- 1 , hb[i,j] <- 0 ) 
  	}
}
rownames(hb) <- hids[,]
colnames(hb) <- hids[,]

colnames(lDA) <- c('x','y','file1','file2','score')
lmDA <- matrix(nrow= max(lDA$y), ncol=max(lDA$y))
lmDA[cbind(lDA$x,lDA$y)] <- as.numeric(lDA$score)
lmDA[cbind(lDA$y,lDA$x)] <- as.numeric(lDA$score)
rownames(lmDA) <- ids[,]
colnames(lmDA) <- ids[,]

colnames(lC) <- c('x','y','file1','file2','score')
lmC <- matrix(nrow= max(lC$y), ncol=max(lC$y))
lmC[cbind(lC$x,lC$y)] <- as.numeric(lC$score)
lmC[cbind(lC$y,lC$x)] <- as.numeric(lC$score)
rownames(lmC) <- ids[,]
colnames(lmC) <- ids[,]

colnames(lL) <- c('x','y','file1','file2','score')
lmL <- matrix(nrow= max(lL$y), ncol=max(lL$y))
lmL[cbind(lL$x,lL$y)] <- as.numeric(lL$score)
lmL[cbind(lL$y,lL$x)] <- as.numeric(lL$score)
rownames(lmL) <- ids[,]
colnames(lmL) <- ids[,]

###
colnames(hDA) <- c('x','y','file1','file2','score')
hmDA <- matrix(nrow= max(hDA$y), ncol=max(hDA$y))
hmDA[cbind(hDA$x,hDA$y)] <- as.numeric(hDA$score)
hmDA[cbind(hDA$y,hDA$x)] <- as.numeric(hDA$score)
rownames(hmDA) <- hids[,]
colnames(hmDA) <- hids[,]

colnames(hC) <- c('x','y','file1','file2','score')
hmC <- matrix(nrow= max(hC$y), ncol=max(hC$y))
hmC[cbind(lC$x,lC$y)] <- as.numeric(hC$score)
hmC[cbind(lC$y,lC$x)] <- as.numeric(hC$score)
rownames(hmC) <- hids[,]
colnames(hmC) <- hids[,]

colnames(hL) <- c('x','y','file1','file2','score')
hmL <- matrix(nrow= max(hL$y), ncol=max(hL$y))
hmL[cbind(hL$x,hL$y)] <- as.numeric(hL$score)
hmL[cbind(hL$y,hL$x)] <- as.numeric(hL$score)
rownames(hmL) <- hids[,]
colnames(hmL) <- hids[,]

#Pearson correlation
lDA.pc <- cor.test(lmDA,b)
lC.pc <- cor.test(lmC,b)
lL.pc <- cor.test(lmL,b) 
# pc$p.value pc$estimate
hDA.pc <- cor.test(hmDA,b)
hC.pc <- cor.test(hmC,b)
hL.pc <- cor.test(hmL,b) 

#Mantel test
lDA.mt <- mantel.test(lmDA,b,nperm=ncol(lmDA)*100) 
lC.mt <- mantel.test(lmC,b,nperm=ncol(lmC)*100)
lL.mt <- mantel.test(lmL,b,nperm=ncol(lmL)*100) 

hDA.mt <- mantel.test(hmDA,b,nperm=ncol(hmDA)*100) 
hC.mt <- mantel.test(hmC,b,nperm=ncol(hmC)*100)
hL.mt <- mantel.test(hmL,b,nperm=ncol(hmL)*100) 
# mt$z.stat, mt$p

###################
# 	PLOTTING 
####################
#get the RFAM ID of all sampled sequences
fam <- substr( ids[,],1,7 )
# get unique families
uniqueFam <- unique(fam)
# get a color palette for this number
#colors <- brewer.pal(length(uniqueFam) ,"Set3")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#merge 
color.map <- cbind( uniqueFam, col_vector[12:(length(uniqueFam)+12-1)] )
fam.colors <- merge(fam, color.map, by=1)

hfam <- substr( hids[,],1,7 )
huniqueFam <- unique(hfam)
hcolor.map <- cbind( huniqueFam, col_vector[12:(length(huniqueFam)+12-1)] )
hfam.colors <- merge(hfam, hcolor.map, by=1)

#titles 
lDAt <- cbind(lDA.pc$method, lDA.pc$estimate)
lLt <- cbind(lL.pc$method, lL.pc$estimate)
lCt <- cbind(lC.pc$method, lC.pc$estimate)

hDAt <- cbind(hDA.pc$method, hDA.pc$estimate)
hLt <- cbind(hL.pc$method, hL.pc$estimate)
hCt <- cbind(hC.pc$method, hC.pc$estimate)

#heatmap.2(mDAo, Rowv=NA, Colv=NA , symm=T, trace='none', col=bluered, RowSideColors=as.character(fam.colors$colors))
#ggplot(DA, aes(file1,file2))+geom_tile(aes(fill=score))+geom_tile(aes(file2,file1,fill=score)) + scale_fill_gradient2()
#heatmap(b, Colv=NA, RowSideColors=as.character(fam.colors$colors), ColSideColors=as.character(fam.colors$colors))
pdf("low_pid_unsorted.pdf")
heatmap.2(b,  Rowv=NA, main="Stochastic RFAM sampling, low PID [0-55]", labCol=fam, symm=T, 
				RowSideColors=as.character(fam.colors$V2), ColSideColors=as.character(fam.colors$V2), trace='none')
legend('bottomleft', legend = data.frame(color.map)$uniqueFam, fill = as.character(data.frame(color.map)$V2), ncol=1 , cex=0.8 )

heatmap.2(lmDA, Colv=NA, Rowv=NA, main=lDAt, labCol=fam, labRow=fam, symm=T,RowSideColors=as.character(fam.colors$V2), ColSideColors=as.character(fam.colors$V2), trace='none')
heatmap.2(lmL, Colv=NA, Rowv=NA, main=lLt, labCol=fam, labRow=fam, symm=T,RowSideColors=as.character(fam.colors$V2), ColSideColors=as.character(fam.colors$V2), trace='none')
heatmap.2(lmC, Colv=NA, Rowv=NA, main=lCt, labCol=fam, labRow=fam, symm=T,RowSideColors=as.character(fam.colors$V2), ColSideColors=as.character(fam.colors$V2), trace='none')
dev.off()

pdf("h_b_heat.pdf")
heatmap.2(hb,  Rowv=NA, main="Stochastic RFAM sampling, high PID [56-95]", labCol=hfam, symm=T, 
				RowSideColors=as.character(hfam.colors$V2), ColSideColors=as.character(hfam.colors$V2), trace='none')
legend('bottomleft', legend = data.frame(hcolor.map)$huniqueFam, fill = as.character(data.frame(hcolor.map)$V2), ncol=1 , cex=0.8 )
dev.off()
pdf("h_DA_heat.pdf")
heatmap.2(hmDA, Rowv=NA,  main=hDAt, labCol=hfam, labRow=hfam, RowSideColors=as.character(hfam.colors$V2), ColSideColors=as.character(hfam.colors$V2), trace='none')
dev.off()
pdf("h_L_heat.pdf")
heatmap.2(hmL, Colv=NA, Rowv=NA, main=hLt, labCol=hfam, labRow=hfam, symm=T,RowSideColors=as.character(hfam.colors$V2), ColSideColors=as.character(hfam.colors$V2), trace='none')
dev.off()
pdf("h_C_heat.pdf")
heatmap.2(hmC, Colv=NA, Rowv=NA, main=hCt, labCol=hfam, labRow=hfam, symm=T,RowSideColors=as.character(hfam.colors$V2), ColSideColors=as.character(hfam.colors$V2), trace='none')
dev.off()

##### PEARSONS CORRELATION PLOT
names <- cbind( c("DotAligner","Locarna","Carna"))
df <- data.frame(rbind( lDA.pc$estimate, lL.pc$estimate, lC.pc$estimate ))
conf95 <- data.frame(rbind( lDA.pc$conf.int[1:2],  lL.pc$conf.int[1:2], lC.pc$conf.int[1:2]) )
lowPC <- cbind( names, df, conf95)
lowPC$pid <- "low_PID"
hdf <- data.frame(rbind( hDA.pc$estimate, hL.pc$estimate, hC.pc$estimate ))
hconf95 <- data.frame(rbind( hDA.pc$conf.int[1:2] , hL.pc$conf.int[1:2], hC.pc$conf.int[1:2]))
highPC <- cbind( names, hdf, hconf95)
highPC$pid <- "high_PID"
PC <- rbind( lowPC, highPC)

dodge <- position_dodge(width=0.9)
p <- ggplot( PC , aes(names, cor, fill=pid))
p + geom_bar(position=dodge,stat="identity") + 
	geom_errorbar(aes(ymin = X1, ymax=X2), position=dodge, width=0.25) +
	theme_bw() +
	ylim(0,1) + 
	ggtitle("Pearson's product-moment correlation with\nbinary RFAM classification matrix") +
	ylab("Pearson's Correlation\n[two-tailed, 95%c.i.]") +
	xlab("Algorithm") +
	scale_fill_brewer(palette="Set1")
	#+	theme(legend.position="bottom")

##### PLOTTING RUNTIME
time.l <- read.delim("locarna_low_time.txt", header=F)
time.c <- read.delim("carna_low_time.txt", header=F)
time.d <- read.delim("dotaligner_low_time.txt", header=F)
time.d$Algorithm <- as.character("DotAligner")
time.c$Algorithm <- as.character("CARNA")
time.l$Algorithm <- as.character("LocaRNA")
times <- rbind( time.l , time.c, time.d )
ggplot(times, aes( factor(Algorithm), log(V1*1000)/log(10))) + geom_violin(aes(fill=Algorithm)) +geom_boxplot(width=0.1, fill="white") + theme_bw()



# ROC 
lroc <- NULL
for (i in 1:nrow(ids)) {
    for (j in i:nrow(ids)) {
    	ifelse( substr(ids[i,],1,7) == substr(ids[j,],1,7) , lroc <- rbind(lroc, 1) , lroc <- rbind(lroc, 0)) 
  	}
}
plot(roc( as.numeric(lroc), lC$score) , main=auc(as.numeric(lroc), lC$score))
plot(roc( as.numeric(lroc), lL$score) , main=auc(as.numeric(lroc), lL$score))
plot(roc( as.numeric(lroc), lDA$score) , main=auc(as.numeric(lroc), lDA$score))




multiplot()
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}