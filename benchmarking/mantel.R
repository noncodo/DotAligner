install.packages("ape")
library(ape)
DA <- read.delim("DA.log", header=F)
colnames(DA) <- c('x','y','score')
mDA <- matrix(nrow= max(DA$y), ncol=max(DA$y))
mDA[cbind(DA$x,DA$y)] <- as.numeric(DA$score)
mDA[cbind(DA$y,DA$x)] <- as.numeric(DA$score)
mDA[cbind(DA$y,DA$y)] <- 1
mDA[cbind(1,1)] <- 1
mt <- mantel.test(mDA,b,nperm=ncol(mDA)*100)
# mt$z.stat, mt$p
pc <- cor.test(mDA,b)
# pc$p.value pc$estimate