#!/bin/bash

#load binaries 
# TO DO: Make more user friendy (i.e. check for presence)
module load gi/ViennaRNA/2.1.3 marsmi/dotaligner/0.3 marsmi/R/3.2.3
DOTALNBIN=$( which DotAligner ) 
RNAFOLDBIN=$( which RNAfold ) 
TIMEFORMAT="%R"
BASEDIR="$(cd "$(dirname "${1}")" && pwd)/$(basename "${1}")"
BASEDIR=${BASEDIR%.*}

#################################################################
#################################################################
######
######	         CUSTOM PARAMETERS BELOW --- EDIT ME
######
#################################################################
#################################################################


MAXGRIDCPU=179 # make an option for non HPC environments
GRID_ACCOUNT="RNABiologyandPlasticity"


#################################################################
#	                        FUNCTIONS
#################################################################
>&2 echo -e "\e[91m             \`\`.:/osyyyys+:\`            "
>&2 echo -e "         \`.--:::::--::+oymNNms:         "
>&2 echo -e "       .-::-.-/sssyyyyyys++odNMm+\`      "
>&2 echo -e "     \`-::--odhsoosssyssooydmhhmNMm:     "
>&2 echo -e "    -:/-.smhyyhhhhhhhhhhhhhydNNNNMM+    "
>&2 echo -e "   :o+--mmhhhhddddddddddddhhhhmNNNMM+   "
>&2 echo -e "  .yy/.NmhdddddddddddddddddddddmNMNMN.  "
>&2 echo -e "  +dh-hNdddddddddd\e[39mBIG\e[91mdddddddddddmNMMMo  "
>&2 echo -e "  hNd/Nmdddddddddd\e[39mRED\e[91mdddddddddddhmMNMh  "
>&2 echo -e "  hNmyNmdddddddddd\e[39mBUTTON\e[91mdddddddhshMNNh  "
>&2 echo -e "  oMNmNNhdddddddddddddddddddddhh:mNNNo  "
>&2 echo -e "  .NMNNNdhhdddddddddddddddddhhho+MNNN.  "
>&2 echo -e "   +MMNMNdhhhdddddddddddddhhhh+/NNNN+   "
>&2 echo -e "    oMMNMNdhhhhhhdddddhhhhhho:yNNNN+    "
>&2 echo -e "     :mMNMMms+syhhhhhhhhyo//yNNNNm:     "
>&2 echo -e "      \`+mMMMMNds++////++ohNNNNNm+\`      "
>&2 echo -e "         :smMNNNNMMNNNNNNNNNds:         "
>&2 echo -e "            \`:+shddmmddhs+:\`            \e[39m"
>&2 echo -e                            

#################################################################
#                 Fold RNA to get dotplots
#################################################################
## TO DO : parallelize to make it faster 
if [[ ! -d $BASEDIR/logs ]]; then mkdir -p $BASEDIR/logs ; fi
#split fasta into smaller files
>&2 echo -ne "[*] Splitting fasta file ..."
if [[ ! -d ${BASEDIR}/split ]]; then mkdir -p $BASEDIR/split ; fi
cd $BASEDIR/split
NUMSEQ=$(wc -l ${BASEDIR}/../$1 | awk '{print $1}')
LINESperFILE=$(( ${NUMSEQ} / ${MAXGRIDCPU} ))
(( LINESperFILE+=1 ))
split -l $LINESperFILE ${BASEDIR%/*}/$1
  >&2 echo -e "[ \e[92mDONE\e[39m ]"  

#prepare folding command
  CMD="qsub -V -terse -sync y -cwd -P ${GRID_ACCOUNT} -N foldit 
 -o $BASEDIR/logs/foldit."$(date +"%d%m%y-%H%M")".out 
 -e $BASEDIR/logs/foldit."$(date +"%d%m%y-%H%M")".err 
 -pe smp 1 
 -t 1:${MAXGRIDCPU}
 -l mem_requested=1G,h_vmem=1G
 -b y -S /bin/bash
 ../foldit.sge"
cd $BASEDIR	&& >&2 echo -e "[*] working in BASEDIR: "$(pwd)
>&2 echo -ne "[*] Evaluating folder contents ..."
counter="-1"
if [[ ! -d ps ]]; then 
	mkdir ps
  	>&2 echo -e "no RNA structures found"
  	>&2 echo -e "[*] Generating pairwise probabilities..."
	$CMD | while read line ; do 
  		counter=$(( $counter + 1 ))
  		progress=$(( 100 * $counter / $MAXGRIDCPU ))
  		echo -ne "\r... "$progress"% complete" 
  	done
 	>&2 echo -e "\r[*] successfully completed command: "
 	>&2 echo -e "\e[31m"$CMD"\e[39m" 
elif [[ $( ls ps | wc -l ) -ne $( grep ">" ../${1} | wc -l ) ]]; then 
	>&2 echo -e "incomplete amount of structures"
	cd $BASEDIR	
	>&2 echo -e "[*] Generating pairwise probabilities..."
	$CMD | while read line ; do 
  		counter=$(( $counter + 1 ))
  		progress=$(( 100 * $counter / $MAXGRIDCPU ))
  		echo -ne "\r... "$progress"% complete" 
  	done
 	>&2 echo -e "\r[*] successfully completed command: "
 	>&2 echo -e "\e[31m"$CMD"\e[39m" 
else
	>&2 echo -e "found existing RNA structures "
  	>&2 echo -ne "    \e[31mDo you want to refold them?\e[39m [y/N] "
	read -n1 redo 
	>&2 echo 
	if [[ $redo != 'n' && $redo != 'N' ]]; then 
 		cd $BASEDIR	
 		>&2 echo -e "[*] Generating pairwise probabilities..."
		$CMD | while read line ; do 
  		counter=$(( $counter + 1 ))
  		progress=$(( 100 * $counter / $MAXGRIDCPU ))
  		echo -ne "\r... "$progress"% complete" 
  	done
 	>&2 echo -e "\r[*] successfully completed command: "
 	>&2 echo -e "\e[31m"$CMD"\e[39m" 
 	fi
fi
cd $BASEDIR
rm -rf split
# Seems like some sequences don't produce .ps files based on 
# line count differences between file_list.txt and fasta header count
# Does not affect clustering, but s few sequences may be missed
# SGE bug or RNAfold error? Folds w/ ViennaRNA-2.2.10

#################################################################
#	                  Convert .ps to .pp 
#################################################################
if [[ ! -d temp ]]; then mkdir temp ; fi
if [[ -e file_list.txt ]]; then rm file_list.txt ; fi
ls ps/ > file_list.txt && >&2 echo -e "[*] saving list of files"
>&2  echo -ne "[*] searching for existing pairwise probability files..."
if [[ -d pp ]]; then 
  if [[ $( ls pp/ | wc -l ) -eq $( ls ps | wc -l ) ]]; then 
    >&2 echo -e "[ \e[92mDONE\e[39m ]"  
    >&2 echo -ne "... detected existing pairwise probability files "
    >&2 echo -ne "    \e[31mDo you want to regenerate them?\e[39m [y/N] "
    read -n1 rePP
    >&2 echo -e
    if [[ $rePP != 'n' && $rePP != 'N' ]]; then 
      NoPP=1
    fi
  fi
else
  >&2 echo -e "[ \e[92mDONE\e[39m ]"  
  mkdir pp
fi
if [[ -z $NoPP ]]; then 
  >&2 echo -e "[*] converting RNA base-pairing probability matrices (patience)..."
  CMD="qsub -V -terse -sync y -cwd -P ${GRID_ACCOUNT} -N makePP 
 -o $BASEDIR/logs/makePP."$(date +"%d%m%y-%H%M")".out 
 -e $BASEDIR/logs/makePP."$(date +"%d%m%y-%H%M")".err 
 -pe smp 1 
 -t 1:${MAXGRIDCPU}
 -b y -S /bin/bash
 ../makePP.sge"
  counter="-1"
  $CMD | while read line ; do 
  	counter=$(( $counter + 1 ))
  	progress=$(( 100 * $counter / $MAXGRIDCPU ))
  	echo -ne "\r... "$progress"% complete" 
  done
fi
## Catch exit status ? 
>&2 echo -e "\r[*] successfully completed command: "
>&2 echo -e "\e[31m"$CMD"\e[39m" 


#################################################################
#	              Generate pairwise comparison list
#################################################################
>&2 echo -en "Generating list of pairwise comparions (patience)...." 
cd ${BASEDIR}/pp
ls *pp | tee ../file_list.txt |\
  awk '
    OFS="\t" { arr[NR]=$0; ++numFiles } 
    END { for ( i=1; i <= numFiles; i++ ) 
      { for ( j=i; j <= numFiles ; j++) 
        print i,j,arr[i],arr[j]
      }
    }
  ' > $BASEDIR/pairwise_list.txt 
>&2 echo -e "... [ \x1B[92mDONE\x1B[39m ]"
cd $BASEDIR


#################################################################
#  SPLIT PAIRWISE LIST INTO $MAXCPU FILES FOR GRID COMPUTING
#################################################################
NUMPW=$(wc -l pairwise_list.txt | awk '{print $1}')
LINESperFILE=$(( ${NUMPW} / ${MAXGRIDCPU} ))
(( LINESperFILE+=1 ))
cd temp
# Clean up previous runs, if present
if [[ ! -d pairwise_groups ]]; then 
  mkdir pairwise_groups  
else 
  rm -rf pairwise_groups  && mkdir pairwise_groups
fi
cd pairwise_groups 
>&2 echo -n "Splitting pairwise comparison list into "
>&2 echo -n $(( 1 + ${NUMPW}/$LINESperFILE ))" files of "
>&2 echo -n $LINESperFILE" lines" 
split -l $LINESperFILE $BASEDIR/pairwise_list.txt 
>&2 echo -e "... [ \x1B[92mDONE\x1B[39m ]"
cd $BASEDIR


################################################################
#                      RUN DOTALIGNER
################################################################
if [[ ! -d hpc ]]; then 
  mkdir hpc 
else 
  rm -rf hpc && mkdir hpc
fi
>&2 echo -e "[*] performing all vs all pairwise comparisons "      
CMD="qsub -V -sync y -cwd -P ${GRID_ACCOUNT} -N DotAlnR -terse
 -o $BASEDIR/logs/launcher.$(date +"%d%m%y-%H%M").out 
 -e $BASEDIR/logs/launcher.$(date +"%d%m%y-%H%M").err 
 -pe smp 1 "
# Dotaligner doesnt need much memory
CMD=${CMD}" -l mem_requested=1256M,h_vmem=1512M
 -t 1:$((1+${NUMPW}/$LINESperFILE)) 
 -b y -S /bin/bash 
 ../worker.sge ${1}"
$CMD | while read line ; do 
  	counter=$(( $counter + 1 ))
  	progress=$(( 100 * $counter / $MAXGRIDCPU ))
  	echo -ne "\r... "$progress"% complete" 
  done
>&2 echo -e "\r[*] completed command: "
>&2 echo -e "\e[31m"$CMD"\e[39m" 



################################################################
#                 FILL IN DISTANCE MATRIX
################################################################
cd $BASEDIR
>&2 echo -n "[*] collating scores..."
cd hpc 
cat *score > ../scores.txt && >&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"

>&2 echo -n "[*] concatenating alignment output..."
cat *out | gzip > ../alignments.out.gz 
>&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"

>&2 echo -n "[*] cleaning up intermediary files..."
cd .. 
rm -rf hpc && >&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"

>&2 echo -n "[*] normalising scores and converting to distance matrix..."
MAXMIN=$(  cat scores.txt | awk 'BEGIN{min = 2 ; max = -1} 
  { if ( $5 > max) {max = $5} else if ( $5 < min ) min = $5 } 
  END {print max" "min}'  )
MAX=${MAXMIN% *}
MIN=${MAXMIN#* }
awk -v min=$MIN -v max=$MAX 'OFS="\t"{ if ( $5 == "-0")  
  print $1,$2,$3,$4,"0" ; 
  else print $1,$2,$3,$4,($5-min)/(max-min) }' scores.txt  \
    > scores_normalized.tsv
awk 'OFS="\t"{ print $1,$2,$3,$4,1-$5 }'  scores_normalized.tsv \
    > dist.tsv
>&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"


################################################################
#                    CLUSTER  ANALYSIS 
################################################################

>&2 echo - "[*] preparing cluster analysis..."
cat > clustering.R << EOF
#!/usr/bin/Rscript
loadPackages <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , quiet = TRUE,
                        repos='http://cran.us.r-project.org' )
      #  Load package after installing
       suppressMessages( suppressWarnings(
                          require( i , character.only = TRUE )))
    }
  }
}
cat("[R] loading required packages...", stderr())
suppressMessages( suppressWarnings(
    loadPackages( c("gplots","dbscan","RcolorBrewer","sparcl") ) ) ) 
cat("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

cat("[R] importing scores...", stderr())
# import dissimilarity scores
d <- read.table( "dist.tsv", header=F )
colnames( d ) <- c( 'x', 'y','x_name','y_name','score' )
cat("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

cat("[R] generating distance matrix...", stderr())
# convert to dissimilarity matrix
md <- matrix( nrow = max( d$y ), ncol=max( d$y ) )
md[ cbind( d$x , d$y ) ] <- as.numeric( d$score )
md[ cbind( d$y , d$x ) ] <- as.numeric( d$score )

# get unique names 
# d$Xname is of format "pp/PROTEIN_CELL_COORDINATES.pp"
# these next 5 lines pull out PROTEIN, which also 
# includes structured RNA controls and decoys
ids <- read.table("ids.txt", header=F)
#ids <- sapply(strsplit(as.character(ids),"/"), "[", 2)
#ids <- sapply(strsplit(as.character(ids),"\\."), "[", 1)

# convert to dist
rownames(md) <- ids$V1
colnames(md) <- ids$V1
D <- as.dist( md )
cat("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

# # hierarchical clustering
# cat("[R] performing hierarchichal clustering ...", stderr())
# hcD <- hclust(D)
# pdf("hclust_dendrogram.pdf")
# 	plot( hcD , cex=0.3, lwd=0.3)
# dev.off()
# cat("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

# cat("[R] extracting hclust clusters ...", stderr())
# clus <- cutree(hcD, h=0.55)
# for  (x in 1:max(clus)) {  
#   c <- labels(which(clus==x)) ; 
#   f <- length( c );
#   u <- length( unique( c ) );  
#   if ( f > 4 && u < 4 &&  u/f <= 0.4 ) {
#       print( x ) ; 
#       print ( c ) ; 
#   }
# }
# cat("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

cat("[R] performing density-based clustering ...", stderr())
#kNNdistplot(D,k=5)
Oc <- opticsXi( optics(D, eps=1, minPts=5, search="dist") , xi = 0.006, minimum=T)

printClust <- function( O ) {
  for (cl in 1:NumClust) {
    print(paste("Cluster ======================== > ",cl)) ; 
    print(labels(D)[ O$cluster == cl ]) ;
  }
}
# printClust(Oc)
# plot(Oc)
# pdf("reachability_plot_optics_xi.pdf")
#   plot(Oc, lwd=0.5)
# dev.off()


# # get a list of optics clusters
# NumClust <- max(Oc$clusters)
# c <- vector("list", NumClust+1 )
# for (cl in 0:NumClust) {
# 	c[[cl+1]] <- labels(D)[ Oc$cluster == cl ]
# }

#===============================================
# 	   		CLUSTER EXTRACTION
# Based on eCLIP-ECS data and merged RFAM benchmark parameters
dir.create( file.path( getwd(), "clusters" ), showWarnings = FALSE)
setwd("clusters")
NumClust <- max(Oc$cluster)
clusters <- 0 #make this 1 if optics(... , minimum=F )
for (cl in 1:NumClust) {
	l <- length(  Oc[Oc$cluster == cl] )
	# extract non-null clusters
	if ( l > 0  ) {
		sequences <- labels(D)[ Oc$cluster == cl ] 
		v <- sapply( strsplit( as.character( sequences ), "_"), "[", 1 ) 
		filename <- paste( "cluster_",cl,".txt", sep="" )
		write.table( sequences, file=filename , col.names=F, row.names=F, quote=F )
		clusters <- clusters + 1
	}
}
cat("Found ",clusters," clusters\n")

cat("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

# #===============================================
# # 	 Plot the reachability w/cluster numbers
# plot(Oc) 
# y <- 0.52 #
# sorted <- sort(unique(Oc$cluster))
# for (x in sorted)  {
# 	if ( x != 0) {
# 		if (  !is.na( as.integer(t[2])) && as.integer( t[2] ) == best ) 
# 			n <- paste( names(t[1]), "-", names(t[2]), sep="" )
# 		else
# 			n <- names(t[1]) 
# 		start <- Oc$clusters_xi$start[  Oc$clusters_xi$cluster_id == x ]	
# 		top <- Oc$reachdist[ start +2 ]
# 		text( start ,  y,  paste(x,n) , col=x+1 , cex=0.9, pos=4, offset=0, srt=25)
# 		segments( start +2 , max( 0.5, top +0.02 ) , 
# 				  start +2, y, 
# 				  col=x+1, lwd=1.5 )
# 		if ( y < 0.61 ) {
# 			y <- y + 0.033
# 		}
# 		else {
# 			y <- 0.52
# 		}
# 		v <- sapply( 
# 				strsplit(
# 				 	as.character( 
# 				 		labels(D)[ Oc$cluster == x ] ), "_"), "[", 1 ) 
# 		t <- sort( table( v ), decreasing=T )
# 		best <- as.integer( t[1] )
# 		# Ignore cluster if it has >1 major representative
# 	}
# }
# #===============================================
# # 	clustering accuracy benchmarking
# mp <- c(3,4,5,6) 
# #Avoid more than 5, unless bypassing the "Minimum" option
# xi <- c(0.005,0.006,0.007,0.008,0.009)
# out <- list
# cat("minPts","Xi","clusters","TP","TN","FP","FN","SENS","SPEC","Minimum","\n",sep="\t", file="cluster_bench.tsv")
# for (MP in mp ) {
# 	for ( Xi in xi ) {
# 		for ( m in c(TRUE,FALSE)) {
# 			#Oc <- opticsXi( optics(D, eps=1, minPts=MP, search="dist") , xi = Xi , minimum=T)
# 			Oc <- opticsXi( optics(D, eps=1, minPts=MP, search="dist") , xi = Xi , minimum = m )
# 			# Get counts of each sequence class
# 			idT <- table( sapply( strsplit( labels( D ), "_"), "[", 1 ) )
# 			C <- list
# 			TP <- list
# 			FP <- list
# 			NumClust <- max(Oc$cluster)
# 			i <- 1 
# 			SENS <- 0
# 			SPEC <- 0 
# 			best <- 0
# 			clusters <- 0
# 			v <- NA
# 			if ( m == TRUE )
# 				first <- 0
# 			else
# 				first <- 1
# 			for (cl in first:NumClust) {
# 				l <- length(  Oc[Oc$cluster == cl] )
# 				# extract non-null clusters
# 				if ( l > 0  ) {
# 					clusters <- clusters + 1
# 					#cat("======= Cluster ",i, "has ",l, "elements =======\n")
# 					C[[ i ]] <- labels(D)[ Oc$cluster == cl ]
# 					# get a vector of names
# 					v <- sapply( strsplit( 
# 							as.character( C[[ i ]] ), "_"), "[", 1 ) 
# 					t <- sort( table( v ), decreasing=T )
# 					best <- as.integer( t[1] )
# 					# Ignore cluster if it has >1 major representative
# 					if ( ( !is.na( best ) && is.na( as.integer(t[2]))) || as.integer( t[2] ) < best ) {
# 						#Get TP and FP
# 						cID <- names( t[ 1 ] )
# 						if ( cID == "shuffled" ) {
# 							FN <- length(v)-best
# 							TN <- best
# 							#cat( cID,TN,FN,"\n",sep="\t")
# 						}
# 						else {
# 							TP[[ i ]] <- best
# 							FP[[ i ]] <- length(v)-best
# 							#cat( cID,best,length(v)-best,"\n",sep="\t")
# 						}
# 						i <- i + 1
# 					}
# 				}
# 			}
# 			SENS = sum( TP, na.rm = T ) / ( sum( TP, na.rm = T ) +FN )
# 			SPEC = TN/(sum( FP, na.rm = T  ) + TN )
# 			cat(MP,Xi,clusters,sum(TP, na.rm = T ),TN,sum(FP,na.rm = T ),FN,SENS,SPEC,m,"\n",sep="\t", file="cluster_bench.tsv", append=T)
# 			#out[[ i ]] <- c(MP, Xi, clusters, SENS, SPEC, m)
# 		}
# 	}
# }
# calculate main cluster component (TP)
# calculate cluster components not in cluster (FN)
# calculate others in cluster (FP)
# calculate others not in cluster (TN)

# gg <- read.table("cluster_bench.tsv",header=T)
# # Youden's J statistic can be calculated in bash with somethig like: 
# # tail -n +2 cluster_bench.tsv | awk 'OFS="\t" {print $1,$2,$10,$3,$8,$9,(($4/($4+$7))+($5/($5+$6))-1)}' | sort -k 7rn
# p <- ggplot(gg) + 
# 	geom_point( data=subset( gg, Minimum == TRUE), aes( 1-SPEC, SENS, color=factor(Xi)
# 		, shape=factor(minPts)), size=2, stroke=1 ) + 
# 	geom_point(  data=subset( gg, Minimum==FALSE), aes( 1-SPEC, SENS, color=factor(Xi), fill=factor(Xi)
# 		, shape=factor(minPts)), size=2, stroke=1, show.legend=F ) + 
# 	scale_y_continuous(lim=c(0,1),breaks=seq(0,1,0.1)) + 
# 	scale_x_continuous(lim=c(0,1),breaks=seq(0,1,0.1)) +
# 	scale_shape_manual(values=c(21,22,23,24,25)) +
# 	scale_color_manual( values=palette(rainbow(6)), name="Steepness\nthreshold Xi") + 
# 	scale_fill_manual(values=palette(rainbow(6))) +
# 	labs( title="Optics clustering accuracy",
# 	x = "False positive rate",
# 	y = "Sensitivty",
# 	color = "Steepness\nthreshold Xi", 
# 	shape = "Minimum\npoints") +
# 	geom_abline(intercept=0, slope=1, linetype=3) 

# l <- ggplot(gg, aes(factor(Xi), clusters, group=interaction(minPts,Minimum), shape=factor(Minimum), color=factor(minPts))) + 
# 	geom_line(aes(linetype=factor(Minimum))) + 
# 	geom_point(size=2, aes(stroke=1)) +
# 	geom_abline(intercept=length(idT), slope=0, linetype=2) + 
# 	guides(shape=F) + 
# 	labs( title="Optics clustering fragmentation",
# 		x = "Steepness threshold Xi",
# 		y = "# of clusters",
# 		color = "Minimum\npoints",
# 		linetype = "Minimum\nCluster\nExtraction") +
# 	theme_bw()


#plot(cut(hcd, h=75)$lower[[2]], 
#     main="Second branch of lower tree with cut at h=75")
# draw the matrix 
#dir.create( "./figs" , showWarnings = FALSE )
# heatmap.2(md, Colv=NA, Rowv=NA, main=hCt, 
#   labCol=D, labRow=D, 
#   symm=T,
#   RowSideColors=as.character(famColours$V2), 
#   ColSideColors=as.character(famColours$V2), 
#   trace='none')
# pdf("./figs/dist_matrix.pdf")
#  heatmap.2( md , 
#       Rowv = as.dendrogram(hcD),
#       Colv="NA", 
#       dendrogram="row", 
#       trace="none", 
#       density="none")
# dev.off()
EOF
>&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"
>&2 echo -ne "[*] launching cluster analysis..."
Rscript --vanilla clustering.R
>&2 echo -e "[*] exiting :) " && exit 0 


# extract fasta records from clusters
for file in cluster_*.txt ; do 
 mkdir ${file%*.txt}
 #change this to mv ? 
 cp $file ${file%*.txt}/
 cd ${file%*.txt}
 #TO DO filter out clusters with identical peaks 
 sed -e 's/_dp\.pp//g' -e 's/_/\.*/g' $file | while read line ; do echo '>'$line ; done > temp 
 grep -A 1 -f temp ../../../../filtered_peaks_RNAbound_ge75_le200_controls.fasta  | grep -v -e '--' > ${file%*.txt}.fasta  
rm temp
~/apps/locarna-1.8.11/bin/mlocarna --probabilistic --iterations=10 \
  --consistency-transformation --threads=6 --noLP ${file%*.txt}.fasta  
 cd ${file%*.txt}.out 
  ~/apps/ViennaRNA-2.2.10/bin/RNAalifold --color --aln -r results/result_prog.aln
  mv alirna.ps ${file%*.txt}_rna.ps
  mv aln.ps ${file%*.txt}_aln.ps
 cd ../..
 rm $file 
done




