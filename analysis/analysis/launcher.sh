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
logo() {
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
                                    
}

foldit () {
 #>&2 echo -e "[DEBUG] Funciton WD = "$( pwd ) 
 >&2 echo -ne "[*] generating RNA base-pairing probability matrices "
 >&2 echo -ne "(\e[31mbe patient\e[39m)..." 
 ## TO DO  split file into $MAXCPU subfiles to speed this up
 $RNAFOLDBIN -p --noPS --noLP < $1 > ../logs/RNAfold.log
 >&2 echo -e "[ \e[92mDONE\e[39m ]"  
}


#################################################################
#                 Fold RNA to get dotplots
#################################################################
## TO DO : parallelize to make it faster 
<<<<<<< HEAD
#>&2 echo -e "[DEBUG] BASEDIR= "$BASEDIR
=======
>&2 echo -e "[DEBUG] BASEDIR= "$BASEDIR
>>>>>>> 8d5862964db2eeac240e465312a44c46c9723858
if [[ ! -d $BASEDIR/logs ]]; then mkdir -p $BASEDIR/logs ; fi
cd $BASEDIR	&& >&2 echo -e "[*] working in BASEDIR: "$(pwd)
>&2 echo -ne "[*] Evaluating folder contents ..."
if [[ ! -d ps ]]; then 
  	>&2 echo -e "no RNA structures found"
	mkdir ps && cd ps
	foldit ${BASEDIR%/*}/$1
elif [[ $( ls ps | wc -l ) -ne $( grep ">" ../${1} | wc -l ) ]]; then 
	>&2 echo -e "incomplete amount of structures"
	cd ps
	foldit ${BASEDIR%/*}/$1
else
	>&2 echo -e "found existing RNA structures "
  	>&2 echo -ne "    \e[31mDo you want to refold them?\e[39m [y/N] "
	read -n1 redo 
	>&2 echo 
	if [[ $redo != 'n' && $redo != 'N' ]]; then 
 		cd ps && foldit ${BASEDIR%/*}/$1
 		cd ..
 	fi
fi
cd $BASEDIR


#################################################################
#	                  Convert .ps to .pp 
#################################################################
if [[ ! -d temp ]]; then mkdir temp ; fi
if [[ -e file_list.txt ]]; then rm file_list.txt ; fi
<<<<<<< HEAD
ls ps/ > file_list.txt && >&2 echo -e "[*] saving list of files"
>&2  echo -ne "[*] searching for existing pairwise probability files..."
=======
#Avoid 'ls' too many times in HPC
ls ps/ > file_list.txt
>>>>>>> 8d5862964db2eeac240e465312a44c46c9723858
if [[ -d pp ]]; then 
  if [[ $( ls pp/ | wc -l ) -eq $( ls ps | wc -l ) ]]; then 
    >&2 echo -e "[ \e[92mDONE\e[39m ]"  
    >&2 echo -ne "... detected existing pairwise probability files "
    >&2 echo -ne "    \e[31mDo you want to regenerate them?\e[39m [y/N] "
    read -n1 rePP
<<<<<<< HEAD
    >&2 echo -e
=======
    echo -e
>>>>>>> 8d5862964db2eeac240e465312a44c46c9723858
    if [[ $rePP != 'n' && $rePP != 'N' ]]; then 
      NoPP=1
    fi
  fi
else
  >&2 echo -e "[ \e[92mDONE\e[39m ]"  
  mkdir pp
fi
if [[ -z $NoPP ]]; then 
<<<<<<< HEAD
  >&2 echo -e "[*] converting RNA base-pairing probability matrices (patience)..."
  CMD="qsub -V -terse -sync y -cwd -P ${GRID_ACCOUNT} -N makePP 
 -o $BASEDIR/logs/makePP."$(date +"%d%m%y-%H%M")".out 
 -e $BASEDIR/logs/makePP."$(date +"%d%m%y-%H%M")".err 
=======
  >&2 echo -e "Converting RNA base-pairing probability matrices "
  # N.B. Java VM requires close to 7.5GB of RAM 
  CMD="qsub -V -sync y -cwd -P ${GRID_ACCOUNT} -N makePP 
 -o $BASEDIR/logs/makePP.$(date +"%d%m%y-%H%M").out 
 -e $BASEDIR/logs/makePP.$(date +"%d%m%y-%H%M").err 
>>>>>>> 8d5862964db2eeac240e465312a44c46c9723858
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
<<<<<<< HEAD
if [[ ! -d hpc ]]; then 
  mkdir hpc 
else 
  rm -rf hpc && mkdir hpc
fi
>&2 echo -e "[*] performing all vs all pairwise comparisons "      
CMD="qsub -V -sync y -cwd -P ${GRID_ACCOUNT} -N DotAlnR -terse
 -o $BASEDIR/logs/launcher.$(date +"%d%m%y-%H%M").out 
 -e $BASEDIR/logs/launcher.$(date +"%d%m%y-%H%M").err 
=======
if [[ ! -d $BASEDIR/hpc ]]; then mkdir $BASEDIR/hpc ; fi      
CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N DotAlnR -terse
 -o ${BASEDIR}/logs/launcher.$(date +"%d%m%y-%H%M").out 
 -e ${BASEDIR}/logs/launcher.$(date +"%d%m%y-%H%M").err 
>>>>>>> 8d5862964db2eeac240e465312a44c46c9723858
 -pe smp 1 "
# Dotaligner doesnt need much memory
CMD=${CMD}" -l mem_requested=1256M,h_vmem=1512M
 -t 1:$((1+${NUMPW}/$LINESperFILE)) 
 -b y -S /bin/bash 
<<<<<<< HEAD
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
loadPackage <- function(x){
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
    loadPackage( c("gplots","dbscan") ) ) ) 
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
# includes known strucutred RNA controls 
ids <- unique(d$x_name)
ids <- sapply(strsplit(as.character(ids),"/"), "[", 2)
ids <- sapply(strsplit(as.character(ids),"\\."), "[", 1)

#get colors
palette(rainbow(length(unique(ids))))
colourMap <- cbind( unique(ids), palette() )
famColours <- merge(d$x_name, colourMap, by=1)

# convert to dist
rownames(md) <- ids
colnames(md) <- ids
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
#DB <- dbscan(D,eps=0.44,minPts=5)
#O <- optics(D, eps=0.5, minPts=5, search="dist")

printClust <- function( O ) {
  for (cl in 1:NumClust) {
    print(paste("Cluster ======================== > ",cl)) ; 
    print(labels(D)[ O$cluster == cl ]) ;
  }
}

Oc <- opticsXi( optics(D, eps=1, minPts=5, search="dist") , 
		xi = 0.01, minimum=T)
printClust(Oc)
plot(Oc)
pdf("reachability_plot_optics_xi.pdf")
  plot(Oc, lwd=0.5)
dev.off()

cat("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

# # get a list of optics clusters
# NumClust <- max(Oc$clusters)
# c <- vector("list", NumClust+1 )
# for (cl in 0:NumClust) {
# 	c[[cl+1]] <- labels(D)[ Oc$cluster == cl ]
# }

# Get counts of each sequence class
idT <- table( sapply( strsplit( labels( D ), "_"), "[", 1 ) )
c <- list
TP <- list
FP <- list
NumClust <- max(Oc$clusters)
i <- 1 
for (cl in 0:NumClust) {
	l <- length(  Oc[Oc$cluster == cl] )
	# extract non-null clusters
	if ( l > 0  ) {
		cat("Cluster ",i-1, "has ",l, "elements; ")
		c[[ i ]] <- labels(D)[ Oc$cluster == cl ]
		# get a vector of names
		v <- sapply( strsplit( 
				as.character( c[[ i ]] ), "_"), "[", 1 ) 
		t <- sort( table( v ), decreasing=T )
		best <- as.integer( t[1] )
		# 1  RF00002    7
		# 2 shuffled    2
		# Ignore cluster if it has >1 major representative
		if ( as.integer( t[2] ) == best ) {
			cat( "but no major representative [ignoring it]")
			i <- i + 1 
			skip
		}
		else {
			#Get TP and FP
			cID <- names( t[ 1 ] )
			if ( cl == 0 ) {
				FN <- length(v)-best
				TN <- best
			}
			else {
				#cat( cID, " " )
				TP[[ i - 1 ]] <- best
				cat(best,"\n")
				FP[[ i -1 ]] <- length(v)-best
			}
		}
		i <- i + 1
	}
}
SENS = sum(TP)/(sum(TP)+FN)
SPEC = TN/(sum(FP)+TN)
#FN <- list(FN, (as.integer( idT[ names(t[1]) ] ) - best ))
			#TN <- TN[[ i ]] <- length(ids) - TP - FP - FN)
			#cat( "TP ",TP," FP ",FP," TN ",TN," FN ",FN )
			#cat( " SENS ",TP/(TP+FN)," SPEC ",TN/(TN+FP),"\n")

sort( as.data.frame( table(  ) )$Freq , decreasing=T)[1]

# calculate main cluster component (TP)
# calculate cluster components not in cluster (FN)
# calculate others in cluster (FP)
# calculate others not in cluster (TN)





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

=======
 ${BASEDIR}/../worker.sge ${1}"
JID=$( $CMD ) 
echo "launched: "$CMD 
JID=$( echo ${JID} | cut -d "." -f 1 )
#qsub -hold_jid $JID -V -o ./$1/logs/dotaligner.$(date +"%d%m%y-%H%M").clean -cwd -N clean_dotalnr -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 dotaligner
>>>>>>> 8d5862964db2eeac240e465312a44c46c9723858
