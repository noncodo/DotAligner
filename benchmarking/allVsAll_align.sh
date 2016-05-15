#!/bin/bash
MAXGRIDCPU=179
TIMEFORMAT="%R"
GRID_ACCOUNT="RNABiologyandPlasticity"

# load modules (for or local server)
module load marsmi/locarna/1.7.13 gi/ViennaRNA/2.1.3 marsmi/dotaligner/0.3

>&2 echo "Generating list of pairwise comparions...." 
k=0
Files=$( ls ${1}/pp/*.pp | tee ${1}/file_list.txt | wc -l )
for (( i=1 ; i <= ${Files} ; i++ )); do 
  for (( j = $i ; j <= ${Files}; j++ )); do 
      F1=$( sed $i'q;d' ${1}/file_list.txt )
      F2=$( sed $j'q;d' ${1}/file_list.txt )
      echo -e $i"\t"$j"\t"`pwd`/${F1//pp/ps}"\t"`pwd`/${F2//pp/ps} 
      ((k+=1))   
      if [[ $(( $k % 50 )) -eq 0 ]]; then 
        >&2 printf "\r%.1f%s" $( echo $k $Files | \
          awk '{printf 100*$1/($2*($2+1)/2)}' ) "% of alignments completed"
      fi
  done
done > $1/pairwise_list.txt
>&2 echo  
cut -d "/" -f 3 $1/file_list.txt | cut -d "_" -f 1-2 > $1/structure_ids.txt

## # #############################################
## # #  MAKE BINARY MATRIX FOR CORRELATION
## # #############################################
## >&2  echo -n "[NOTE] generating matrix..."
## LINES=$( wc -l $1/DA.log )
## k=1 
## ###   Print out column names 
## (sed 's/^.*\///g' low_pi/file_list.txt | cut -d "_" -f 1,2 | while read line ; 
##   do  echo -n $line" " 
##   done ; echo ) | sed 's/.$//' > $1/binary.matrix
#
## for (( i=1 ; i <= $Files ; i++ )); do 
##   S1=$(sed $i'q;d' $1/file_list.txt | sed 's/^.*\///' | cut -d "_" -f 1 )     
##   for (( j=1 ; j <= $Files; j++ )); do 
##     S2=$(sed $j'q;d' $1/file_list.txt | sed 's/^.*\///' | cut -d "_" -f 1 ) 
##     ###   Print out row names 
##     #if [[ $j -eq 1 ]]; then 
##     #    echo -n $S1" " >> $1/binary.matrix
##     if [ $S1 == $S2 ]; then 
##       echo -n "1 " 
##     else
##       echo -n "0 " 
##     fi
##   done
##   echo 
## done | sed 's/.$//' >>  $1/binary.matrix
## >&2 echo "DONE"

#################################################################
# SPLIT PAIRWISE LIST INTO $MAXCPU FILES FOR GRID COMPUTING
#################################################################
NUMPW=$(wc -l $1/pairwise_list.txt | awk '{print $1}')
LINESperFILE=$(( ${NUMPW} / ${MAXGRIDCPU} ))
(( LINESperFILE+=1 ))
cd $1
if [[ ! -d pairwise_groups ]]; then 
  mkdir pairwise_groups  
else rm pairwise_groups/* 
fi
 cd pairwise_groups 
 >&2 echo "Splitting pairwise comparison list into "$((1+${NUMPW}/$LINESperFILE))" files of "$LINESperFILE" lines" \
  && split -l $LINESperFILE ../pairwise_list.txt 
 cd ..
cd ..

if [[ ! -d $1/logs ]]; then mkdir $1/logs ; fi
#################################################################
#                       RUN DOTALIGNER
#################################################################
if [[ ! -d $1/dotaligner ]]; then mkdir $1/dotaligner ; fi      # pre-existing files get removed in the sub-commands
if [[ -e ${1}/qsub.dotalnr.out ]] ; then rm ${1}/qsub.dotalnr.out ; fi 
if [[ -e ${1}/qsub.dotalnr.err ]] ; then rm ${1}/qsub.dotalnr.err ; fi 
CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N DotAlnR -o ${1}/logs/dotaligner.$(date +"%d%m%y-%H%M").out -e ${1}/logs/dotaligner.$(date +"%d%m%y-%H%M").err -pe smp 1"
CMD=${CMD}" -l mem_requested=1256M,h_vmem=1512M,h=!(epsilon-0-24.local|epsilon-0-25.local) -t 1:$((1+${NUMPW}/$LINESperFILE)) -b y ./dotaligner.sge ${1}"
DOTALNR=$( $CMD ) 
echo "launching: "$CMD 
echo $DOTALNR
JID=$( echo ${DOTALNR} | cut -d " " -f 3 | cut -d "." -f 1 )
qsub -hold_jid $JID -V -o ./$1/logs/dotaligner.$(date +"%d%m%y-%H%M").clean -cwd -N clean_dotalnr -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 dotaligner


#################################################################
#                         RUN LOCARNA
#################################################################
if [[ ! -d $1/locarna ]]; then mkdir $1/locarna ; fi            # pre-existing files get removed in the sub-commands
if [[ -e ${1}/qsub.locarna.out ]] ; then rm ${1}/qsub.locarna.out ; fi 
if [[ -e ${1}/qsub.locarna.err ]] ; then rm ${1}/qsub.locarna.err ; fi 
CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N LOCARNA -o ${1}/logs/locarna.$(date +"%d%m%y-%H%M").out -e ${1}/logs/locarna.$(date +"%d%m%y-%H%M").err -pe smp 1"
CMD=${CMD}" -l mem_requested=1792M,h_vmem=2048M -t 1:$((1+${NUMPW}/$LINESperFILE)) -b y ./locarna.sge ${1}"
LOCARNA=$( $CMD ) 
echo "launching: "$CMD 
echo $LOCARNA
JID=$( echo ${LOCARNA} | cut -d " " -f 3 | cut -d "." -f 1 )
qsub -hold_jid $JID -V -o ./$1/logs/locarna.$(date +"%d%m%y-%H%M").clean -cwd -N clean_locarna -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 locarna


#################################################################
#                         RUN CARNA
#################################################################
 if [[ ! -d $1/carna ]]; then mkdir $1/carna ; fi                 # pre-existing files get removed in the sub-commands
 if [[ -e ${1}/qsub.carna.out ]] ; then rm ${1}/qsub.carna.out ; fi 
 if [[ -e ${1}/qsub.carna.err ]] ; then rm ${1}/qsub.carna.err ; fi 
 CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N CARNA -o ${1}/logs/carna.$(date +"%d%m%y-%H%M").out -e ${1}/logs/carna.$(date +"%d%m%y-%H%M").err -pe smp 1"
 CMD=${CMD}" -l  mem_requested=2792M,h_vmem=3048M,h=!(epsilon-0-24.local|epsilon-0-25.local) -t 1:$((1+${NUMPW}/${LINESperFILE})) -b y ./carna.sge ${1}"
 CARNA=$( $CMD ) 
 echo "launching: "$CMD
 echo $CARNA
 #Example output: 
 #Your job-array 4425908.1-192:1 ("CARNA") has been submitted
 # N.B. qsub option "-terse" only outputs the job number
 JID=$( echo ${CARNA} | cut -d " " -f 3 | cut -d "." -f 1 )
 qsub -hold_jid $JID -V -o ./$1/logs/carna.$(date +"%d%m%y-%H%M").clean -cwd -N clean_carna -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 carna


# #################################################################
# #                             R
# #################################################################

# # R CMD BATCH [options] my_script.R [outfile]
# # Rscript myScript.R 5 100
# ## myScript.R
# ## args <- commandArgs(trailingOnly = TRUE)
# ## rnorm(n=as.numeric(args[1]), mean=as.numeric(args[2]))


