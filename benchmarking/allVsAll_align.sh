#!/bin/bash
MAXGRIDCPU=192
TIMEFORMAT="%R"
GRID_ACCOUNT="RNABiologyandPlasticity"

# load modules (for or local server)
module load marsmi/locarna/1.7.13 gi/ViennaRNA/2.1.3 marsmi/dotaligner/0.3

if [[ -e $1/DA.log ]]; then rm $1/DA.log ; fi 
if [[ -e $1/DA_time.log ]]; then rm $1/DA_time.log ; fi 
DABIN=$( which DotAligner ) && >&2 echo "Running: $DABIN" 
k=0
Files=$( ls ${1}/pp/*.pp | tee ${1}/file_list.txt | wc -l )
for (( i=1 ; i < ${Files} ; i++ )); do 
  for (( j = $i+1 ; j <= ${Files}; j++ )); do 
      F1=$( sed $i'q;d' ${1}/file_list.txt )
      F2=$( sed $j'q;d' ${1}/file_list.txt )
      echo -e $i"\t"$j"\t"`pwd`/${F1//pp/ps}"\t"`pwd`/${F2//pp/ps} 
      echo -ne $i"\t"$j"\t" >> $1/DA.log 
      ( time $DABIN -d $F1 -d $F2 -o 0.4 -e 0.4 -t 0.2 >> $1/DA.log ) 2>> $1/DA_time.log
      ((k+=1))   
      if [[ $(( $k % 50 )) -eq 0 ]]; then 
        >&2 printf "\r%.1f%s" $( echo $k $Files | \
          awk '{printf 100*$1/($2*($2-1)/2)}' ) "% of alignments completed"
      fi
  done
done > $1/pairwise_list.txt
echo

#############################################
#  MAKE BINARY MATRIX FOR CORRELATION
#############################################
>&2  echo -n "[NOTE] generating matrix..."
LINES=$( wc -l $1/DA.log )
k=1 
###   Print out column names 
(sed 's/^.*\///g' low_pi/file_list.txt | cut -d "_" -f 1,2 | while read line ; 
  do  echo -n $line" " 
  done ; echo ) | sed 's/.$//' > $1/binary.matrix

for (( i=1 ; i <= $Files ; i++ )); do 
  S1=$(sed $i'q;d' $1/file_list.txt | sed 's/^.*\///' | cut -d "_" -f 1 )     
  for (( j=1 ; j <= $Files; j++ )); do 
    S2=$(sed $j'q;d' $1/file_list.txt | sed 's/^.*\///' | cut -d "_" -f 1 ) 
    ###   Print out row names 
    #if [[ $j -eq 1 ]]; then 
    #    echo -n $S1" " >> $1/binary.matrix
    if [ $S1 == $S2 ]; then 
      echo -n "1 " 
    else
      echo -n "0 " 
    fi
  done
  echo 
done | sed 's/.$//' >>  $1/binary.matrix
>&2 echo "DONE"

#################################################################
# SPLIT PAIRWISE LIST INTO $MAXCPU FILES FOR GRID COMPUTING
#################################################################
NUMPW=$(wc -l $1/pairwise_list.txt | awk '{print $1}')
MAXGRIDCPU2=$(( $MAXGRIDCPU +1 )) 
NUMFILES=$(( ${NUMPW} / ${MAXGRIDCPU2} + 1 ))
cd $1
if [[ ! -d pairwise_groups ]]; then mkdir pairwise_groups ; fi
 cd pairwise_groups 
 split -l $NUMFILES ../pairwise_list.txt 
 cd ..
if [[ ! -d carna ]]; then mkdir carna ; else rm -f carna/* ; fi
if [[ ! -d locarna ]]; then mkdir locarna ; else rm -f locarna/* ; fi
cd ..
#################################################################
#                         RUN CARNA
#################################################################
CARNA=$( qsub -V -o $1/carna/qsub.log -cwd -P $GRID_ACCOUNT -N CARNA -j y -pe smp 1 -l mem_requested=2G,h_vmem=2G -t 1:${MAXGRIDCPU} -b y ./carna.sge $1 )
echo $CARNA
#Example output: 
#Your job-array 4425908.1-192:1 ("CARNA") has been submitted
# N.B. qsub option "-terse" only outputs the job number
JID=$( echo $CARNA | cut -d " " -f 3 | cut -d "." -f 1 )
qsub -hold_jid $JID -V -cwd -N clean_carna -j y -pe smp 1 -l mem_requested=1G,h_vmem=1G -b y ./cleanup.sge $1 carna

#################################################################
#                         RUN LOCARNA
#################################################################
LOCARNA=$( qsub -V -o $1/locarna/qsub.log -cwd -P $GRID_ACCOUNT -N LOCARNA -j y -pe smp 1 -l mem_requested=2G,h_vmem=2G -t 1:${MAXGRIDCPU} -b y ./locarna.sge $1 ) 
echo $LOCARNA
JID=$( echo $LOCARNA | cut -d " " -f 3 | cut -d "." -f 1 )
qsub -hold_jid $JID -V -cwd -N clean_locarna -j y -pe smp 1 -l mem_requested=1G,h_vmem=1G -b y ./cleanup.sge $1 locarna

#################################################################
#                         CLEAN UP 
#################################################################





