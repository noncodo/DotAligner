#if [[ -e $1/DA.log ]]; then rm $1/DA.log ; fi 
Files=$( ls $1/pp/*.pp | wc -l )
TIMEFORMAT="%R"
# k=0
# for (( i=1 ; i < $Files ; i++ )); do 
#   for (( j = $i+1 ; j <= $Files; j++ )); do 
#       F1=$( ls $1/pp/*.pp | sed $i'q;d' )
#       F2=$( ls $1/pp/*.pp | sed $j'q;d' )
#       echo -ne $i"\t"$j"\t" >> $1/DA.log 
#       ( time /Users/martinsmith/apps/dotaligner_v0_3/src/DotAligner -d $F1 -d $F2 -o 0.4 -e 0.4 -t 0.2 >> $1/DA.log ) 2>> DA_time.log
#       #echo "{ time ${EXE} -d $FAMSEQ1 -d $FAMSEQ2 -k 0.5 -a -0.05 -b -0.05 -r 5 -t 0.2 -l 0.5 -s 0 -m 15; } &>> ${DATA}.${PROG};"
#       # DotAligber output format
#       # seq_score probability norm_seq_score  norm_str_score  combined_score  length  alignment
#       ((k+=1))   
#       if [[ $(( $k % 10 )) -eq 0 ]]; then 
#         >&2 printf "\r%.1f%s" $( echo $k $Files | awk '{printf 100*$1/($2*($2-1)/2)}' ) "% of alignments completed"
#       fi
#   done
# done
# echo

#############################################
#  MAKE BINARY MATRIX FOR CORRELATION
#############################################
>&2  echo -n "[NOTE] generating matrix..."
LINES=$( wc -l $1/DA.log )
k=1 
###   Print out column names 
(sed 's/^.*\///g' low_pi/file_list.txt | cut -d "_" -f 1,2 | while read line ; do echo -n $line" " ; done ; echo ) |\
 sed 's/.$//' > $1/binary.matrix
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