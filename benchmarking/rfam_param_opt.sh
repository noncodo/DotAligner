#!/bin/bash
MAXGRIDCPU=179
TIMEFORMAT="%R"
GRID_ACCOUNT="RNABiologyandPlasticity"

# load modules (for or local server)
module load gi/ViennaRNA/2.1.3 marsmi/dotaligner/0.3

# >&2 echo "Generating list of pairwise comparions...." 
# k=0
# Files=$( ls ${1}/pp/*.pp | tee ${1}/file_list.txt | wc -l )
# for (( i=1 ; i <= ${Files} ; i++ )); do 
#   for (( j = $i ; j <= ${Files}; j++ )); do 
#       F1=$( sed $i'q;d' ${1}/file_list.txt )
#       F2=$( sed $j'q;d' ${1}/file_list.txt )
#       echo -e $i"\t"$j"\t"`pwd`/${F1//pp/ps}"\t"`pwd`/${F2//pp/ps} 
#       ((k+=1))   
#       if [[ $(( $k % 50 )) -eq 0 ]]; then 
#         >&2 printf "\r%.1f%s" $( echo $k $Files | \
#           awk '{printf 100*$1/($2*($2+1)/2)}' ) "% of alignments completed"
#       fi
#   done
# done > $1/pairwise_list.txt
# >&2 echo  
# cut -d "/" -f 3 $1/file_list.txt | cut -d "_" -f 1-2 > $1/structure_ids.txt
if [[ ! -d $1/logs ]]; then mkdir $1/logs ; fi
if [[ ! -d $1/bench ]]; then mkdir $1/bench ; else rm $1/bench/* ; fi

for k in `seq 0.25 0.25 0.75`; do  
 for t in `seq 0.25 0.25 0.75`; do 
  for s in 5000 1000 200; do  
   for T in 10 5 1; do
    for o in 0.75 1; do  
     for e in 0.05 0.2; do
      CMD="DotAligner -T ${T} -s ${s} -k ${k} -t ${t} -o ${o} -e ${e} "
      PARAMS="T-${T}_s-${s}_k-${k}_t-${t}_o-${o}_e-${e}"
      echo $CMD >> $1/bench/cmds.list
	  echo $PARAMS >> $1/bench/out.list
	 done
	done
   done
  done
 done
done 
echo "Generated "`wc -l $1/bench/cmds.list`" commands"

timestamp=$( date +"%d%m%y-%H%M" )
#CMD=$( echo qsub -cwd -V -S /bin/bash  -o $1/bench_qsub.out.${timestamp} -e $1/bench_qsub.err.${timestamp} -P RNABiologyandPlasticity -t 1:`wc -l ${1}/bench/cmds.list` -V -cwd -N bench_low -l mem_requested=1G,h_vmem=2G -j y -b y ./bench_it_high.sge  )
cmds=$(wc -l ${1}/bench/cmds.list | awk '{print $1}' )
CMD="qsub -cwd -V -S /bin/bash \
	-o $1/logs/qsub_bench.out.${timestamp} \
	-e $1/logs/qsub_bench.err.${timestamp} \
	-P RNABiologyandPlasticity \
	-t 1:${cmds} \
	-V -cwd \
	-N bench_${1} \
	-l mem_requested=1.5G,h_vmem=2G \
	-j y -b y \
	./da_bench_${1}.sge" 
#	-t 1:`wc -l ${1}/bench/cmds.list`
(echo "[COMMAND]"; echo $CMD) && $CMD 
#echo $CMD 

#qsub  -t 1:`wc -l ${1}/bench/cmds.list` -V -cwd -N bench_high -l mem_requested=1G,h_vmem=2G -j y -b y ./bench_it_high.sge

# 
# if [[ -e ${1}/qsub.dotalnr.out ]] ; then rm ${1}/qsub.dotalnr.out ; fi 
# if [[ -e ${1}/qsub.dotalnr.err ]] ; then rm ${1}/qsub.dotalnr.err ; fi 
# CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N DotAlnR -o ${1}/logs/dotaligner.$(date +"%d%m%y-%H%M").out -e ${1}/logs/dotaligner.$(date +"%d%m%y-%H%M").err -pe smp 1"
# CMD=${CMD}" -l mem_requested=1256M,h_vmem=1512M,h=!(epsilon-0-24.local|epsilon-0-25.local) -t 1:$((1+${NUMPW}/$LINESperFILE)) -b y ./dotaligner.sge ${1}"
# DOTALNR=$( $CMD ) 
# echo "launching: "$CMD 
# echo $DOTALNR
# JID=$( echo ${DOTALNR} | cut -d " " -f 3 | cut -d "." -f 1 )
# qsub -hold_jid $JID -V -o ./$1/logs/dotaligner.$(date +"%d%m%y-%H%M").clean -cwd -N clean_dotalnr -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 dotaligner

