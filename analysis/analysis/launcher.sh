#!/bin/bash

#load binaries 
# TO DO: Make more user friendy (i.e. check for presence)
module load gi/ViennaRNA/2.1.3 marsmi/dotaligner/0.3 
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
#	Fold RNA to get dotplots
#################################################################
foldit () {
 echo -e "[DEBUG] Funciton WD = "$( pwd ) 
 >&2 echo -ne "Generating RNA base-pairing probability matrices "
 >&2 echo -ne "[ this takes about 1m per 500 seqs ] ..." 
 ## TO DO  split file into $MAXCPU subfiles to speed this up
 $RNAFOLDBIN -p --noPS --noLP < $1 > ../logs/RNAfold.log
 >&2 echo -e "DONE"  
}
## TO DO : parallelize to make it faster 
echo -e "[DEBUG] BASEDIR= "$BASEDIR
if [[ ! -d $BASEDIR/logs ]]; then mkdir -p $BASEDIR/logs ; fi
cd $BASEDIR	
if [[ ! -d ps ]]; then 
	mkdir ps && cd ps
	echo -e "[DEBUG] making PS folder. Running foldit "${BASEDIR%/*}/$1
	foldit ${BASEDIR%/*}/$1
elif [[ $( ls ps | wc -l ) -eq $( grep ">" ../${1} | wc -l ) ]]; then 
	echo -e "\e[33mDetected existing RNA structures."
  echo -ne "Do you want to refold them? (y/N)\e[0m"
	read -n1 redo 
	echo -e
	if [[ $redo != 'n' && $redo != 'N' ]]; then 
 		cd ps && foldit ${BASEDIR%/*}/$1
 		cd ..
 	fi
fi
cd $BASEDIR


#################################################################
#	Convert .ps to .pp 
#################################################################
if [[ ! -d temp ]]; then mkdir temp ; fi
if [[ -e file_list.txt ]]; then rm file_list.txt ; fi
if [[ -d pp ]]; then 
  if [[ $( ls pp | wc -l ) -eq $( ls ps | wc -l ) ]]; then 
    echo -e "\e[33mDetected existing pairwise probability files."
    echo -ne "Do you want to regenerate them? (y/N)\e[0m "
    read -n1 rePP
    echo -e
    if [[ $rePP != 'n' && $rePP != 'N' ]]; then 
      NoPP=1
    fi
  fi
else
  mkdir pp
fi
if [[ -z $NoPP ]]; then 
  >&2 echo -e "Converting RNA base-pairing probability matrices "
  CMD="qsub -V -sync y -cwd -P ${GRID_ACCOUNT} -N makePP 
 -o $BASEDIR/logs/makePP.$(date +"%d%m%y-%H%M").out 
 -e $BASEDIR/logs/makePP.$(date +"%d%m%y-%H%M").err 
 -pe smp 1 
 -l mem_requested=2G,h_vmem=2G
 -t 1:${MAXGRIDCPU}
 -b y -S /bin/bash
 ../makePP.sge"
  counter=0
  echo "[ CMD ] "$CMD 
  echo "[ be patient ]"
  $CMD | head -n 1 
  for dir in temp/* ; do
	 mv $dir/*pp pp/ ; 
   rm -rf $dir
  done
fi


#################################################################
#	Generate pairwise comparison list
#################################################################
>&2 echo -n "Generating list of pairwise comparions...." 
k=0
cd pp
ls *pp | tee ../file_list.txt |\
  awk '
    OFS="\t" { arr[NR]=$0; ++numFiles } 
    END { for ( i=1; i <= numFiles; i++ ) 
      { for ( j=i; j <= numFiles ; j++) 
        print i,j,arr[i],arr[j]
      }
    }
  ' > $BASEDIR/pairwise_list.txt 
>&2 echo "DONE"
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
else rm pairwise_groups/* 
fi
cd pairwise_groups 
>&2 echo "Splitting pairwise comparison list into \
	"$((1+${NUMPW}/$LINESperFILE))" files of "$LINESperFILE" lines" \
	&& split -l $LINESperFILE $BASEDIR/pairwise_list.txt 
cd $BASEDIR


################################################################
#                      RUN DOTALIGNER
################################################################
if [[ ! -d $BASEDIR/hpc ]]; then mkdir $BASEDIR/hpc ; fi      
CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N DotAlnR -terse
 -o ${BASEDIR}/logs/launcher.$(date +"%d%m%y-%H%M").out 
 -e ${BASEDIR}/logs/launcher.$(date +"%d%m%y-%H%M").err 
 -pe smp 1 "
# Dotaligner doesnt need much memory
CMD=${CMD}" -l mem_requested=1256M,h_vmem=1512M
 -t 1:$((1+${NUMPW}/$LINESperFILE)) 
 -b y -S /bin/bash 
 ${BASEDIT}/worker.sge ${1}"
#JID=$( $CMD ) 
echo "launching: "$CMD 
JID=$( echo ${JID} | cut -d "." -f 1 )
#qsub -hold_jid $JID -V -o ./$1/logs/dotaligner.$(date +"%d%m%y-%H%M").clean -cwd -N clean_dotalnr -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 dotaligner
