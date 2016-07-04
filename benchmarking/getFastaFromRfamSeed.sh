#!/bin/bash
############################################################
# (c) Martin A. Smith 2016
# 	Downloads, splits, and converts RFAM seed alignments to 
# 	fasta format, then stochatically samples sequences 
#	within a certian pairwise identity range ( < and > 55%PID ), 
#	Then generates RNA base pairing dotplot from the sequences.
#
#	Requires ViennaRNA package and GNU coreutils > 2.5
#
############################################################

# load modules (if using Environment Modules for OpenMPI - Rocks)
module load gi/ViennaRNA/2.1.3

########################################
#     	Download RFAM seed 
########################################
if [[ ! -e Rfam.seed ]]; then 
	wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.seed.gz
	gunzip Rfam.seed.gz
fi

########################################
# 	Split alignments into families
########################################
if [[ ! -d ./split ]]; then 
	mkdir split 
fi 
cd split
>&2 echo -n "==> Splitting seed file..."
csplit -n 4 ../Rfam.seed '/^\/\/$/' {*} > /dev/null
# for macs using older BSD. Install homebrew then "brew install coreutils"
# then uncomment this line (commenting out the previous) 
# gcsplit -n 4 ../Rfam.seed '/^\/\/$/' {*} > /dev/null
ls -t | head -n 1 | xargs rm 
>&2 echo "DONE"
cd ..


########################################
#	Convert Stockholm to Fasta
########################################
if [[ ! -d ./fastas ]]; then 
	mkdir fastas
	if [[ -e SelectedRfamIDs.txt ]]; then 
		rm SelectedRfamIDs.txt 
	fi
	>&2 echo "==> Processing Stockholm into fastas"
	fileNum=0
	numFiles=$( ls ./split/x* | wc -l | awk '{print $1}')
	for file in ./split/x* ; do 
		ID=$( grep -m 1 " AC " $file | awk '{print $3}' ) 
		PUB=$( grep -m 1 "Published" $file | wc -l) 
		SS=1
		PSDK=$( grep -m 1 "seudoknot" $file | wc -l )
		NUMS=$( grep -m 1 "GF SQ " $file | \
			awk '{ print $3}' )
		DESC=$( grep -m 1 "GF DE " $file | \
			awk '{ for (i=3;i<=NF;i++) { if (i == NF ) printf $i"\n"; else printf $i" "}}' )
		if [[ $PUB -gt 0 && $SS -gt 0 && $PSDK -eq 0 && $NUMS -gt 10 ]]; then 
			egrep -v -e STOCK -e "//" -e "#=GF" -e RF $file |\
				sed -e 's/#=GC //' -e '/^$/d' |\
				awk '{ printf ">"$1"\n"$2"\n"}'> ./fastas/$ID.fasta
			(echo -ne $ID"\t" ; echo $DESC )>> SelectedRfamIDs.txt
		fi 
		rm $file 
		fileNum=$(( fileNum + 1 ))
		if [[ $( expr $fileNum % 50 ) -eq 0 ]]; then 
			>&2 printf "\r$fileNum of ${numFiles} Files processed" 
		fi
	done
	echo "...DONE"
fi

########################################
>&2 echo -n "==> Extracting less similar [10-55%] sequences for benchmarking..."
#######################################
#### Make a .jar file instead?

java GenerateRFAMsubsets -min_win 70 -max_win 170 -min_pi 10 -max_pi 55 -o ./low_pi_test -i fastas -strip 2> ./low_pi_test/sampling.log
>&2 echo "DONE"
>&2 echo -n "==> Extracting more similar [56-95%]sequences for benchmarking..."
java GenerateRFAMsubsets -min_win 70 -max_win 170 -min_pi 56 -max_pi 95 -o ./high_pi_test -i fastas -strip 2> ./high_pi_test/sampling.log
>&2 echo "DONE"


########################################
#	get one sequence per file. Some tools need this (foldalign)
########################################
>&2 echo  "==> Generating RNA structure matrices for low pairwise identity samples"
cd low_pi_test
if [[ ! -d ./fasta ]]; then mkdir fasta ; else find fasta/* -exec rm {} \; ; fi 
i=0; while read line ; do 
	i=$(( $i + 1 )) 
	case $(( $i % 3 )) in 
		1) file=$( echo $line | cut -c 2-  ) ; echo $line > fasta/${file}.fasta ;; 
		2) echo $line >> fasta/${file}.fasta ;; 
	esac 
done  < output.fasta


########################################
#	Fold to get dotplots
########################################

# use --noLP option? 
RNAfold -p --noPS < output.fasta > RNAfold.log


########################################
#	Convert .ps to .pp 
########################################

if [[ ! -d ./pp ]]; then mkdir pp; fi
if [[ ! -d ./ps ]]; then mkdir ps; fi
if [[ -e ./file_list.txt ]]; then rm file_list.txt ; fi
fileNum=0
for file in *dp.ps ; do 
	java -classpath ../ ps2pp $file > ./pp/${file%*.ps}.pp 
	(( fileNum += 1 ))
	mv ${file} ./ps
	echo `pwd`/pp/${file%*.ps}.pp >> file_list.txt
	>&2 printf "\r$fileNum structures processed into base pairing matrices" 
done
echo "...DONE"


########################################
#	Rinse and Repeat
########################################

>&2 echo "==> Generating RNA structure matrices for high pairwise identity samples"
cd ../high_pi_test 
if [[ ! -d ./fasta ]]; then mkdir fasta ; else find fasta/* -exec rm {} \; ; fi 
# get one sequence per entry
### awk would speed this up
i=0; while read line ; do 
	i=$(( $i + 1 )) 
	case $(( $i % 3 )) in 
		1) file=$( echo $line | cut -c 2-  ) ; echo $line > fasta/${file}.fasta ;; 
		2) echo $line >> fasta/${file}.fasta ;; 
	esac 
done  < output.fasta

########################################

RNAfold -p --noPS < output.fasta > RNAfold.log

########################################

if [[ ! -d ./pp ]]; then mkdir pp; fi
if [[ ! -d ./ps ]]; then mkdir ps; fi
if [[ -e ./file_list.txt ]]; then rm file_list.txt ; fi
fileNum=0
for file in *dp.ps ; do 
	### Modify the source code to run on a folder
	java -classpath ../ ps2pp $file > ./pp/${file%*.ps}.pp 
	(( fileNum += 1 ))
	mv ${file} ./ps
	echo `pwd`/pp/${file%*.ps}.pp >> file_list.txt
	>&2 printf "\r$fileNum structures processed into base pairing matrices" 
done
echo "...DONE"
