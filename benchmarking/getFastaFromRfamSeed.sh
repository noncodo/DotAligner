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

# load modules (for or local server)
module load gi/ViennaRNA/2.1.3

if [[ ! -e Rfam.seed ]]; then 
	wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.seed.gz
	gunzip Rfam.seed.gz
fi
if [[ ! -d ./split ]]; then 
	mkdir split && cd split
	>&2 echo -n "==> Splitting seed file..."
	csplit -n 4 ../Rfam.seed '/^\/\/$/' {*} > /dev/null
	ls -t | head -n 1 | xargs rm 
	>&2 echo "DONE"
	cd ..
fi
if [[ ! -d ./fastas ]]; then 
	mkdir fastas
	if [[ -e SelectedRfamIDs.txt ]]; then 
		rm SelectedRfamIDs.txt 
	fi
	>&2 echo "==> Processing Stockholm into fastas"
	fileNum=0
	for file in ./split/x* ; do 
		ID=$( grep -m 1 " AC " $file | awk '{ print $3}' ) 
		PUB=$( grep -m 1 "Published" $file | wc -l) 
	#	SS=$( grep -m 1 -e "rystal" -e "ertiary" $file | wc -l )
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
			>&2 printf "\r$fileNum Files processed" 
		fi
	done
	echo
fi

## pccasional GenerateRFAMsubsets error on last sequence (rfam seqs = 1) 
>&2 echo  "==> Extracting less similar [10-55%] sequences for benchmarking..."
java GenerateRFAMsubsets -min_win 70 -max_win 170 -min_pi 10 -max_pi 55 -o ./low_pi -i fastas -strip
>&2 echo "DONE"
>&2 echo "==> Extracting more similar [56-95%]sequences for benchmarking..."
java GenerateRFAMsubsets -min_win 70 -max_win 170 -min_pi 56 -max_pi 95 -o ./high_pi -i fastas -strip
>&2 echo "DONE"

>&2 echo  "==> Generating RNA structure matrices for low pairwise identity samples"

cd low_pi && RNAfold -p --noPS < output.fasta > RNAfold.log
if [[ ! -d ./pp ]]; then mkdir pp; fi
if [[ ! -d ./ps ]]; then mkdir ps; fi
if [[ -e ./file_list.txt ]]; then rm file_list.txt ; fi
fileNum=0
for file in *dp.ps ; do 
	java -classpath ../ ps2pp $file > ./pp/${file%*.ps}.pp 
	(( fileNum += 1 ))
	mv ${file} ./ps
	echo `pwd`/pp/${file%*.ps}.pp >> file_list.txt
	>&2 printf "\r$fileNum lower PI structures processed into base pairing matrices" 
done
echo

>&2 echo "==> Generating RNA structure matrices for high pairwise identity samples"
cd ../high_pi && RNAfold -p --noPS < output.fasta > RNAfold.log
if [[ ! -d ./pp ]]; then mkdir pp; fi
if [[ ! -d ./ps ]]; then mkdir ps; fi
if [[ -e ./file_list.txt ]]; then rm file_list.txt ; fi
fileNum=0
for file in *dp.ps ; do 
	java -classpath ../ ps2pp $file > ./pp/${file%*.ps}.pp 
	(( fileNum += 1 ))
	mv ${file} ./ps
	echo `pwd`/pp/${file%*.ps}.pp >> file_list.txt
	>&2 printf "\r$fileNum higher PI structures processed into base pairing matrices" 
done
echo