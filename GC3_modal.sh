#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: GC3_modal.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nCompute GC3 of modal sequences and generates a GC3 vs ancestry plot.\nThis script scans all the child directories for CC1 fasta files, generates the FUC matrix for all the CCi genes and computes the GC percent of the third codon base.\nGenerates a GC3 vs Evolutionary distance graphic.\n\n'
exit 1
}

while getopts yh option; do
    case "${option}" in
        y) SKIP=0;;
        h) show_help;;
    esac
done

if [ ${SKIP} != 1 ]
then
    echo -e 'Beginning calculations.\n'
else
    echo -e 'PRESS ANY KEY TO CONTINUE.\n'
    read -n 1 -s
fi


ls -d */ | while read F; do 
    echo Processing "${F}"
    cd "./${F}/"
    
    if [[ ! -f dist.txt ]] ; then
        echo 'File "dist.txt" is not there, skipping.'
        cd ..
        continue
    fi
    
   	#GC3
    touch sdata.txt
	for fa in $(ls *.modal_freq.fas | sort -V); do base=$(basename "$fa" '.modal_freq.fas'); sed -e "s/^>.*/>${base}/" -e 's/_[Rr][Ee][Ff]$//' -e 's/_[Nn][Ee][Ww]$//' -e 's/__/_/' "$fa" >> sdata.txt; done
    grep '^>' sdata.txt | sed -e 's/^>//' > cdata.txt
    sed -e 's/\(^C.*\)_\{0,1\}woPHE/\1_woPHE/' -e 's/__/_/' cdata.txt > cdata.txted
	compile_codon_counts < sdata.txt > fcmdata.txt
	sed -e 's/,/ /g' -e 's/  / /g' fcmdata.txt > fcmdata.txted
    cut -f 1-59 -d" " fcmdata.txted > fcmdata.txt
	paste -d' ' cdata.txted fcmdata.txt > fcmdata.txted
	mv fcmdata.txted fcmdata.txt
	rm sdata.txt cdata.txt cdata.txted 

    cd ..

    Rscript --vanilla GC3_modal.R "./${F}"

done
