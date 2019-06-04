#!/bin/bash

clear

SKIP=1

show_help()
{
echo "
        Usage: GC3.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"
echo -e '\nCompute the average GC3 and generates a GC3 vs ancestry plot.\nThis script scans all the child directories for CC1 fasta files, generates the FUC matrix for all the CCi genes and computes the GC percent of the third codon base.\nGenerates a GC3 vs Evolutionary distance graphic.\n\n'
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


    #echo $(pwd) #check correct directory
    touch sdata.txt
    #for fa in *.fa not correct order
    for fa in $(ls *.fa | grep 'GNM' -v |sort -V); do base=$(basename "$fa" '.fa'); sed -e "s/^>.*/>${base}/" -e 's/_[Rr][Ee][Ff]$//' -e 's/_[Nn][Ee][Ww]$//' -e 's/__/_/' "$fa" >> sdata.txt; done
	compile_codon_counts < sdata.txt > fcdata.txt
	sed -e 's/,/ /g' -e 's/  / /g' fcdata.txt > fcdata.txted
	cut -f 2 fcdata.txted > cdata.txt
    sed -e 's/\(^C.*\)_\{0,1\}woPHE/\1_woPHE/' -e 's/__/_/' cdata.txt > cdata.txted
    cut -f 1-59 -d" " fcdata.txted > fcdata.txt
	paste -d' ' cdata.txted fcdata.txt > fcdata.txted
	mv fcdata.txted fcdata.txt
	rm sdata.txt cdata.txt cdata.txted
    cd ..

    Rscript --vanilla GC3.R "./${F}"

done
