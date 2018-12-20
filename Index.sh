#!/bin/bash

clear

SKIP=1

show_help()
{
echo "
        Usage: Index.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nCompute the average Adaptaion Index.\nThis script scans all the child directories for C1 fasta files, generates the FUC matrix for all the Ci genes and computes the average adaptation index I.\nGenerates a I vs Evolutionary distance plot.\n\n'
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

    if [[ ! -f rgf.txt ]] ; then
        echo 'File "rgf.txt" is not there, skipping.'
        cd ..
        continue
    fi

    #echo $(pwd) #check correct directory
    touch sdata.txt
    #for fa in *.fa not correct order
    for fa in $(ls *.fa | sort -V); do base=$(basename "$fa" '.fa'); sed -e "s/^>.*/>${base}/" -e 's/_[Rr][Ee][Ff]$//' -e 's/_[Nn][Ee][Ww]$//' -e 's/__/_/' "$fa" >> sdata.txt; done
    grep '^>' sdata.txt | sed -e 's/^>//' > cdata.txt
    sed -e 's/\(^C.*\)_\{0,1\}woPHE/\1_woPHE/' -e 's/__/_/' cdata.txt > cdata.txted
    compile_codon_freqs < sdata.txt > fdata.txt
    sed -e 's/,/ /g' -e 's/|/ /g' fdata.txt > fdata.txted
    paste -d' ' cdata.txted fdata.txted > fdata.txt
	rm sdata.txt cdata.txt cdata.txted fdata.txted
    cd ..

    Rscript --vanilla Index.R "./${F}"

done
