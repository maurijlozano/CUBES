#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: Index_modal.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"


echo -e '\nCompute the average Adaptaion Index of Modal frequencies.\nThis script scans all the child directories for C1*.modal_freq files and computes de adaptation index I.\nGenerates a I vs Evolutionary distance plot.\n\n'

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
    
    cat $(ls *.modal_freq | sort -V) | sed -e 's/,/ /g' -e 's/\t/ /g' -e 's/|/ /g' > C_modal_freqs.txt
    ls *.modal_freq | sort -V | sed -e 's/_\{0,2\}\(woPHE\).*/:\1/' -e 's/_.*//' -e 's/\..*//' -e 's/:/_/' > Cindex.txt
    paste -d' ' Cindex.txt C_modal_freqs.txt > fmdata.txt
    rm C_modal_freqs.txt Cindex.txt

	

    cd ..

    Rscript --vanilla Index_modal.R "./${F}"

done
