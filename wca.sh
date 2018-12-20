#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: wca.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nCompute Within-groups Corresppondence analysis.\nThis script scans all the child directories for C1 fasta files and computes the wca.\nThree graphics are generated. Fig1/2: C1 vs C2 - Genes. Fig3: C1 vs C2 - Codons.\n\n'

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
    #echo $(pwd) #check correct directory

    cat $(ls C1*.fa | grep -E '^C1_?[NRnr]?.?.?\.fa') > genes.txt
    cat $(ls *.modal_freq.fas | sort -V) > modal.txt
    grep '^>' genes.txt | sed -e 's/^>.*/genes/' > gnames.txt
    grep '^>' modal.txt | sed -e 's/^>.*/modal/' > mnames.txt
    cat gnames.txt mnames.txt > tname.txt
    cat genes.txt modal.txt > wca.txt
    
    compile_codon_counts  < wca.txt > wca.counts
    sed -e 's/  / /g' wca.counts > wca.counted
    cut -f 1 wca.counted > wca.countedited

    grep -E '^>.' genes.txt | sed -e 's/>//' -e 's/\s\{1,3\}.*//' > rgname.txt
    grep -E '^>.' modal.txt | sed -e 's/^>\.\/.*\//>/' -e 's/_[Rr][Ee][Ff]$//' -e 's/_[Nn][Ee][Ww]$//' -e 's/_\{0,2\}woPHE/_woPHE/' -e 's/>//' -e 's/\s\{1,3\}.*//' > rmname.txt
    cat rgname.txt rmname.txt > rnames.txt
    
    paste -d' ' rnames.txt tname.txt wca.countedited > wca.counts 

    rm genes.txt mnames.txt gnames.txt modal.txt tname.txt wca.txt rgname.txt rmname.txt rnames.txt wca.counted wca.countedited

    cd ..
    Rscript --vanilla graph_wca.R "./${F}" wca.counts
    
done
