#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: coa_GNM_rv.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nCompute Corresppondence analysis of RSCU.\nThis script scans all the child directories for CC1 and singleton fasta files and computes the COA using codonw software.\nThree graphics are generated. Fig1/2: C1 vs C2 - Genes. Fig3: C1 vs C2 - Codons.\n\n'

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
    SINGLE=$(ls single*.modal_freq.fas)
    sed -e 's/>.*/>Single/' ${SINGLE} > ${SINGLE}ed
    mv ${SINGLE}ed ${SINGLE}
    cat GNM.fa $(ls *.modal_freq.fas | sort -V) > coa_GNM.txt #-V sorts numerical values in the correct order 1 -> 10
    sed -e 's/^>\.\/.*\//>/' -e 's/_[Rr][Ee][Ff]$//' -e 's/_[Nn][Ee][Ww]$//'  -e 's/_\{0,2\}woPHE/_woPHE/' -e 's/\s/_/g' coa_GNM.txt > coa_GNM.txted
    mv coa_GNM.txted coa_GNM.txt
	if [[ ! -f coa_GNM.rscu ]]
	then
		cd ..
		./seq2rscu.pl "./${F}/"coa_GNM.txt "./${F}/"coa_GNM.rscu
		cd "./${F}/"
	fi

	rm coa_GNM.txt coa_GNM.counts coa_GNM.freqs &>/dev/null
    cd ..
    Rscript --vanilla graph_coaRSCU_GNM.R "./${F}"
done
