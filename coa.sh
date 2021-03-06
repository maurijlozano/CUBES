#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: coa.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nCompute Corresppondence analysis of RSCU.\nThis script scans all the child directories for CC1 fasta files and computes the COA using codonw software.\nThree graphics are generated. Fig1/2: C1 vs C2 - Genes. Fig3: C1 vs C2 - Codons.\n\n'

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
    cat $(ls C1*.fa | grep -E '^C1_?[NRnr]?.?.?\.fa' | grep -E '^C1[0-9]\.fa' -v) $(ls *.modal_freq.fas | sort -V) > coa.txt #-V sorts numerical values in the correct order 1 -> 10
    sed -e 's/^>\.\/.*\//>/' -e 's/_[Rr][Ee][Ff]$//' -e 's/_[Nn][Ee][Ww]$//'  -e 's/_\{0,2\}woPHE/_woPHE/' coa.txt > coa.txted
    codonw coa.txted -nomenu -nowarn -silent -coa_rscu -noblk
    rm hilo.coa fop.coa cusort.coa coa.out coa.blk cbi.coa cai.coa eigen.coa coa_raw coa.txt coa.txted
    sed -e 's/\(^[0-9]*\)_/\1 /' -e 's/\(C[0-9]*\)_*woPHE/\1woPHE/' -e 's/[a-zA-Z0-9\|]*__*[a-zA-Z0-9\.]*_*[a-zA-Z0-9]*_*[a-zA-Z0-9]*/genes/' -e 's/label/label class/' -e 's/ _/ /' genes.coa | cut -f 1-4 -d' ' > genes_C.coa
	rm genes.coa
	mv codon.coa codon_C.coa
	mv summary.coa summary_C.coa
	touch coa_var_axis1_2_C.txt
	echo 'varC1 varC2' > coa_var_axis1_2_C.txt
	sed '17q;d' summary_C.coa | sed -e 's/   / /g' | cut -f3,7 -d' ' >> coa_var_axis1_2_C.txt
	sed -e 's/\+//g' coa_var_axis1_2_C.txt > coa_var_axis1_2_C.txted
	mv coa_var_axis1_2_C.txted coa_var_axis1_2_C.txt
    cd ..
    Rscript --vanilla graph_coa.R "./${F}" genes_C.coa
    Rscript --vanilla graph_codons.R "./${F}" codon_C.coa
done
