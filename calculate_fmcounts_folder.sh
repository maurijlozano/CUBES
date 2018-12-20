#!/bin/bash
clear

SKIP=1
show_help()
{
echo "
        Usage: calculate_fmcounts_folder.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help

"

echo -e '\nCompute counts from Fasta files.\nThis script scans the directorie for sequences representing modal codon usage and computes counts using G. Olsen software.\nThen, using codonW a table of Effective Number of Codons (ENc), G+C content of gene (all 3 codon positions) and GC of synonymous codons 3rd positions is generated.\n\n'

exit 1
}

while getopts ":yh" option; do
    case "${option}" in
        y) SKIP=0;;
        h) show_help
            exit 1
            ;;
       esac
done


if [[ $# -eq 0 ]] ; then
    echo 'Parameters required, run -h for help'
    exit 1
fi


if [ ${SKIP} != 1 ]
then
    echo -e 'Beginning calculations.\n'
else
    echo -e 'PRESS ANY KEY TO CONTINUE.\n'
    read -n 1 -s
fi

#generates the codon counts. 

ls -d */ | while read D; do
    echo Processing "${D}"
    cd "./${D}/"
    cat $(ls *.modal_freq.fas | sort -V) > modal.counts
    compile_codon_counts  < modal.counts > modal.count
    sed -e 's/  / /g' modal.count > modal.counts
    cut -f 1 modal.counts > modal.count
    echo -e "$(ls *.modal_freq.fas | sort -V)" | sed -e 's/.modal_freq.fas//' -e 's/[_ ][RrNn][Ee][WwFf]$//' > modal.rname
    paste -d' ' modal.rname modal.count | sed -e 's/\^\n//' > modal.counts
    rm *.count modal.rname
    cd ..
done
