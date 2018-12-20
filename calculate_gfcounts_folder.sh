#!/bin/bash
clear

SKIP=1
show_help()
{
echo "
        Usage: calculate_gfcounts_folder.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help

"
echo -e '\nCompute counts from Genome Fasta files.\nThis script scans the directorie for GNM.fa file and computes counts using G. Olsen software.\n\n'

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
    echo Processing "GENOME fasta file"
    compile_codon_counts  < GNM.fa > GNM.count
    sed -e 's/  / /g' GNM.count > genome.freq.counts
    cut -f 1 genome.freq.counts > genome.freq.count
    cut -f 2 genome.freq.counts > genome.freq.rname
    paste -d' ' genome.freq.rname genome.freq.count > genome.freq.counts
    rm *.count genome.freq.rname
	cd ..
done
